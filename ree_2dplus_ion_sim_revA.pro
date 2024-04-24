; 2D+ ION SIMULATION CORE
; NOT FOR GENERAL USE

pro ree_2dplus_ion_sim, dt, it, jt, kt, npart, imax, nnew, nh, box_str, $
  B_str, Inj_str, x, y, z, px, py, pz, $
  Ex, Ey, Ez, Ey0, Epll, E1, E2, dsE, nEtimes, Eplevel, $
  Bx, By, Bz, Bmag, dBx, dBy, dBz, dB1, dB2, dB3, dsB, nBtimes, Bplevel, $
  pxn, pyn, pzn, pxold, pyold, pzold, gam, rat, plot_step, plot_mod, $
  Xpwr, Ypwr, Zpwr, psym, vth, Pdx, ind=ind, yplotrange = yplotrange, radiation=radiation, $ 
  maskzero=maskzero, zerocount=zerocount, mz3=mz3, mz4=mz4, maskleave = maskleave, wpi=wpi, $
  import_index = import_index, trajectories=trajectories, highp = highp, maskhot = maskhot, masklowp = masklowp, vinicial=vinicial, track = track, temps=temps
  
; NATURAL NUMBERS
ee    = 1.6d-19
me    = 0.911d-30
twopi = 2D*!dpi
mu0   = 4d-7*!dpi
mi    = 1.67d-27
Re    = 6374.d3
c     = 3d8
c2    = 9d16
sqrt2 = sqrt(2D)
W0    = mi*c2/ee ; REST ENERGY IN eV

; CONSTANTS FOR SIMULATION
qmdt   = ee/mi*dt/1000.0 ; E in mV/m
bqmdt  = ee/mi*dt ; FOR dB

flag  = 0 ; to exit a loop
highleave = 0
lowbleave = 0
lowbz_highpleave = 0
; UNPACK BOX STRUCTURE
xmin     = box_str.xmin
xmax     = box_str.xmax
xexitmin = box_str.xexitmin
xexitmax = box_str.xexitmax
xturbmin = box_str.xturbmin
xturbmax = box_str.xturbmax

zmin     = box_str.zmin
zmax     = box_str.zmax
zexitmin = box_str.zexitmin
zexitmax = box_str.zexitmax
zturbmin = box_str.zturbmin
zturbmax = box_str.zturbmax

ymin     = box_str.ymin
ymax     = box_str.ymax
yexitmin = box_str.yexitmin
yexitmax = box_str.yexitmax
yturbmin = box_str.yturbmin
yturbmax = box_str.yturbmax

; UNPACK B_str ### WILL REMOVE WITH NEW BCALC ROUTINE #####
Bx0 = B_str.Bx0
By0 = B_str.By0
Bz0 = B_str.Bz0
X0  = B_str.X0
Y0  = B_str.Y0
Z0  = B_str.Z0

;MAKE ARRAY to record x coordinates of test particles
;xtraj1 = MAKE_ARRAY(imax+1)
;ytraj1 = MAKE_ARRAY(imax+1)
;ztraj1 = MAKE_ARRAY(imax+1)

;MAKKE ARRAY TO RECORD IMPORTANT PARTICLES
steps = 100
if keyword_set(import_index) then begin
  restore, 'import_indext.sav'
  trajectories = temporary(make_array(imax/30 + 1, 12, n_elements(import_index)))
  ;vinicial = make_array(3, n_elements(import_index))
  track = 1
 endif
 
 ;;;;; INITIALIZE THE FLUX STRUCTURE TO RECORD THE FLUX OF THE LEAVING PARTICLES
;flux = {xminus:make_array(imax+1), xplus:make_array(imax+1), yminus:make_array(imax+1),$
 ; yplus:make_array(imax+1),zminus:make_array(imax+1), zplus:make_array(imax+1)}
  
;nleave = nnew*imax + npart 
 nleave = npart 
leave_x   = dblarr(nleave)
leave_y   = dblarr(nleave)
leave_z   = dblarr(nleave) 
leave_px  = dblarr(nleave)
leave_py  = dblarr(nleave)
leave_pz  = dblarr(nleave)

;maskhot   = boolarr(npart)

maskleave = boolarr(npart)
;maskLowBz = boolarr(npart)
maskzero = boolarr(npart)
;masklowp = boolarr(npart)

mz2 = boolarr(npart)
mz3 = boolarr(npart)
mz4 = boolarr(npart)
zerocount = intarr(npart)
highp = boolarr(npart)

ileave = 0l
iz = 0
; MAIN LOOP

FOR i = 0d, imax do BEGIN & $

  ; REMOVE OUT OF BOUNDARY PARTICLES. GAM=0 INDICATES NO PARTICLE
  indbad = where( (x[ind] LT xexitmin) OR (x[ind] GT xexitmax) OR $
                  (y[ind] LT yexitmin) OR (y[ind] GT yexitmax) OR $
                  (z[ind] LT zexitmin) OR (z[ind] GT zexitmax), nindbad)
  
  if nindbad GT 0 then begin
    leave_x[ind[indbad]]   =  x[ind[indbad]]
    leave_y[ind[indbad]]   =  y[ind[indbad]]
    leave_z[ind[indbad]]   =  z[ind[indbad]]
    leave_px[ind[indbad]]   = px[ind[indbad]]
    leave_py[ind[indbad]]   = py[ind[indbad]]
    leave_pz[ind[indbad]]   = pz[ind[indbad]]
    maskleave[ind[indbad]]   = 1
 
 ;ileave : ileave + nindbad - 1
 
    ileave = ileave + nindbad
    gam[ind[indbad]] = 0D
  endif
  ;RECORD HOW MANY PARTICLE CROSS THE PLANE OF NEUTRALITY
  ninzero = 0
  maskzero[ind[where((abs(z[ind]) LT 50*wpi), ninzero)]] = 1
  if ninzero GT 0 then begin
    zerocount[ind[where((abs(z[ind]) LT 2*wpi), ninzero2)]] += 1
   ;print,abs(z[ind[where(abs(z[ind]) LT 2*wpi)]])
   endif
  if i GT imax/4. then mz3[ind[where((abs(z[ind]) LT 2*wpi), ninzero3)]] = 1
        
  if i GT imax/2. then mz4[ind[where((abs(z[ind]) LT 2*wpi), ninzero4)]] = 1
   
      
      
  
   if i gt 20000 then begin
    nnew = 1
;    print,'yes'
    endif
  ;INJECT PARTICLES ADDS PARTICLES AT BOUNDARIES
  inject_particles3D, nnew, npart, i, jt=jt, x=x, y=y, z=z, px=px, py=py, pz=pz, $
    gam=gam, dsE=dsE, dsB=dsB, Inj_str=Inj_str, indnew = indnew
    
    ind = where(gam NE 0D, nind)  ;indbad was removed
  ;  print,'indnew =',py[indnew]

;RECORD TRAYECTORIES FOR LYAPUNOV EXPONENT CALCULATION
;

if track eq 0 and keyword_set(import_index) then begin
match, indnew, import_index,suba,subb
if n_elements(subb) gt 0 then begin
  vinicial[0, subb]  = px[import_index[subb]]
  vinicial[1, subb]  = py[import_index[subb]]
  vinicial[2, subb]  = pz[import_index[subb]]
 endif
  endif
   

if keyword_set(import_index) and track and i mod steps eq 0 then begin
  trajectories[uint(i/steps), 0, *]  = x[import_index] 
  trajectories[uint(i/steps), 1, *]  = y[import_index]
  trajectories[uint(i/steps), 2, *]  = z[import_index]
  trajectories[uint(i/steps), 3, *]  = px[import_index]
  trajectories[uint(i/steps), 4, *]  = py[import_index]
  trajectories[uint(i/steps), 5, *]  = pz[import_index]
  trajectories[uint(i/steps), 6, *]  = Ex[import_index]
  trajectories[uint(i/steps), 7, *]  = Ey[import_index]
  trajectories[uint(i/steps), 8, *]  = Ez[import_index]
  trajectories[uint(i/steps), 9, *]  = Bx[import_index] + dBx[import_index]
  trajectories[uint(i/steps), 10, *] = By[import_index] + dBy[import_index]
  trajectories[uint(i/steps), 11, *] = Bz[import_index] + dBz[import_index]
 ; if i lt 1 then begin 
 ;   HighMaxBoltz = where(trajectories[1, 3, *]^2+trajectories[1, 4, *]^2+trajectories[1, 5, *]^2 gt $
 ;     mean(trajectories[1, 3, *]^2+trajectories[1, 4, *]^2+trajectories[1, 5, *]^2 ))
 ;     maskhot[HighMaxBoltz] = 1
 ;     endif
;  print,z[import_index]
 ;print,ex[import_index]
endif
 

  
  ;Isolate the faster particles
  
  highp1 =  ind[where(abs(py[ind]) GT 2*mean(abs(py[ind])))]
 ;  lowp1 =  indnew[where(py[indnew]^2 + px[indnew]^2 + pz[indnew]^2 lT .67*mean(py[ind]^2 +px[ind]^2+pz[ind]^2),n)]
 ;  if n gt 0 then masklowp[lowp1] = 1
 ;  lowp1 =  indnew[where(py[indnew]^2 + px[indnew]^2 + pz[indnew]^2 gT  2*mean(py[ind]^2 +px[ind]^2+pz[ind]^2),n)]
 ;  if n gt 0 then maskhot[lowp1] = 1
 ;  
 ;  lowp1 =  indnew[where(py[indnew]^2 + px[indnew]^2 + pz[indnew]^2 lT .67*vth^2)]
 ;  if n gt 0 then masklowp[lowp1] = 1
 ;  lowp1 =  indnew[where(py[indnew]^2 + px[indnew]^2 + pz[indnew]^2 gT  2*vth^2)]
 ;  if n gt 0 then maskhot[lowp1] = 1
temps[i] = (mean(px[ind]^2*mi/ee)+mean(px[ind]^2*mi/ee)+mean(py[ind]^2*mi/ee))*.5
   
  ;LowBz =   ind[where(abs(Bz[ind]) LT 1e-11, nlowBz)] ;this doesn't remember what had lobB
  
  ;masklowbz[ind[[where(abs(Bz[ind]) LT 1e-11 and gam ne 0, nlowBz)]]] = 1
  ;highp[ind[where(py[ind] GT 3*mean(abs(py[ind])))]] = 1
  ;
  
  if indbad[0] ne -1 then begin
  match,ind[indbad],highp1,suba, subb
  ind2=ind[indbad]
  highp[ind2[suba]] = 1
  
 ; match,ind[indbad],where(masklowbz),suba, subb
 ; if (suba[0] ne -1) then begin
 ;   ; print, 'High velocity particles', ind[indbad[suba]], 'have left at i=', i
 ;   LowBleave+= n_elements(suba)
 ; endif
  
  ; highpleave_ind = ind2[suba]
 ;  match, highpleave_ind, LowBz,suba, subb ;see which particles that have had lowBw leave with highp
  ;   if keyword_set(suba) then begin
  ; if suba[0] ne -1 then lowbz_highpleave += n_elements(subb)
   endif
  ; endif
   

   
   

  
 
      
  ; MOVE PARTICLES
  x(ind)     = x(ind) + px(ind)*dt/gam(ind)
  y(ind)     = y(ind) + py(ind)*dt/gam(ind)
  z(ind)     = z(ind) + pz(ind)*dt/gam(ind)

  ; UPDATE B;
  bcalc, x, z, ind, B_str, Bx=Bx, Bz=Bz, Bmag=Bmag, $
    sechx=pxn, tanhx=pyn, sechz=pzn, atan_sinhz=rat		; DUMMY ARRAYS
  
  ; UPDATE RANDOM E
  Xpwr(ind) = -0.5D*cos( ( ((x(ind)+xturbmax)/Pdx) < (!dpi) ) > 0D ) - $
               0.5D*cos( ( ((xturbmax-x(ind))/Pdx) < (!dpi) ) > 0D ) & $
  Ypwr(ind) = -0.5D*cos( ( ((y(ind)+yturbmax)/Pdx) < (!dpi) ) > 0D ) - $
               0.5D*cos( ( ((yturbmax-y(ind))/Pdx) < (!dpi) ) > 0D ) & $
  Zpwr(ind) = -0.5D*cos( ( ((z(ind)+zturbmax)/Pdx) < (!dpi) ) > 0D ) - $
               0.5D*cos( ( ((zturbmax-z(ind))/Pdx) < (!dpi) ) > 0D ) & $
  rat(ind) = Xpwr(ind)*Ypwr(ind)*Zpwr(ind) ; NOTE RAT IS A DUMMY ARRAY
  Ex(ind)  = (Epll(dsE(ind)+it)*Bx(ind) - E2(dsE(ind)+it)*Bz(ind)) * $
              Eplevel*rat(ind)/Bmag(ind)
  Ey(ind)  = E1(dsE(ind)+it)*Eplevel*rat(ind) + Ey0
  Ez(ind)  = (Epll(dsE(ind)+it)*Bz(ind) + E2(dsE(ind)+it)*Bx(ind)) * $
              Eplevel*rat(ind)/Bmag(ind)
  indroll  = where( (dsE(ind) + it) GE (nEtimes-1), nindroll)
  if nindroll GT 0 then dsE(ind(indroll)) = dsE(ind(indroll))-nEtimes
  it = it + 1L

  ; UPDATE dB 
  dbx(ind)   = dB1(dsB(ind)+kt)*Bplevel*rat(ind) & $
  dby(ind)   = dB2(dsB(ind)+kt)*Bplevel*rat(ind) & $
  dbz(ind)   = dB3(dsB(ind)+kt)*Bplevel*rat(ind) & $
  indroll  = where( (dsB(ind) + it) GE (nBtimes-1), nindroll)
  if nindroll GT 0 then dsB(ind(indroll)) = dsB(ind(indroll))-nBtimes
  kt = kt + 1L

; print,ex[ind]
  ; ADVANCE: NOTE THAT pxn, pxold, rat ARE DUMMMY ARRAYS
  boris_advance, px, py, pz, gam, Bx+dBx, By+dBy, Bz+dbz, Ex, Ey, Ez, $
    ind, bqmdt, qmdt, pxn, pyn, pzn, pxold, pyold, pzold, rat

  ;UPDATE GAM
  gam(ind)  = sqrt(px(ind)*px(ind)+py(ind)*py(ind)+pz(ind)*pz(ind)+c2)/c & $
   

   ;OUTPUT ***** PLOT *****
  ;if ((i mod plot_step) EQ plot_mod) AND (plot_step GT 0) then $
   ; ree_run_time_plot, x, y, z, i, ind, nind, Box_Str
ENDFOR


radiation = make_array(n_elements(leave_x), 7)

radiation[*,0] = leave_x;[0:ileave -1]
radiation[*,1] = leave_y;[0:ileave -1]
radiation[*,2] = leave_z;[0:ileave -1]
radiation[*,3] = leave_px;[0:ileave -1]
radiation[*,4] = leave_py;[0:ileave -1]
radiation[*,5] = leave_pz;[0:ileave -1]


;print, 'Probability of a left particles being high velocity', float(highleave)/float(ileave)

;junk = where(px[ind] GT 2*mean(abs(px[ind])), nhigh)

;junk = where(gam ne 0, particlespresent)

;print, 'Frequency of high velocity particles', float(nhigh)/particlespresent

;print, 'Probabiliy of a left particle being LowBz', float(lowbleave)/ileave

;junk = where(gam ne 0 and maskLowBz, nLowBzpresent)

;print, 'Frequency of LowBz particles', float(nLowBzpresent)/float(particlespresent)


RETURN
end
