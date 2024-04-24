; 2D+ ION SIMULATION CORE
; NOT FOR GENERAL USE

pro ree_2dplus_ion_sim, dt, it, jt, kt, npart, imax, nnew, nh, box_str, $
  B_str, Inj_str, x, y, z, px, py, pz, $
  Ex, Ey, Ez, Ey0, Epll, E1, E2, dsE, nEtimes, Eplevel, $
  Bx, By, Bz, Bmag, dBx, dBy, dBz, dB1, dB2, dB3, dsB, nBtimes, Bplevel, $
  pxn, pyn, pzn, pxold, pyold, pzold, gam, rat, plot_step, plot_mod, $
  Xpwr, Ypwr, Zpwr, psym, vth, Pdx, ind=ind, yplotrange = yplotrange, $
  record_flux=record_flux, Flux_str=Flux_Str, $
  record_temp=record_temp, Temp_Str=Temp_Str, $
  record_dist=record_dist, Dist_Str=Dist_Str, $
  record_den=record_den, Den_Str=Den_Str, $
  record_track=record_track, Track_Str=Track_Str

  
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


; MAIN LOOP
FOR i = 0L, imax-1 do BEGIN & $

  ; REMOVE OUT OF BOUNDARY PARTICLES. GAM=0 INDICATES NO PARTICLE
  indexit = where( (x(ind) LT xexitmin) OR (x(ind) GT xexitmax) OR $
                   (y(ind) LT yexitmin) OR (y(ind) GT yexitmax) OR $
                   (z(ind) LT zexitmin) OR (z(ind) GT zexitmax), nindexit)
  if (nindexit GT 0) and keyword_set(record_flux) then $
    ree_record_flux, x, y, z, gam, ind(indexit), nindexit, W0, box_str, it, $
      Flux_str=Flux_Str
  if nindexit GT 0 then gam(ind(indexit)) = 0D


  ; INJECT PARTICLES ADDS PARTICLES AT BOUNDARIES
  inject_particles3D, nnew, jt=jt, x=x, y=y, z=z, px=px, py=py, pz=pz, $
    gam=gam, dsE=dsE, dsB=dsB, Inj_str=Inj_str
  ind = where(gam NE 0D, nind)
    
  ; RECORD DESNITY 
  if keyword_set(record_den) then $
    ree_record_density, x, y, z, ind, box_str, Den_Str=Den_Str
    
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

  ; ADVANCE: NOTE THAT pxn, pxold, rat ARE DUMMMY ARRAYS
  boris_advance, px, py, pz, gam, Bx+dBx, By+dBy, Bz+dbz, Ex, Ey, Ez, $
    ind, bqmdt, qmdt, pxn, pyn, pzn, pxold, pyold, pzold, rat

  ;UPDATE GAM
  gam(ind)  = sqrt(px(ind)*px(ind)+py(ind)*py(ind)+pz(ind)*pz(ind)+c2)/c & $

  ; OUTPUT ***** PLOT *****
  if ((i mod plot_step) EQ plot_mod) AND (plot_step GT 0) then $
    ree_run_time_plot, x, y, z, i, ind, nind, Box_Str

  ; RECORD TEMPERATURES
  if keyword_set(record_temp) then $
    if ((i mod Temp_Str.Temp_Step) EQ 0) then $
      ree_record_temps, dt, x, y, z, px, py, pz, Bx, By, Bz, ind, Box_str, $
        Temp_Str=Temp_Str

  ; RECORD DISTRIBUTIONS
  if keyword_set(record_dist) then $
    if ((i mod Dist_STR.Dist_Step) EQ 0) then $
      ree_record_dist, x, y, z, px, py, pz, gam, ind, W0, box_str, $
        Dist_Str=Dist_Str
        
  ; TRACK PARtICLES
  if keyword_set(record_track) then $
    if ((i mod Track_Str.Track_Step) EQ 0) then $
      ree_track_particles, x, y, z, gam, Bx+dBx, By+dBy, Bz+dbz, W0, $
        Track_Str=Track_Str

ENDFOR

RETURN
end
