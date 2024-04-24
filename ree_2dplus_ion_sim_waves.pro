; %%%%%%%%%%%%%%%%%%%%%%%%% B GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%

; FUNCTION TO CALCULATE B SPEEDS FOR JULY 26 2017 EVENT ONLY
; SEE SEGA PAPER
; FIT FORMULA: 436D + alog(F/0.0221D)*215D

function b_speed_july26, Freqs
Vs = (400D3 + alog(Freqs/0.0221D)*215D3) > 400D3
return,Vs
end

; PRO TO GENERATE 3D B WITH CORRELATION LOSS
pro db_set, fmin, fmax, fbreak, nfreqs, alpha1, alpha2, Vmin, Vmax, CorPer, $
    EplRat, Brms, dBgStr=dBgStr, dAs=dAs, dFreqs=dFreqs, useformula=useformula

    seeds = [756,5675,55672,908,10890,13894,423434595,6900-966,61231239,86789,10080, 99123249, 7346484]


twopi   = 2D*!dpi
if not keyword_set(dAs) then dAs = 0.2D
if not keyword_set(dFreqs) then dFreqs = 0.2D

; SET FREQUENCIES
beta   = (fmax/fmin)^(1D/double(nfreqs))
freqr  = fmin*beta^dindgen(nfreqs) 
freqs  = freqr * (1D + (beta-1D)*(dFreqs/2D - dFreqs*randomu(seeds[0], nfreqs, /double)))

; VELOCITIES
if keyword_set(useformula) then $
  Vs = b_speed_july26(Freqs) * (0.7D + 0.6D*randomn(seeds[1], nfreqs, /double)) else $
  Vs     = Vmin + (Vmax-Vmin)*randomu(seeds[2], nfreqs, /double)

; SET RANDOM PHASES, VELOCITIES AND FREQUENCIES
phis   = randomu(seeds[3], nfreqs, /double)*twopi
tcorrs = CorPer*(0.94D + 0.12D*randomn(seeds[4], nfreqs, /double))/freqs 
dfs    = freqr*(beta-1D)
Ks     = freqs*twopi/Vs

; SET AMPLITUDES
As     = freqr^(-alpha1/2D)*sqrt(dfs)
ind    = where(freqr GE fbreak)
As(ind)= As(ind(0)) *  ( freqr(ind)/freqr(ind(0)) ) ^ (-alpha2/2D) * $
         sqrt(dfs(ind)/dfs(ind(0)))
As = As * (1D + dAs*randomu(seeds[5], nfreqs, /double))
Norm   = sqrt(total(As*As))/sqrt(2D)

; E DIRECTION ELECTROMAGNETIC
dirx = (2D*randomu(seeds[6], nfreqs, /double) - 1D)*EplRat
diry = 2D*randomu(seeds[7], nfreqs, /double) - 1D
dirz = 2D*randomu(seeds[8], nfreqs, /double) - 1D
damp  = sqrt(dirx*dirx + diry*diry + dirz*dirz)
Exs = dirx/damp
Eys = diry/damp
Ezs = dirz/damp

; RANDOM CROSS DIRECTION
crsx = 2D*randomu(seeds[9], nfreqs, /double) - 1D
crsy = 2D*randomu(seeds[10], nfreqs, /double) - 1D
crsz = 2D*randomu(seeds[11], nfreqs, /double) - 1D
camp  = sqrt(crsx*crsx + crsy*crsy + crsz*crsz)

; SET B PERPENDICULAR TO E
Bdirx = (Eys*crsz - Ezs*crsy)
Bdiry = (Ezs*crsx - Exs*crsz)
Bdirz = (Exs*crsy - Eys*crsx)
Bamp = sqrt(Bdirx*Bdirx + Bdiry*Bdiry + Bdirz*Bdirz)
Bxs = Bdirx/Bamp
Bys = Bdiry/Bamp
Bzs = Bdirz/Bamp

; SET UP K
Kxs = (Bys*Ezs - Bzs*Eys)*Ks
Kys = (Bzs*Exs - Bxs*Ezs)*Ks
Kzs = (Bxs*Eys - Bys*Exs)*Ks

; MAKE B STRUCTURE
ws  = freqs*twopi
dBgStr = {Bxs:Bxs, Bys:Bys, Bzs:Bzs, ws:ws, Ks:Ks, Kxs:Kxs, Kys:Kys, Kzs:Kzs, $
          As:As, Norm:Norm, dfs:dfs, Vs:Vs, phis:phis, tcorrs:tcorrs, $
          Exs:Exs, Eys:Eys, Ezs:Ezs, Brms:Brms} 

; CHECK
;print, (Bxs*Kxs + Bys*Kys + Bzs*Kzs)/Ks
;print, (Exs*Kxs + Eys*Kys + Ezs*Kzs)/Ks
;print, (Exs*Bxs + Eys*Bys + Ezs*Bzs)

RETURN
end



; ************** 
; PRO TO GENERATE B AND EM PART OF E FOR SIMULATION RUN
; T IS A SINGLE VALUE!
; X, Y, Z 

; *****
pro bcalc_wave, t, dt, x, y, z, ind, nind, dBgStr, $
  dBx=dBx, dBy=dBy, dBz=dBz, dEx=dEx, dEy=dEy, dEz=dEz, i=i

gain = dBgStr.Brms/dBgStr.Norm 
FOR j=0L, nind-1 do BEGIN
  cargs   = cos(dBgStr.Kxs*x(ind(j)) + dBgStr.Kys*y(ind(j)) + $
                dBgStr.Kzs*z(ind(j)) - dBgStr.ws*t + dBgStr.phis)
  dbx(ind(j)) = total(dBgStr.As*dBgStr.Bxs*cargs)*gain
  dby(ind(j)) = total(dBgStr.As*dBgStr.Bys*cargs)*gain
  dbz(ind(j)) = total(dBgStr.As*dBgStr.Bzs*cargs)*gain
  dex(ind(j)) = total(dBgStr.As*dBgStr.Exs*dBgStr.Vs * cargs)*gain*1d-6
  dey(ind(j)) = total(dBgStr.As*dBgStr.Eys*dBgStr.Vs * cargs)*gain*1d-6
  dez(ind(j)) = total(dBgStr.As*dBgStr.Ezs*dBgStr.Vs * cargs)*gain*1d-6
ENDFOR
seeds=[i]
; DE-CORRELATION
dBgStr.phis = dBgStr.phis + (6D*dt/dBgStr.tcorrs)*randomn(seeds[0], /double)

RETURN
end

; ************** 
; PRO TO TEST B AND EM BY GENERATING TIME SERIES
; nt is number of points
pro btest_time_series, nt, dt, dBgStr, kern=kern, $
  t=t, dBx=dBx, dBy=dBy, dBz=dBz, dEx=dEx, dEy=dEy, dEz=dEz

seeds=[23674]
; T MUST BE AN ARRAY
if not keyword_set(kern) then kern = [0.18D, 0.64D, 0.18D] 
t = dindgen(nt)*dt
dBxr = dblarr(nt)
dByr = dblarr(nt)
dBzr = dblarr(nt)
dExr = dblarr(nt)
dEyr = dblarr(nt)
dEzr = dblarr(nt)

gain = dBgStr.Brms/dBgStr.Norm
FOR j=0L, nt-1 do BEGIN & $
  args = - dBgStr.ws*t(j) + dBgStr.phis & $
  dBxr(j) = total(dBgStr.As*dBgStr.Bxs*cos(args))*gain & $
  dByr(j) = total(dBgStr.As*dBgStr.Bys*cos(args))*gain & $
  dBzr(j) = total(dBgStr.As*dBgStr.Bzs*cos(args))*gain & $
  dExr(j) = total(dBgStr.As*dBgStr.Exs*dBgStr.Vs*cos(args))*gain*1d-6 & $
  dEyr(j) = total(dBgStr.As*dBgStr.Eys*dBgStr.Vs*cos(args))*gain*1d-6 & $
  dEzr(j) = total(dBgStr.As*dBgStr.Ezs*dBgStr.Vs*cos(args))*gain*1d-6 & $
  dBgStr.phis = dBgStr.phis + (6D*dt/dBgStr.tcorrs)*randomn(seed, /double)
ENDFOR

; FILTER
dBx = convol(dBxr, kern)
dBy = convol(dByr, kern)
dBz = convol(dBzr, kern)
dEx = convol(dExr, kern)
dEy = convol(dEyr, kern)
dEz = convol(dEzr, kern)

RETURN
end




; ************** 
; PRO TO MAKE A SINGLE SPECTRA OUT OF A TIME SERIES
; nt is number of points
pro btest_spectra, npts, dt, Din, beta, $
             freq=freq, spec=spec, flog=flog, slog=slog

; SET UP
slide = 0.5
nspec = floor(n_elements(Din)/(npts*slide))-1
result = fltarr(npts/2)

; FFT AND AVERAGE
FOR jj=0L, nspec-1 do BEGIN & $
  ree_fft_power, Din(jj*npts*slide:jj*npts*slide+npts-1), $
    npts, dt, spec=spec, freq=freq & $
  result(*) = result(*) + spec(1:npts/2) & $
ENDFOR
freq = freq(1:*)
spec = result/nspec

; LOG SPECTRA
ree_v06_make_log_space, freq, spec, beta, Xlog=flog, Ylog=slog

RETURN
end



; %%%%%%%%%%%%%%%%%%%%%%%%% E GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%
; FUNCTION TO CALCULATE B SPEEDS FOR JULY 26 2017 EVENT ONLY
; SEE SEGA PAPER
; B FORMULA: 436D + alog(F/0.0221D)*215D

function e_speed_july26, Freqs
Vs = 50D3 + alog(Freqs/0.005)*215D3
return,Vs
end


; ELECTROSTATIC
pro de_set, fmin, fmax, fbreak, nfreqs, alpha1, alpha2, Vmin, Vmax, CorPer, $
    EplRat, Erms, dEeStr=dEeStr, dAs=dAs, dFreqs=dFreqs, useformula=useformula
    
    seeds = [320810420110, 971392730719, 630435691399, 463021129392, 140412630836, 235640050319, 127415909340, 958471853222, 714794892117, 099838618922, 838042928698, 256923978833] 

twopi   = 2D*!dpi
if not keyword_set(dAs) then dAs = 0.2D
if not keyword_set(dFreqs) then dFreqs = 0.25D

; SET FREQUENCIES
beta   = (fmax/fmin)^(1D/double(nfreqs))
freqr  = fmin*beta^dindgen(nfreqs) 
freqs  = freqr * (1D + (beta-1D)*(dFreqs/2D - dFreqs*randomu(seeds[0], nfreqs, /double)))

; VELOCITIES
if keyword_set(useformula) then $
  Vs = e_speed_july26(Freqs) * (0.7D + 0.6D*randomn(seeds[1], nfreqs, /double)) else $
  Vs     = Vmin + (Vmax-Vmin)*randomu(seeds[2], nfreqs, /double)

; SET RANDOM PHASES, VELOCITIES AND FREQUENCIES
phis   = randomu(seeds[3], nfreqs, /double)*twopi
tcorrs = CorPer*(0.94D + 0.12D*randomn(seeds[4], nfreqs, /double))/freqs 
dfs    = freqr*(beta-1D)
Ks     = freqs*twopi/Vs

; E DIRECTION ELECTROSTATIC
dirx = (2D*randomu(seeds[5], nfreqs, /double) - 1D)*EplRat
diry = 2D*randomu(seeds[6], nfreqs, /double) - 1D
dirz = 2D*randomu(seeds[7], nfreqs, /double) - 1D
damp  = sqrt(dirx*dirx + diry*diry + dirz*dirz)
Exs = dirx/damp
Eys = diry/damp
Ezs = dirz/damp

; SET AMPLITUDES
As     = freqr^(-alpha1/2D)*sqrt(dfs)
ind    = where(freqr GE fbreak)
As(ind)= As(ind(0)) *  ( freqr(ind)/freqr(ind(0)) ) ^ (-alpha2/2D) * $
         sqrt(dfs(ind)/dfs(ind(0)))
As = As * (1D + dAs*randomu(seeds[8], nfreqs, /double))
Norm   = sqrt(total(As*As))/sqrt(2D)

; K ELECTROSTATIC - NEED TO INPUT dt TO REPLACE 1D-2
; DE-CORRELATE
Kxs = Ks * (Exs + randomn(seeds[9], nfreqs, /double)*1D-2/CorPer)
Kys = Ks * (Eys + randomn(seeds[10], nfreqs, /double)*1D-2/CorPer)
Kzs = Ks * (Ezs + randomn(seeds[11], nfreqs, /double)*1D-2/CorPer)
ws  = freqs*twopi

; MAKE E STRUCTURE
dEeStr = {Exs:Exs, Eys:Eys, Ezs:Ezs, ws:ws, Ks:Ks, Kxs:Kxs, Kys:Kys, Kzs:Kzs, $
          As:As, Norm:Norm, dfs:dfs, Vs:Vs, phis:phis, tcorrs:tcorrs, $
          Erms:Erms} 

RETURN
end

; ************** 
; PRO TO GENERATE ES PART OF E FOR SIMULATION RUN
; T IS A SINGLE VALUE!
; X, Y, Z 


; *****
pro ecalc_wave, t, dt, x, y, z, ind, nind, dEeStr,  $
  dEx=dEx, dEy=dEy, dEz=dEz, i=i

gain = dEeStr.Erms/dEeStr.Norm
FOR j=0L, nind-1 do BEGIN
  cargs   = cos(dEeStr.Kxs*x(ind(j)) + dEeStr.Kys*y(ind(j)) + $
                dEeStr.Kzs*z(ind(j)) - dEeStr.ws*t + dEeStr.phis)
  dex(ind(j)) = total(dEeStr.As*dEeStr.Exs*cargs)*gain
  dey(ind(j)) = total(dEeStr.As*dEeStr.Eys*cargs)*gain
  dez(ind(j)) = total(dEeStr.As*dEeStr.Ezs*cargs)*gain
ENDFOR

seeds=[i]
; DE-CORRELATION
dEeStr.phis = dEeStr.phis + (6D*dt/dEeStr.tcorrs)*randomn(seeds[0], /double)

RETURN
end


; ************** 
; PRO TO TEST ES BY GENERATING TIME SERIES
; nt is number of points
pro etest_time_series, nt, dt, dEeStr, $
  t=t, dEx=dEx, dEy=dEy, dEz=dEz, kern=kern

if not keyword_set(kern) then kern = [0.18D, 0.64D, 0.18D] 

; T MUST BE AN ARRAY
t = dindgen(nt)*dt
dExr = dblarr(nt)
dEyr = dblarr(nt)
dEzr = dblarr(nt)

gain = dEeStr.Erms/dEeStr.Norm

FOR j=0L, nt-1 do BEGIN & $
  args = - dEeStr.ws*t(j) + dEeStr.phis & $
  dExr(j) = total(dEeStr.As*dEeStr.Exs*cos(args))*gain & $
  dEyr(j) = total(dEeStr.As*dEeStr.Eys*cos(args))*gain & $
  dEzr(j) = total(dEeStr.As*dEeStr.Ezs*cos(args))*gain & $
  dEeStr.phis = dEeStr.phis + (6D*dt/dEeStr.tcorrs)*randomn(seed, /double)
ENDFOR

; FILTER
dEx = convol(dExr, kern)
dEy = convol(dEyr, kern)
dEz = convol(dEzr, kern)

RETURN
end


; %%%%%%%%%%%%%%%%%%%%%%%% EXIT/ENTER FLUX %%%%%%%%%%%%%%%%%%%%%%%%%
pro make_envel, box_str, di, nRE=nRE, envel=envel, ds=ds 

Re    = 6374.d3


RETURN
end



; %%%%%%%%%%%%%%%%%%%%%%%% CALCULATION OF TURBULENCE %%%%%%%%%%%%%%%%%%%%%%%%%
pro EMcalc, t, dt, x, y, z, gam, ind, nind, dBgStr, dEeStr, dEeStrP, $
  dumx, dumy, dumz, dBx=dBx, dBy=dBy, dBz=dBz, dEx=dEx, dEy=dEy, dEz=dEz, i=i

; DO EM WAVES
bcalc_wave, t, dt, x, y, z, ind, nind, dBgStr, $
  dBx=dBx, dBy=dBy, dBz=dBz, dEx=dEx, dEy=dEy, dEz=dEz, i=i
dBx(ind) = dBx(ind)*1d-9
dBy(ind) = dBy(ind)*1d-9
dBz(ind) = dBz(ind)*1d-9

; DO PERP ES WAVES
ecalc_wave, t, dt, x, y, z, ind, nind, dEeStr,  $
  dEx=dumx, dEy=dumy, dEz=dumz, i=i

dEx(ind) = dEx(ind) + dumx(ind)
dEy(ind) = dEy(ind) + dumy(ind)
dEz(ind) = dEz(ind) + dumz(ind)
  

; DO PAR ES WAVES
ecalc_wave, t, dt, x, y, z, ind, nind, dEeStrP,  $
  dEx=dumx, dEy=dumy, dEz=dumz, i=i

dEx(ind) = dEx(ind) + dumx(ind)
dEy(ind) = dEy(ind) + dumy(ind)
dEz(ind) = dEz(ind) + dumz(ind)
 
RETURN
end



; %%%%%%%%%%%%%%%%%%%%%%%%% INJECT PARTICLES 3D %%%%%%%%%%%%%%%%%%%%%%%%%%
; REMOVED DsE and DsB - NOT NEEDED

pro inject_particles3D_W, nnew, i, jt=jt, x=x, y=y, z=z, px=px, py=py, pz=pz, $
  gam=gam, Inj_str=Inj_str, npart=npart

  ; GENERATE INDOPEN
;indopen = where(gam EQ 0)
;indnew = indopen(0:nnew-1)

indices = lindgen(npart)
indnew  = indices[i:nnew + i -1]

   
Ijmax = jt+nnew-1

; ROLL OVER IF jt EXCEEDS ARRAY SIZE
IF (Ijmax GE Inj_str.ninjarr) then BEGIN
  jt = 0L
  Ijmax = jt+nnew-1
ENDIF

x(indnew)   = Inj_str.x(jt:Ijmax)
y(indnew)   = Inj_str.y(jt:Ijmax)
z(indnew)   = Inj_str.z(jt:Ijmax)
px(indnew)  = Inj_str.px(jt:Ijmax)
py(indnew)  = Inj_str.py(jt:Ijmax)
pz(indnew)  = Inj_str.pz(jt:Ijmax)
gam(indnew) = Inj_str.gam(jt:Ijmax)

jt = jt + nnew

RETURN
end




; %%%%%%%%%%%%%%%%%%%%%%%% MAIN PRO %%%%%%%%%%%%%%%%%%%%%%%%%
; 2D+ ION SIMULATION CORE
; NOT FOR GENERAL USE

pro ree_2dplus_ion_sim_waves, dt, it, jt, imax, nnew, box_str, $
  B_str, Inj_str, x, y, z, px, py, pz, $
  Ex, Ey, Ez, Ey0, dEeStr, dEeStrP, Eplevel, $
  Bx, By, Bz, Bmag, dBx, dBy, dBz, dBgStr, Bplevel, $
  pxn, pyn, pzn, pxold, pyold, pzold, gam, rat, plot_step, plot_mod, $
  Xpwr, Ypwr, Zpwr, psym, vth, Pdx, ind=ind, yplotrange = yplotrange, $
  record_flux=record_flux, Flux_str=Flux_Str, $
  record_temp=record_temp, Temp_Str=Temp_Str, $
  record_dist=record_dist, Dist_Str=Dist_Str, $
  record_den=record_den, Den_Str=Den_Str, $
  record_track=record_track, Track_Str=Track_Str, $
  radiation=radiation, npart=npart, $
  maskzero=maskzero, zerocount=zerocount, mz3=mz3, mz4=mz4, maskleave = maskleave, wpi=wpi, $
  import_index = import_index, trajectories=trajectories, highp = highp, maskhot = maskhot, masklowp = masklowp, off=off

  
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

;MAKE ARRAY TO RECORD THE TRAJECTORY OF NORMAL PARTICLES
track=0
steps = 20
if keyword_set(import_index) then begin
  
  restore, 'import_indexwall1.sav'
if keyword_set(off) then restore, 'import_indexwall0.sav'
  ;  restore, 'import_index4r.sav'
  ;  restore, 'import_indexlonge0b0.sav' ;this is maskzero only
  ; restore, 'import_indexbobleave.sav'
  trajectories = temporary(make_array(imax/steps + 1 , 12, n_elements(import_index)))
  track = 1
endif

;;;;; INITIALIZE THE FLUX STRUCTURE TO RECORD THE FLUX OF THE LEAVING PARTICLES


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
FOR i = 0L, imax-1 do BEGIN & $

  ; REMOVE OUT OF BOUNDARY PARTICLES. GAM=0 INDICATES NO PARTICLE
  indexit = where( (x(ind) LT xexitmin) OR (x(ind) GT xexitmax) OR $
                   (y(ind) LT yexitmin) OR (y(ind) GT yexitmax) OR $
                   (z(ind) LT zexitmin) OR (z(ind) GT zexitmax), nindexit)
  
  if (nindexit GT 0) and keyword_set(record_flux) and i gt 1000 then $
    ree_record_flux, x, y, z, gam, ind(indexit), nindexit, W0, box_str, it, $
    Flux_str=Flux_Str

  leave_x[ind[indexit]]   =  x[ind[indexit]]
  leave_y[ind[indexit]]   =  y[ind[indexit]] ;the ys of the things that leave
  leave_z[ind[indexit]]   =  z[ind[indexit]]
  leave_px[ind[indexit]]   = px[ind[indexit]]
  leave_py[ind[indexit]]   = py[ind[indexit]]
  leave_pz[ind[indexit]]   = pz[ind[indexit]]
  maskleave[ind[indexit]]   = 1
  indleave = ind[indexit]
  ileave = ileave + nindexit

  if nindexit GT 0 then     gam[ind[indexit]] = 0D

  ;RECORD HOW MANY PARTICLE CROSS THE PLANE OF NEUTRALITY
  ninzero = 0
  maskzero[ind[where((abs(z[ind]) LT 50*wpi), ninzero)]] = 1
  if ninzero GT 0 then begin
    zerocount[ind[where((abs(z[ind]) LT 30*wpi), ninzero2)]] += 1
    ;print,abs(z[ind[where(abs(z[ind]) LT 2*wpi)]])
  endif
  if i GT imax/4. then mz3[ind[where((abs(z[ind]) LT 2*wpi), ninzero3)]] = 1
  if i GT imax/2. then mz4[ind[where((abs(z[ind]) LT 2*wpi), ninzero4)]] = 1



  ; INJECT PARTICLES ADDS PARTICLES AT BOUNDARIES
  ;RECORD TRAYECTORIES FOR LYAPUNOV EXPONENT CALCULATION
  ;


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
  endif



  ; INJECT PARTICLES ADDS PARTICLES AT BOUNDARIES
  inject_particles3D_W, nnew, i, jt=jt, x=x, y=y, z=z, px=px, py=py, pz=pz, $
    gam=gam, Inj_str=Inj_str, npart=npart
  ind = where(gam NE 0D, nind)
    
  ; RECORD DESNITY 
  if keyword_set(record_den) then $
    ree_record_density, x, y, z, ind, box_str, Den_Str=Den_Str
    
  ; MOVE PARTICLES
  x(ind)     = x(ind) + px(ind)*dt/gam(ind)
  y(ind)     = y(ind) + py(ind)*dt/gam(ind)
  z(ind)     = z(ind) + pz(ind)*dt/gam(ind)

  ; UPDATE B
  bcalc, x, z, ind, B_str, Bx=Bx, Bz=Bz, Bmag=Bmag, $
    sechx=pxn, tanhx=pyn, sechz=pzn, atan_sinhz=rat		; DUMMY ARRAYS
  
  ; UPDATE TURBULENCE
  IF (Eplevel GT 1d-3) then $ 
    EMcalc, it*dt, dt, x, y, z, gam, ind, nind, dBgStr, dEeStr, dEeStrP, $
      pxn, pyn, pzn, dBx=dBx, dBy=dBy, dBz=dBz, dEx=Ex, dEy=Ey, dEz=Ez, i=i
  
  ; POWER ENVELOPE
  Xpwr(ind) = -0.5D*cos( ( ((x(ind)+xturbmax)/Pdx) < (!dpi) ) > 0D ) - $
               0.5D*cos( ( ((xturbmax-x(ind))/Pdx) < (!dpi) ) > 0D ) & $
  Ypwr(ind) = -0.5D*cos( ( ((y(ind)+yturbmax)/Pdx) < (!dpi) ) > 0D ) - $
               0.5D*cos( ( ((yturbmax-y(ind))/Pdx) < (!dpi) ) > 0D ) & $
  Zpwr(ind) = -0.5D*cos( ( ((z(ind)+zturbmax)/Pdx) < (!dpi) ) > 0D ) - $
               0.5D*cos( ( ((zturbmax-z(ind))/Pdx) < (!dpi) ) > 0D ) & $
  rat(ind) = Xpwr(ind)*Ypwr(ind)*Zpwr(ind) ; NOTE RAT IS A DUMMY ARRAY

  Ex(ind)  = Ex(ind)*Eplevel*rat(ind)
  Ey(ind)  = Ey(ind)*Eplevel*rat(ind) + Ey0*(0.75D + 0.25D*Ypwr(ind))
  Ez(ind)  = Ez(ind)*Eplevel*rat(ind)
  it = it + 1L

  ; UPDATE dB 
  dbx(ind)   = dbx(ind)*Bplevel*rat(ind) & $
  dby(ind)   = dby(ind)*Bplevel*rat(ind) & $
  dbz(ind)   = dbz(ind)*Bplevel*rat(ind) & $

  ; ADVANCE: NOTE THAT pxn, pxold, rat ARE DUMMMY ARRAYS
  boris_advance, px, py, pz, gam, Bx+dBx, By+dBy, Bz+dbz, Ex, Ey, Ez, $
    ind, bqmdt, qmdt, pxn, pyn, pzn, pxold, pyold, pzold, rat

  ;UPDATE GAM
  gam(*) = 0D
  gam(ind)  = sqrt(px(ind)*px(ind)+py(ind)*py(ind)+pz(ind)*pz(ind)+c2)/c

  ; OUTPUT ***** PLOT *****
;  if ((i mod plot_step) EQ plot_mod) AND (plot_step GT 0) then $
 ;   ree_run_time_plot, x, y, z, i, ind, nind, Box_Str

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
        
  ; TRACK PARTICLES
  if keyword_set(record_track) then $
    if ((i mod Track_Str.Track_Step) EQ 0) then $
      ree_track_particles, x, y, z, gam, Bx+dBx, By+dBy, Bz+dbz, W0, $
        Track_Str=Track_Str
        
        ;Record position and velocity out of all faces

      radiation[indleave,0] = leave_x[indleave];[0:ileave -1]
      radiation[indleave,1] = leave_y[indleave];[0:ileave -1]
      radiation[indleave,2] = leave_z[indleave];[0:ileave -1]
      radiation[indleave,3] = leave_px[indleave];[0:ileave -1]
      radiation[indleave,4] = leave_py[indleave];[0:ileave -1]
      radiation[indleave,5] = leave_pz[indleave];[0:ileave -1]


ENDFOR

RETURN
end








