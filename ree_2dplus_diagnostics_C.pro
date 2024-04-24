; %%%%%%%%%%%%%%%%%%%%%%% PARTICLE PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%
pro ree_run_time_plot, x, y, z, i, ind, nind, Box_Str 

Re    = 6374.d3
ree_sym, 1.0, /fill

plot, x(ind)/Re, z(ind)/Re, psym=3, xran=[Box_Str.xmin, Box_Str.xmax]/Re, $
    yran = [Box_Str.zmin,Box_Str.zmax]/Re, xtit = 'X (R!DE!N)', $
    ytit = 'Z (R!DE!N)'
oplot, x(100:101)/Re, z(100:101)/Re, psym=8, col=2
oplot, x(500:501)/Re, z(500:501)/Re, psym=8, col=4
oplot, x(900:901)/Re, z(900:901)/Re, psym=8, col=6 
  
wi, 1
plot, y(ind)/Re, z(ind)/Re, psym=3, xran=[Box_Str.ymin, Box_Str.ymax]/Re, $
      yran = [Box_Str.zmin,Box_Str.zmax]/Re, xtit = 'Y (R!DE!N)', $
      ytit = 'Z (R!DE!N)', /xstyle
oplot, y(100:101)/Re, z(100:101)/Re, psym=8, col=2
oplot, y(500:501)/Re, z(500:501)/Re, psym=8, col=4
oplot, y(900:901)/Re, z(900:901)/Re, psym=8, col=6
wi, 0

print, i, nind
RETURN
end


; %%%%%%%%%%%%%%%%%%%%%%% DENSIY %%%%%%%%%%%%%%%%%%%%%%%%%%
pro ree_density_basic, x, y, z, px, py, pz, ind, Box_str, turb=turb, $
  plot=plot, ddist=ddist, nx=nx, ny=ny, nz=nz, xax=xax, yax=yax, zax=zax

; CHECK KEYWORDS
if not keyword_set(ddist) then ddist=1d6
Re    = 6374.d3

; SET UP HISTOGRAMS - X
npts  = round((Box_Str.xmax - Box_Str.xmin)/ddist)
binx  = (Box_Str.xmax - Box_Str.xmin)/npts
xax   = dindgen(npts)*binx + Box_Str.xmin
nx    = histogram(x(ind), bin=binx, min=Box_Str.xmin, max=Box_Str.xmax-1D)
if keyword_set(plot) then plot, xax/Re, double(nx)/max(nx), $
  xtit='X: blue, Y: green, Z:red (R!DE!N)', ytit = 'Density (Normalized)'
if keyword_set(plot) then oplot, xax/Re, double(nx)/max(nx), col=2

; Y
npts  = round((Box_Str.ymax - Box_Str.ymin)/ddist)
biny  = (Box_Str.ymax - Box_Str.ymin)/npts
yax   = dindgen(npts)*biny + Box_Str.ymin
ny    = histogram(y(ind), bin=biny, min=Box_Str.ymin, max=Box_Str.ymax-1D)
if keyword_set(plot) then oplot, yax/Re, double(ny)/max(ny), col=4

; Z
npts  = round((Box_Str.zmax - Box_Str.zmin)/ddist)
binz  = (Box_Str.zmax - Box_Str.zmin)/npts
zax   = dindgen(npts)*binz + Box_Str.zmin
nz    = histogram(z(ind), bin=binz, min=Box_Str.zmin, max=Box_Str.zmax-1D)
if keyword_set(plot) then oplot, zax/Re, double(nz)/max(nz), col=6

RETURN
end

; %%%%%%%%%%%%%%%%%%%%%%% DISTiBUTION FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%
pro ree_distribution_basic, x, y, z, px, py, pz, ind, Box_str, turb=turb, plot=plot, $ 
  Pmax=Pmax, Bin=Bin, fpx=fpx, fpy=fpy, fpz=fpz, pax=pax, fave=fave, Wax=Wax, $
  fabspx=fabspx, fabspy=fabspy, fabspz=fabspz, pabs=pabs

; NATURAL NUMBERS
ee    = 1.6d-19
mi    = 1.67d-27
c     = 3d8
c2    = 9d16
W0    = mi*c2/ee ; REST ENERGY IN eV

; SET UP INDUSE
IF keyword_set(turb) then BEGIN & $
  indturb = where( $
    (x(ind) GT Box_str.xturbmin) AND (x(ind) LT Box_str.xturbmax) AND $
    (y(ind) GT Box_str.yturbmin) AND (y(ind) LT Box_str.yturbmax) AND $
    (z(ind) GT Box_str.zturbmin) AND (z(ind) LT Box_str.zturbmax) ) & $
  induse = ind(indturb) & $
ENDIF else induse = ind

; KEYWORDS
if not keyword_set(Pmax) then Pmax = 9.99999d7
if not keyword_set(Bin) then Bin = 1d5

; DO ABS FIRST
npts   = ceil(Pmax/Bin)
pabs   = dindgen(npts)*Bin + Bin/2D
Wax    = (sqrt(pabs*pabs + c2) - c)*mi*c/ee
fabspx = histogram(abs(px(induse)), bin=bin, min=0.0D, max = Pmax)
fabspy = histogram(abs(py(induse)), bin=bin, min=0.0D, max = Pmax)
fabspz = histogram(abs(pz(induse)), bin=bin, min=0.0D, max = Pmax)
fave   = float(fabspx+fabspy+fabspz)/3.0

; PLOT
if keyword_set(plot) then plot, Wax, float(fave)>0.5, /xlog, /ylog, $
  xtit='Energy (eV)', ytit = 'Counts'

; DO FULL DISTRIBUTIONS
npts   = ceil(2*Pmax/Bin) + 1
pax   = dindgen(npts)*Bin - Bin*(npts-1)/2
fpx = histogram((px(induse)), bin=bin, min=-Pmax, max = Pmax)
fpy = histogram((py(induse)), bin=bin, min=-Pmax, max = Pmax)
fpz = histogram((pz(induse)), bin=bin, min=-Pmax, max = Pmax)

RETURN
end


; %%%%%%%%%%%%%%%%%%%%%%% TEMPERATURE DIAGNOSTICS %%%%%%%%%%%%%%%%%%%%%%%%%%

pro ree_print_temps, x, y, z, px, py, pz, Bx, By, Bz, ind, Box_str

; NATURAL NUMBERS
ee    = 1.6d-19
mi    = 1.67d-27
c     = 3d8
c2    = 9d16

; SKIP
print, ' '

; AVERAGE OVER BOX
TallB  = average(sqrt(px(ind)*px(ind)+py(ind)*py(ind)+pz(ind)*pz(ind)+c2)-c) * $
                 mi*c/ee*2D/3D
Btot = sqrt(Bx(ind)*Bx(ind) + By(ind)*By(ind) + Bz(ind)*Bz(ind))
ppll = (px(ind)*Bx(ind) + py(ind)*By(ind) + pz(ind)*Bz(ind) ) / Btot 
TparB  = average(sqrt(ppll*ppll+c2)-c)*2D*mi*c/ee
TperpB = (3D*TallB - TparB)/2D

print, TparB/1000D, TperpB/1000D, TallB/1000D, format =  $
 '("  TBoxPar  = ", F9.3,"  TBoxPerp  = ", F9.3, "  TBoxAve  = ", F9.3, " keV")'

indd = where( $
  (x(ind) GT Box_str.xturbmin) AND (x(ind) LT Box_str.xturbmax) AND $
  (y(ind) GT Box_str.yturbmin) AND (y(ind) LT Box_str.yturbmax) AND $
  (z(ind) GT Box_str.zturbmin) AND (z(ind) LT Box_str.zturbmax) ) & $

indt = ind(indd)

TallT  = average(sqrt(px(indt)*px(indt)+py(indt)*py(indt)+ $
                      pz(indt)*pz(indt)+c2)-c) * mi*c/ee*2D/3D
Btot = sqrt(Bx(indt)*Bx(indt) + By(indt)*By(indt) + Bz(indt)*Bz(indt))
ppll = (px(indt)*Bx(indt) + py(indt)*By(indt) + pz(indt)*Bz(indt) ) / Btot 
TparT  = average(sqrt(ppll*ppll+c2)-c)*2D*mi*c/ee
TperpT = (3D*TallT - TparT)/2D

print, TparT/1000D, TperpT/1000D, TallT/1000D, format =  $
 '("  TTurbPar  =", F9.3,"  TTurbPerp  =", F9.3, "  TTurbAve  =", F9.3, " keV")'

; SKIP
print, ' '

RETURN
end

; %%%%%%%%%%%%%%%%%%%%%%% TEMPERATURE DIAGNOSTICS %%%%%%%%%%%%%%%%%%%%%%%%%%

pro ree_record_temps, dt, x, y, z, px, py, pz, Bx, By, Bz, ind, Box_str, $
  Temp_Str=Temp_Str

; NATURAL NUMBERS
ee    = 1.6d-19
mi    = 1.67d-27
c     = 3d8
c2    = 9d16
Re    = 6374.d3

;UNPACK STRUCTURE
Buff  = Temp_Str.Buff*Re
kTemp = Temp_Str.kTemp
if (kTemp GT 0L) then $
  Temp_Str.t(ktemp) = Temp_Str.t(ktemp-1) + Temp_Str.Temp_Step*dt else $
    Temp_Str.t(ktemp) = 0D

; AVERAGE OVER BOX
TallB  = average(sqrt(px(ind)*px(ind)+py(ind)*py(ind)+pz(ind)*pz(ind)+c2)-c) * $
                 mi*c/ee*2D/3D
Btot = sqrt(Bx(ind)*Bx(ind) + By(ind)*By(ind) + Bz(ind)*Bz(ind))
ppll = (px(ind)*Bx(ind) + py(ind)*By(ind) + pz(ind)*Bz(ind) ) / Btot 
TparB  = average(sqrt(ppll*ppll+c2)-c)*2D*mi*c/ee
TperpB = (3D*TallB - TparB)/2D

Temp_Str.TBoxPar(ktemp)  = TparB/1000D
Temp_Str.TBoxPerp(ktemp) = TperpB/1000D
Temp_Str.TBoxAve(ktemp)  = TallB/1000D

; TURBULENT REGION
indd = where( $
  (x(ind) GT (Box_str.xturbmin+Buff)) AND $
  (x(ind) LT (Box_str.xturbmax-Buff)) AND $
  (y(ind) GT (Box_str.yturbmin+Buff)) AND $
  (y(ind) LT (Box_str.yturbmax-Buff)) AND $
  (z(ind) GT (Box_str.zturbmin+Buff)) AND $
  (z(ind) LT (Box_str.zturbmax-Buff)) ) & $

indt = ind(indd)

TallT  = average(sqrt(px(indt)*px(indt)+py(indt)*py(indt)+ $
                      pz(indt)*pz(indt)+c2)-c) * mi*c/ee*2D/3D
Btot = sqrt(Bx(indt)*Bx(indt) + By(indt)*By(indt) + Bz(indt)*Bz(indt))
ppll = (px(indt)*Bx(indt) + py(indt)*By(indt) + pz(indt)*Bz(indt) ) / Btot 
TparT  = average(sqrt(ppll*ppll+c2)-c)*2D*mi*c/ee
TperpT = (3D*TallT - TparT)/2D

Temp_Str.TturbPar(ktemp)  = TparT/1000D
Temp_Str.TturbPerp(ktemp) = TperpT/1000D
Temp_Str.TturbAve(ktemp)  = TallT/1000D

Temp_Str.kTemp = Temp_Str.kTemp + 1L

RETURN
end


; %%%%%%%%%%%%%%%%%%%%%%% RECORD EXITING FLUX MEASURMENT %%%%%%%%%%%%%%%%%%%%%

pro ree_record_flux, x, y, z, gam, Iexit, nindexit, W0, box_str, it, $
  Flux_str=Flux_Str

; Iexit = ind(indexit)
dw = Flux_str.dW
Flux_str.it = it

; CALCULATE ENERGY
W = W0*(gam(Iexit) - 1D)
iW = floor(W/dW) 
induse = where( (iW GE 0) AND (iW LT Flux_Str.Nflux), nuse)
Flux_str.Nexit(Flux_str.it) = nuse
if (nuse LE 0) then RETURN
iuse = Iexit(induse)

; FIRST RECORD BINS - HARD CODED FOR NOW - NEED TO FIX
jW = floor(W(induse)/10d3) < 10L
FOR i = 0, nuse-1 do Flux_str.FlxVtime(Flux_str.it, jW(i)) = $
                     Flux_str.FlxVtime(Flux_str.it, jW(i)) + W(induse(i))
Flux_str.FlxVtime(Flux_str.it,11) = total(W(induse))

; START LOOP
FOR i = 0, nuse-1 do BEGIN
  j = iW(induse(i))
  ; CHECK LOCATION
  if (x(iuse(i)) LT box_str.xexitmin) then Flux_str.mx(j) = Flux_str.mx(j)+1
  if (x(iuse(i)) GT box_str.xexitmax) then Flux_str.px(j) = Flux_str.px(j)+1
  if (y(iuse(i)) LT box_str.yexitmin) then Flux_str.my(j) = Flux_str.my(j)+1
  if (y(iuse(i)) GT box_str.yexitmax) then Flux_str.py(j) = Flux_str.py(j)+1
  if (z(iuse(i)) LT box_str.zexitmin) then Flux_str.mz(j) = Flux_str.mz(j)+1
  if (z(iuse(i)) GT box_str.zexitmax) then Flux_str.pz(j) = Flux_str.pz(j)+1
ENDFOR

RETURN
end

; %%%%%%%%%%%%%%%%%%%%%%% RECORD DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%
pro ree_record_dist, x, y, z, px, py, pz, gam, ind, W0, box_str, $
  Dist_Str=Dist_Str

; NATURAL NUMBERS
ee    = 1.6d-19
mi    = 1.67d-27
c     = 3d8
c2    = 9d16
Re    = 6374.d3

;UNPACK STRUCTURE
Pmax  = Dist_Str.Pmax
Bin   = Dist_Str.Bin
Turb  = Dist_Str.Turb
Wmax  = Dist_Str.Wmax
Wbin  = Dist_Str.Wbin
BuffSize = Dist_Str.BuffSize*Re
YStart = Dist_Str.YStart*Re
Ystop  = Dist_Str.Ystop*Re

; SET UP INDUSE
IF keyword_set(turb) then BEGIN & $
  indturb = where( (y(ind) GE YStart) AND (y(ind) LE YStop) AND $
    (x(ind) GT (Box_str.xturbmin+BuffSize)) AND $
    (x(ind) LT (Box_str.xturbmax-BuffSize)) AND $
    (z(ind) GT (Box_str.zturbmin+BuffSize)) AND $
    (z(ind) LT (Box_str.zturbmax-BuffSize)) ) & $
  induse = ind(indturb) & $
ENDIF else induse = ind

; DO FULL DISTRIBUTIONS
Dist_Str.fpx(*) = Dist_Str.fpx(*) + $
  histogram((px(induse)), bin=bin, min=-Pmax, max = Pmax)
Dist_Str.fpy(*) = Dist_Str.fpy(*) + $
  histogram((py(induse)), bin=bin, min=-Pmax, max = Pmax)
Dist_Str.fpz(*) = Dist_Str.fpz(*) + $
  histogram((pz(induse)), bin=bin, min=-Pmax, max = Pmax)
Dist_Str.fw(*) = Dist_Str.fw(*) + $
  histogram((gam(induse)-1D)*W0, bin=Wbin, min=0D, max = Wmax)

RETURN
end

; %%%%%%%%%%%%%%%%%%%%%%% RECORD DENSITY %%%%%%%%%%%%%%%%%%%%%
; NOTE THAT PXN, PYN, AND PZN ARE DUMMY ARRAYS
pro ree_record_density, x, y, z, ind, box_str, Den_Str=Den_Str

; SET UP INDUSE
indden = where( (y(ind) GT Den_Str.ymin) AND (y(ind) LT Den_Str.ymax) )
induse = ind(indden)

; ROUND
Den_Str.DenMap(floor( (x(induse)-box_str.xmin)/Den_Str.dx ) + $
               floor( (z(induse)-box_str.zmin)/Den_Str.dx ) * Den_Str.nx )++

RETURN
end


; %%%%%%%%%%%%%%%%%%%%%%% TRACK INDIVIDUAL PARTICLES %%%%%%%%%%%%%%%%%%%%%
; CREATES A TIME SERIES OF PARTICLE P, POSITION, B, and E
; USE WITH VERY FEW PARTICLES!!!!!
pro ree_track_particles, x, y, z, gam, BtX, Bty, BTz, W0, $
                        Track_Str=Track_Str

; SET UP
itrk = Track_Str.itrk
nm1  = Track_Str.ntrack - 1L
Btrk = sqrt(Btx(0:nm1)*Btx(0:nm1) + Bty(0:nm1)*Bty(0:nm1) +Btz(0:nm1)*BTz(0:nm1))

; RECORD
Track_Str.xtrk(itrk, *) = x(0:nm1)
Track_Str.ytrk(itrk, *) = y(0:nm1)
Track_Str.ztrk(itrk, *) = z(0:nm1)
Track_Str.Wtrk(itrk, *) = (gam(0:nm1)-1D)*W0
Track_Str.Btrk(itrk, *) = Btrk
Track_Str.itrk++

RETURN
end



; %%%%%%%%%%%%%%%%%%%%%%% PLOT DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%
pro ree_plot_dist, Dist_str, Tuse=Tuse, Beta=Beta, Ncal=Ncal, maxW=maxW, $
  minW=minW, xran=xran, yran=yran, Wfit1=Wfit1, Wfit2=Wfit2, $
  Wfpi=Wfpi, Ffpi=Ffpi, Wfeeps=Wfeeps, Ffeeps=Ffeeps, Wmaxwell=Wmaxwell, $
  Wax=Wax, fw=fw, Wlog=Wlog, flog=flog, flux=flux, fluxL=fluxL

; NATURAL NUMBERS
ee    = 1.6d-19
mi    = 1.67d-27
Re    = 6374.d3

; KEYWORDS
if not keyword_set(Tuse) then Tuse = 4d3
if not keyword_set(Beta) then Beta = 1.05D
if not keyword_set(maxW) then maxW = 500D
if not keyword_set(minW) then minW = 0.5D

; UNPACK STRUCTURE
fw   = Dist_Str.fw > 1d-6
fw   = fw/max(fw)
WaxW = Dist_Str.Wax

; MAKE LOG-SPACED DISTRIBUTION
ree_v06_make_log_space, WaxW, fw, beta, Xlog=Wlog, Ylog=flog

; CALIBRATE
VaxW = sqrt(2D*WaxW*ee/mi)
flux = Dist_Str.fw*VaxW/max(Dist_Str.fw)
dW = WaxW(1)-WaxW(0)

IF keyword_set(Ncal) then BEGIN
  Ntest = total(fw)*dW * 4D * !dpi * 1d-3 ; EV to KEV
  Cal = Ncal/Ntest
ENDIF ELSE BEGIN
  Cal = 1D/max(flux)
ENDELSE

; SET FLUX
flux = flux*Cal
fluxL = flog*Cal*sqrt(2D*Wlog*ee/mi)
Wlog  = Wlog*1e-3

; SET YRAN AND XRAN
IF not keyword_set(yran) then BEGIN
  Yran = fltarr(2)
  yran(1) = max(fluxL)*3.0
  yran(0) = yran(1)*1e-5
ENDIF

if not keyword_set(xran) then xran = [minW/2.0, maxW*2.0]

; DO FITS
IF keyword_set(Wfit1) then BEGIN
  indfit1 = where( (Wlog GE Wfit1(0)) AND (Wlog LE Wfit1(1))  )
  xax     = alog(Wlog(indfit1))
  yax     = alog(fluxL(indfit1))
  C1      = ladfit(xax, yax)
  dum     =  strcompress(string(round(-C1(1) * 100.0)), /rem)
  C1str   = '-' + strmid(dum,0,1) + '.' + strmid(dum,1,2)
  if abs(C1(1)) LT 1.0 then C1str  = '-0.' + strmid(dum,0,2)
  print, 'Fit 1: ' + C1str
ENDIF

IF keyword_set(Wfit2) then BEGIN
  indfit2 = where( (Wlog GE Wfit2(0)) AND (Wlog LE Wfit2(1)) )
  xax     = alog(Wlog(indfit2))
  yax     = alog(fluxL(indfit2))
  C2      = ladfit(xax, yax)
  dum     =  strcompress(string(round(-C2(1) * 100.0)), /rem)
  C2str   = '-' + strmid(dum,0,1) + '.' + strmid(dum,1,2)
  if abs(C2(1)) LT 1.0 then C2str  = '-0.' + strmid(dum,0,2)
  print, 'Fit 2: ' + C2str
ENDIF

; SET UP MAXWELLIAN
IF keyword_set(Tuse) AND keyword_set(Wmaxwell) then BEGIN
  indMxwl = where( (WaxW GE Wmaxwell(0)*1D3) AND (WaxW LE Wmaxwell(1)*1D3) )
  WMxwl = WaxW(indMxwl)
  Fraw = exp(-WMxwl/Tuse) * WMxwl/Tuse
  NFcal = total(flux/WaxW)
  NFtest = total(Fraw/WaxW(indMxwl))
  Mcal=NFcal/NFtest
  Mflux =  Fraw*Mcal
ENDIF

; ****** PLOT ******
indW = where( (Wlog GE minW) AND (Wlog LE maxW) )
ree_sym, 0.25
plot, Wlog(indW), fluxL(indW), psym=8, /xlog, /ylog, $
  xtitle = 'Energy (keV)', ytitle = 'Flux (m!U2!N sr s keV)!U-1!N', yran=yran, $
  xran = [minW/2.0, maxW*2.0], /xstyle, /ystyle, title= 'Test_Particle Simulation';, $

; PLOT FITS
IF keyword_set(Wfit1) then BEGIN
  Warr = Wlog(indfit1)
  oplot, Warr, exp(C1(0) + alog(Warr)*C1(1)), col=4
ENDIF

IF keyword_set(Wfit2) then BEGIN
  Warr = Wlog(indfit2)
  oplot, Warr, exp(C2(0) + alog(Warr)*C2(1)), col=190
ENDIF

; PLOT MEASURED DATA
ree_sym, 0.25, /sq
if keyword_set(Wfpi) AND keyword_set(Ffpi) then $
  oplot, Wfpi, Ffpi, col=85, psym=-8
if keyword_set(Wfeeps) AND keyword_set(Ffeeps) then $
  oplot, Wfeeps, Ffeeps, col=6, psym=-8

; PLOT MAXWELLIAN
if keyword_set(Tuse) AND keyword_set(Wmaxwell) then $
  oplot, WMxwl*1e-3, Mflux, col=4
  
  RETURN
end


function ddsega_plot_trajectory, trajectories, number, box_str, xy=xy, xz=xz, yz=yz, ddd = ddd, vx = vx, vy = vy, vz = vz, fields = fields, overplot=overplot, scalebar=scalebar
  Re    = 6374.d3

  ;

  ;Plots the trajectories produced my a modefied version of run_2dplus05A.pro
  ;It cleans up the arrays where the trajectories are allocated

  ; EXAMPLE:
  ; ddsega_plot_trajectory(trajectories, 21, box_str, yz=1, fields = 1)

  ;print,number
  !p.multi = [0,1,1]
  ;UNPACK BOX STRUCTURE

  xexitmin = box_str.xexitmin
  xexitmax = box_str.xexitmax

  zexitmin = box_str.zexitmin
  zexitmax = box_str.zexitmax

  yexitmin = box_str.yexitmin
  yexitmax = box_str.yexitmax

  goodx = where(trajectories[*, 3, number] NE 0)
  goody = where(trajectories[*, 4, number] NE 0)
  goodz = where(trajectories[*, 5, number] NE 0)

  dpx = trajectories[goodx,3,number]-shift(trajectories[goodx,3,number],1)
  ddpx = dpx-shift(dpx,1)
  dddpx = ddpx-shift(ddpx,1)
  ddddpx = dddpx-shift(dddpx,1)


  fin = where(dpx eq 0 and ddpx eq 0and dddpx eq 0and ddddpx eq 0,m)
  if m gt 1 then fin = fin[0] else fin=-1
  ini=10
  ini=0


  x = trajectories[goodx[ini:fin], 0, number]
  y = trajectories[goody[ini:fin], 1, number]
  z = trajectories[goodz[ini:fin], 2, number]
  px = trajectories[goodx[ini:fin], 3, number]
  py = trajectories[goody[ini:fin], 4, number]
  pz = trajectories[goodz[ini:fin], 5, number]

  Ex = trajectories[goodx[ini:fin], 6, number]
  Ey = trajectories[goody[ini:fin], 7, number]
  Ez = trajectories[goodz[ini:fin], 8, number]
  Bx = trajectories[goodx[ini:fin], 9, number]
  By = trajectories[goody[ini:fin], 10, number]
  Bz = trajectories[goodz[ini:fin], 11, number]

  ;x = trajectories[goodx[ini]:-1, 0, number]
  ;y = trajectories[goody[ini]:-1, 1, number]
  ;z = trajectories[goodz[ini]:-1, 2, number]
  ;px = trajectories[goodx[ini]:-1, 3, number]
  ;py = trajectories[goody[ini]:-1, 4, number]
  ;pz = trajectories[goodz[ini]:-1, 5, number]

  ;Ex = trajectories[goodx[ini]:-1, 6, number]
  ;Ey = trajectories[goody[ini]:-1, 7, number]
  ;Ez = trajectories[goodz[ini]:-1, 8, number]
  ;Bx = trajectories[goodx[ini]:-1, 9, number]
  ;By = trajectories[goody[ini]:-1, 10, number]
  ;Bz = trajectories[goodz[ini]:-1, 11, number]

  if keyword_set(xy) then begin
    wi,1
    plot, x[1:-1], y[1:-1], xrange = [xexitmin,xexitmax], yrange = [yexitmin,yexitmax]
    y_coor = arrgen(yexitmin,yexitmax,10000)
    xplus = dblarr(10000) + xexitmax
    xminus = dblarr(10000) + xexitmin

    oplot,xplus,y_coor;,linestyle = 'dot'
    oplot,xminus,y_coor;, linestyle = 'dot'

    x_coor = arrgen(xexitmin,xexitmax,10000)
    yplus = dblarr(10000) + yexitmax
    yminus = dblarr(10000) + yexitmin

    oplot,x_coor, yplus
    oplot,x_coor, yminus



  endif
  if keyword_set(xz) then plot, x[1:-1], z[1:-1], xrange = [xexitmin,xexitmax], yrange = [zexitmin,zexitmax]
  if keyword_set(yz) then begin
    plot, y[1:-1], z[1:-1], xrange = [yexitmin,yexitmax], yrange = [zexitmin,zexitmax]

    y_coor = arrgen(yexitmin,yexitmax,10000)
    zplus = dblarr(10000) + zexitmax
    zminus = dblarr(10000) + zexitmin

    oplot,y_coor, zplus
    oplot,y_coor, zminus
  endif
  if keyword_Set(ddd) then begin
    ; wi,10,wsize=[1600,1200]
    factor=1
    if keyword_Set(overplot) and keyword_Set(scalebar) then factor=4.
    if keyword_Set(overplot) then begin
      graphic1 = plot3d(-x[1:-1]/re, y[1:-1]/re, -z[1:-1]/re,depth_cue=[0,4]$
        ,dimensions=[1100,1100],margin=[0.2,0.33,0.2,0.2],/perspective,buffer=0,AXIS_STYLE=2,thick=4,rgb_table=2,VERT_COLORS=BYTSCL(findgen(n_elements(x[1:-1])*factor)), $
        SHADOW_COLOR="red",xrange = [xexitmin/re,xexitmax/re], yrange = [yexitmin/re,yexitmax/re], zrange = [zexitmin/re,zexitmax/re],clip=0,overplot=1)

      c = COLORBAR(range= [0,max(indgen(n_elements(x[1:-1])))/(100*60.)], ORIENTATION=1, POSITION=[0.972,0.16,0.992,0.48],rgb_table=2,ticklayout=1,font_size=0)

      ax = graphic1.AXES
      ax[2].HIDE = 1
      ax[6].HIDE = 1
      ax[7].HIDE = 1

      graphic1.XY_SHADOW=1
      ; graphic1.Xz_SHADOW=1
      ;graphic1.YZ_SHADOW=1

      Zaxis = AXIS('Z', LOCATION=['left','top'], TITLE='$Z [R_E]$', AXIS_RANGE=[zexitmin/re,zexitmax/re])

      graphic1.Font_size=21


    endif else begin

      ;graphic1 = plot3d(x[1:-1]/re, y[1:-1]/re, -z[1:-1]/re, xtitle='$X [R_E]$',ytitle='$Y [R_E]$',ztitle='$Z [R_E]$',depth_cue=[0,4]$
      ;  ,dimensions=[1100,1100],margin=[0.2,0.33,0.2,0.2],/perspective,buffer=0,AXIS_STYLE=2,thick=3,rgb_table=16,VERT_COLORS=BYTSCL(findgen(n_elements(x[1:-1]))), $
      ;  SHADOW_COLOR="blue",xrange = [xexitmin/re,xexitmax/re], yrange = [yexitmin/re,yexitmax/re], zrange = [zexitmin/re,zexitmax/re],clip=0)


      graphic1 = plot3d(-x[1:-1]/re, y[1:-1]/re, z[1:-1]/re, xtitle='$X [R_E]$',ytitle='$Y [R_E]$',ztitle='$Z [R_E]$',depth_cue=[0,4]$
        ,dimensions=[1100,1100],margin=[0.2,0.33,0.2,0.2],/perspective,buffer=0,AXIS_STYLE=2,thick=4,rgb_table=16,VERT_COLORS=BYTSCL(findgen(n_elements(x[1:-1]))), $
        SHADOW_COLOR="blue",xrange = [xexitmin/re,xexitmax/re], yrange = [yexitmin/re,yexitmax/re], zrange = [zexitmin/re,zexitmax/re],clip=0)
      ax = graphic1.AXES
      ax[2].HIDE = 1
      ax[6].HIDE = 1
      ax[7].HIDE = 1

      ;graphic1.XY_SHADOW=1
      ;graphic1.Xz_SHADOW=1
      graphic1.YZ_SHADOW=1

      Zaxis = AXIS('Z', LOCATION=['left','top'], TITLE='$Z [R_E]$', AXIS_RANGE=[zexitmin/re,zexitmax/re])

      graphic1.Font_size=21

      c = COLORBAR(range= [0,max(indgen(n_elements(x[1:-1])))/(100*60.)], ORIENTATION=1, TITLE='Time [min]',font_size=18, POSITION=[0.95,0.16,0.97,0.48],ticklayout=1)
    endelse

    ;  graphic1 = plot3d(x[1:-1]/xexitmax, y[1:-1]/yexitmax, z[1:-1]/zexitmax, xtitle='X',ytitle='Y',ztitle='Z',depth_cue=[0,4]$
    ;    ,dimensions=[1820,1100],margin=[0.2,0.33,0.2,0.2],/perspective,buffer=0,AXIS_STYLE=2,thick=3,rgb_table=33,VERT_COLORS=BYTSCL(n_elements(x[1:-1])), $
    ;    SHADOW_COLOR="deep sky blue",xrange = [-1,1], yrange = [-1,1], zrange = [-1,1],clip=0)

    cd,'dump3d'
    graphic1.save, 'trajectories3D'+string(number)+'.png',resolution=1000,bitmap=1
    cd,'..'

  endif

  if keyword_set(vx) then begin
    wi,0
    plot, px, charsize = 2
  endif

  if keyword_set(vy) then begin
    wi,2
    plot, py, charsize = 2
  endif

  if keyword_set(vz) then begin
    wi,3
    plot, pz, charsize = 2
  endif

  if keyword_set(fields) then begin
    wi, 2,wsize = [1320,720]

    !P.MULTI=[0,1,3]


    cgplot, (px^2 + py^2 + pz^2)/2., charsize = 3.5, yrange = [-8e12, 8e12] ,thick = 2
    ;oplot, py, color=cgcolor('green'),thick = 2
    ;oplot, pz, color=cgcolor('blue'),thick = 2

    cgplot, sqrt(Ex^2 + Ey^2 + Ez^2), nsum=2, thick = 2, charsize = 3.5, yrange = [-30, 30]
    ;cgoplot, Ey, nsum=50, color=cgcolor('green'), thick = 2
    ;cgoplot, Ez, nsum=50, color=cgcolor('blue') ,thick = 2

    cgplot, Bx, thick =2, charsize = 3.5
    cgoplot, By, color=cgcolor('green'), thick =2
    cgoplot, Bz, color=cgcolor('blue'), thick =2
    ;plot,ex

    cgText, 0.988, 0.94, 'px', ALIGNMENT=0.5, /NORMAL, charsize=2.5
    cgText, 0.988, 0.92, 'py', ALIGNMENT=0.5, /NORMAL, color= cgcolor('green'), charsize=2.5
    cgText, 0.988, 0.9, 'pz', ALIGNMENT=0.5, /NORMAL, color= cgcolor('blue'), charsize=2.5

    cgText, 0.988, 0.58, 'Ex', ALIGNMENT=0.5, /NORMAL, charsize=2.5
    cgText, 0.988, 0.56, 'Ey', ALIGNMENT=0.5, /NORMAL, color= cgcolor('green'), charsize=2.5
    cgText, 0.988, 0.54, 'Ez', ALIGNMENT=0.5, /NORMAL, color= cgcolor('blue'), charsize=2.5

    cgText, 0.988, 0.26, 'Bx', ALIGNMENT=0.5, /NORMAL, charsize=2.5
    cgText, 0.988, 0.24, 'By', ALIGNMENT=0.5, /NORMAL, color= cgcolor('green'), charsize=2.5
    cgText, 0.988, 0.22, 'Bz', ALIGNMENT=0.5, /NORMAL, color= cgcolor('blue'), charsize=2.5

    ; px = trajectories[goodx[0:fin], 3, number]
    ; py = trajectories[goody[0:fin], 4, number]
    ; pz = trajectories[goodz[0:fin], 5, number]


    Energy = px[where(abs(Bz) lt 1e-9)]^2 +   py[where(abs(Bz) lt 1e-9)]^2 +   pz[where(abs(Bz) lt 1e-9)]^2

    Energy1 = (px[0]^2 +   py[0]^2 +   pz[0]^2)/2
    Energy2 = (px[-1]^2 +   py[-1]^2 +   pz[-1]^2)/2


    ;derivative = (energy - shift(Energy,1))/0.01
    ;derivative = derivative[2:-1]

    ;print,energy
    ;print,'derivative',derivative

    print,'Initial Energy  =',energy1
    print,'Final Energy  =',energy2





    ;cgLegend, Color=['red','green', 'blue'], Location=[0.73, 0.85],  Titles=['X-band', 'K-band','S-band'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color='white',charsize=2.5
    ;cgLegend, Color=['red','magenta', 'blue'], Location=[0.73, 0.85],  Titles=['X-band', 'K-band','S-band'], Length=0.075, /Box, VSpace=2.75,/Background, bg_color='white',charsize=2.5


    ;plot, py, charsize = 2, color=0
    ;oplot, Ex*1e5, color=cgcolor('red'), nsum=50, thick = 2
    ;oplot, Ey*1e5, color=cgcolor('green'), nsum=50, thick = 2
    ;oplot, Ez*1e5, color=cgcolor('blue'), nsum=50, thick = 2
    ;oplot, Bx*1e14, color=cgcolor('cyan'), thick =2
    ; oplot, By*1e14, color=cgcolor('Magenta'), thick =2
    ;  oplot, Bz*1e14, color=cgcolor('orange'), thick =2


  endif
  if keyword_set(graphic1) then return, graphic1
  end

function dsega_energy_flux, radiation, box_str

  ;UNPACK BOX STRUCTURE

  xexitmin = box_str.xexitmin
  xexitmax = box_str.xexitmax

  zexitmin = box_str.zexitmin
  zexitmax = box_str.zexitmax

  yexitmin = box_str.yexitmin
  yexitmax = box_str.yexitmax

  ixplus  = where(radiation[*, 0] GT xexitmax, nxp)
  ixminus = where(radiation[*, 0] LT xexitmin, nxm)
  iyplus  = where(radiation[*, 1] GT yexitmax, nyp)
  iyminus = where(radiation[*, 1] LT yexitmin, nym)
  izplus  = where(radiation[*, 2] GT zexitmax, nzp)
  izminus = where(radiation[*, 2] LT zexitmin, nzm)

  flux = {xminus:make_array(nxm+1), xplus:make_array(nxp+1), yminus:make_array(nym+1),$
    yplus:make_array(nyp+1),zminus:make_array(nzm+1), zplus:make_array(nzp+1)}


  flux.xminus  = radiation[ixminus, 3]
  flux.xplus   = radiation[ixplus,  3]
  flux.yminus  = radiation[iyminus, 4]
  flux.yplus   = radiation[iyplus,  4]
  flux.zminus  = radiation[izminus, 5]
  flux.zplus   = radiation[izplus,  5]




  ;endfor

  RETURN, flux



end