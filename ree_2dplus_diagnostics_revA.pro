; %%%%%%%%%%%%%%%%%%%%%%% CURRENT MEASURMENT %%%%%%%%%%%%%%%%%%%%%%%%%%
pro ree_y_current, z, py, gam, ind, Box_Str, plot=plot, ddist=ddist, $
  Jy=Jy, zax=za
  
; CHECK KEYWORDS
if not keyword_set(ddist) then ddist=1d6
Re    = 6374.d3
ee    = 1.6d-19
n0    = 0.05d6
mu0   = 4d-7*!dpi

nind  = n_elements(ind)
npts  = round((Box_Str.zmax - Box_Str.zmin)/ddist)
binz  = (Box_Str.zmax - Box_Str.zmin)/npts
zax   = dindgen(npts)*binz + Box_Str.zmin
zind  = (round( (z(ind) - Box_Str.zmin)/binz ) > 0 ) < (npts-1)
nz    = histogram(z(ind), bin=binz, min=Box_Str.zmin, max=Box_Str.zmax-1D)
cal   = ee*n0/max(nz)
Jy    = dblarr(npts)
for ii=0, nind-1 do jy(zind(ii)) = jy(zind(ii)) + cal*py(ind(ii))/gam(ind(ii))

if keyword_set(plot) then plot, zax/Re, Jy*1d9, $
  xtit='Z (R!DE!N)', ytit = 'Current (nA/m!U2!N)'
if keyword_set(plot) then print, 'DeltaB = ', $
  strcompress(string(total(mu0*total(Jy)*binz*1d9))), ' nT'

RETURN
end

; %%%%%%%%%%%%%%%%%%%%%%% PARTICLE PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%
pro ree_run_time_plot, x, y, z, i, ind, nind, Box_Str 

Re    = 6374.d3
ree_sym, 1.0, /fill

plot, x(ind)/Re, z(ind)/Re, psym=3, xran=[Box_Str.xmin, Box_Str.xmax]/Re, $
    yran = [Box_Str.zmin,Box_Str.zmax]/Re, xtit = 'X (R!DE!N)', $
    ytit = 'Z (R!DE!N)'
oplot, x(1000:1001)/Re, z(1000:1001)/Re, psym=8, col=2
oplot, x(1002:1003)/Re, z(1002:1003)/Re, psym=8, col=4
oplot, x(9000:9001)/Re, z(9000:9001)/Re, psym=8, col=6 
  
wi, 1
plot, y(ind)/Re, z(ind)/Re, psym=3, xran=[Box_Str.ymin, Box_Str.ymax]/Re, $
      yran = [Box_Str.zmin,Box_Str.zmax]/Re, xtit = 'Y (R!DE!N)', $
      ytit = 'Z (R!DE!N)', /xstyle
oplot, y(1000:1001)/Re, z(1000:1001)/Re, psym=8, col=2
oplot, y(1002:1003)/Re, z(1002:1003)/Re, psym=8, col=4
oplot, y(9000:9001)/Re, z(9000:9001)/Re, psym=8, col=6
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
TallB  = mean(sqrt(px(ind)*px(ind)+py(ind)*py(ind)+pz(ind)*pz(ind)+c2)-c) * $
                 mi*c/ee*2D/3D
Btot = sqrt(Bx(ind)*Bx(ind) + By(ind)*By(ind) + Bz(ind)*Bz(ind))
ppll = (px(ind)*Bx(ind) + py(ind)*By(ind) + pz(ind)*Bz(ind) ) / Btot 
TparB  = mean(sqrt(ppll*ppll+c2)-c)*2D*mi*c/ee
TperpB = (3D*TallB - TparB)/2D

print, TparB/1000D, TperpB/1000D, TallB/1000D, format =  $
 '("  TBoxPar  = ", F9.3,"  TBoxPerp  = ", F9.3, "  TBoxAve  = ", F9.3, " keV")'

indd = where( $
  (x(ind) GT Box_str.xturbmin) AND (x(ind) LT Box_str.xturbmax) AND $
  (y(ind) GT Box_str.yturbmin) AND (y(ind) LT Box_str.yturbmax) AND $
  (z(ind) GT Box_str.zturbmin) AND (z(ind) LT Box_str.zturbmax) ) & $

indt = ind(indd)

TallT  = mean(sqrt(px(indt)*px(indt)+py(indt)*py(indt)+ $
                      pz(indt)*pz(indt)+c2)-c) * mi*c/ee*2D/3D
Btot = sqrt(Bx(indt)*Bx(indt) + By(indt)*By(indt) + Bz(indt)*Bz(indt))
ppll = (px(indt)*Bx(indt) + py(indt)*By(indt) + pz(indt)*Bz(indt) ) / Btot 
TparT  = mean(sqrt(ppll*ppll+c2)-c)*2D*mi*c/ee
TperpT = (3D*TallT - TparT)/2D

print, TparT/1000D, TperpT/1000D, TallT/1000D, format =  $
 '("  TTurbPar  =", F9.3,"  TTurbPerp  =", F9.3, "  TTurbAve  =", F9.3, " keV")'

; SKIP
print, ' '

RETURN
end
;%%%%%%%%%%%%%%%%%%%%%%% Record Flux %%%%%%%%%%%%%%%%%%%%
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

function ddsega_plot_trajectory, trajectories, number, box_str, xy=xy, xz=xz, yz=yz, ddd = ddd, vx = vx, vy = vy, vz = vz, fields = fields

;Plots the trakectories produced by a modefied version of run_2dplus05A.pro
;It cleans up the arrays where the trajectories are allocated

;Examples:
;
;ddsega_plot_trajectory(trajectories,1,box_str,xy=1)
;ddsega_plot_trajectory(trajectories,55,box_str,xy=1,fields=1)

!p.multi = [0,1,1]
;UNPACK BOX STRUCTURE

xexitmin = box_str.xexitmin
xexitmax = box_str.xexitmax

zexitmin = box_str.zexitmin
zexitmax = box_str.zexitmax

yexitmin = box_str.yexitmin
yexitmax = box_str.yexitmax

goodx = where(trajectories[*, 0, number] NE 0)
goody = where(trajectories[*, 1, number] NE 0)
goodz = where(trajectories[*, 2, number] NE 900)

dpx = trajectories[goodx,3,number]-shift(trajectories[goodx,3,number],1)
fin = where(dpx eq 0)
fin = fin[1]


x = trajectories[goodx[0:fin], 0, number]
y = trajectories[goody[0:fin], 1, number]
z = trajectories[goodz[0:fin], 2, number]
px = trajectories[goodx[0:fin], 3, number]
py = trajectories[goody[0:fin], 4, number]
pz = trajectories[goodz[0:fin], 5, number]

Ex = trajectories[goodx[0:fin], 6, number]
Ey = trajectories[goody[0:fin], 7, number]
Ez = trajectories[goodz[0:fin], 8, number]
Bx = trajectories[goodx[0:fin], 9, number]
By = trajectories[goody[0:fin], 10, number]
Bz = trajectories[goodz[0:fin], 11, number]

wi,1
if keyword_set(xy) then begin
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
  graphic = plot3d(x[1:-1], y[1:-1], z[1:-1])
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
  wi, 2, wsize = [1500,1200]
  
  !P.MULTI=[0,1,3]

  
   plot, px, charsize = 3.5, yrange = [-8e6, 8e6] ,thick = 2
  oplot, py, color=cgcolor('green'),thick = 2
  oplot, pz, color=cgcolor('blue'),thick = 2
  
  cgplot, Ex, nsum=50, thick = 2, charsize = 3.5, yrange = [-30, 30]
  oplot, Ey, nsum=50, color=cgcolor('green'), thick = 2
  oplot, Ez, nsum=50, color=cgcolor('blue') ,thick = 2
  
  plot, Bx, thick =2, charsize = 3.5
  oplot, By, color=cgcolor('green'), thick =2
  oplot, Bz, color=cgcolor('blue'), thick =2
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

;return, px




end







