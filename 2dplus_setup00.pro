
; SET UP BOX
xmin     = -6.25d7 ; LOCATION OF INJECTION BOUNDARY (PARTICLE ENTRY)
xmax     =  6.25d7 ; ABOUT 10 RE
xexitmin = -6.25d7 ; LOCATION PARTICLE EXIT BOUNDARY
xexitmax =  6.25d7
xturbmin = -5.25d7 ; LOCATION OF TURBULENCE REGION
xturbmax = 5.25d7

zmin     = -1.75d7  ; ABOUT 3 RE
zmax     =  1.75d7
zexitmin =  zmin
zexitmax =  zmax
zturbmin = -1.25d7  ; ABOUT 2 RE
zturbmax =  1.25d7

ymin     = -2.0d7
ymax     =  2.0d7
yexitmin = ymin
yexitmax = ymax
yturbmin = -1.5d7
yturbmax =  1.5d7

; BOX STRUCTURE CONTAINS NEEDED INFORMATION OF THE BOX
box_str = {xmin:xmin, xmax:xmax, xexitmin:xexitmin, xexitmax:xexitmax, $
           xturbmin:xturbmin, xturbmax:xturbmax, zmin:zmin, zmax:zmax, $
           zexitmin:zexitmin, zexitmax:zexitmax, zturbmin:zturbmin, $
           zturbmax:zturbmax, ymin:ymin, ymax:ymax, yexitmin:yexitmin, $
           yexitmax:yexitmax, yturbmin:yturbmin, yturbmax:yturbmax}
           
; SET UP MAGNETIC FIELD
Bx0 = 20d-9
Bz0 = 0d-9

;Bz0 = 2.5d-9  ; 8:1 RATIO
X0 = 1.00d7
;Z0 = 1.25d6 ; IMPORTANT! Z0 IS ALSO USED IN DENSITY PROFILE
Z0  = 6.5d6
;Bz = Bz0*tanh(x/x0) & $
;Bx = Bx0*tanh(z/z0) & $

; DELTA B - DISABLED
db0 = 32D-9
;db0 = 32D-9
Bbeta  = !dpi*dt/4D
Balpha = 1D - Bbeta

;  PUT IN A MAGNETIC HOLE NOT ACTIVATED HERE
zH0 = 0.35D7
zHW = 0.25D7
xHW = 2.5D7
;BH0 = 16D-9
BH0 = 0D-9 ; NO MAGNTIC HOLE

B_str = {Bx0:Bx0, By0:0D, Bz0:Bz0, X0:X0, Y0:0D, Z0:Z0, $
         BH0:BH0, xHw:xHw, zHw:xHw, zH0:xHw, $
         db0:db0, Bbeta:Bbeta, Balpha:Balpha}

; SET UP EFIELD
Pdx   = 0.5D*Re

; SETUP - TAKES OVER A MINUTE TO GENERATE RANDOM E WITH MEASURED SPECTRA 
nEtimes = 100000L
ion_e_set_up, 0.01, 50D, nfreqs=5000L, Erms=10.5D, $
  alpha=1.0D, ntimes=nEtimes, EIonT=E1
ion_e_set_up, 0.01, 50D, nfreqs=5000L, Erms=10.5D, $
  alpha=1.0D, ntimes=nEtimes, EIonT=E2
ion_e_set_up, 0.01, 50D, nfreqs=5000L, Erms=5D, $
  alpha=0.3D, ntimes=nEtimes, EIonT=Epll
Ey0 = 2.7D
it  = 0L ; TRACKS POSITION IN ELECTRIC FIELD ARRAYS

; SET UP PARTICLE INJECTION - VARIES WITH Z  
ninjarr  = 100000L
DenInjZ  = 0.2D ; ESTABLISH A BACKGROUND DENSITY AT Zmax & Zmin
DenZ0    = Z0 ; NEEDED TO BALANCE B PRESSURE

; PRE-MAKE INJECTION ARRAYS
inject_particles3D_setup, ninjarr, vth, DenZ0, DenInjZ, nEtimes, $
   Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Inj_str=Inj_str
jt = 0L ; TRACKS POSITION IN PARTICLE INJECTION ARRAYS

; TEST HISTROGRAM 
dz       = (zmax-zmin)/100.0
zax      = (dindgen(100)-49.5)*dz
hist = histogram(Inj_str.z, min=zmin, max=zmax*0.9999999999999D, bin = dz)
plot, zax, hist
DenExtra = ( DenInjZ - 1D/(cosh(Zmax/DenZ0)^2) ) > 0D
Zinside = zax(1:98)
plot, Zinside, float(hist(1:98))/max(hist(1:98))
oplot, Zinside, (DenExtra + (1D - DenExtra)/(cosh(Zinside/DenZ0)^2) ), col=2


; TEST PRESSURE BALANCE
ZdistEff = 2D*tanh(Zmax/DenZ0)*DenZ0 ; INTEGRATION OF SECH^2
thresh = DenExtra*(Zmax-Zmin)/(ZdistEff*(1D - DenExtra) + $
DenExtra*(Zmax-Zmin))
n0 = 0.25d6/(1D - denextra)
BxTest = B_str.Bx0*tanh(Zinside/B_str.Z0)
plot, Zinside, BxTest
PressureB = BxTest*Bxtest/(2D*mu0) * 1d9
PressureN = n0*TempStart*ee/(cosh(Zinside/Z0)^2) * (1D - denextra) * 1d9 + $
             n0*denextra*TempStart*ee * 1d9
PressureNs = n0*TempStart*ee*hist(1:98)/max(hist(1:98)) * 1d9
plot, Zinside, PressureB, yran = [0, 0.3]
oplot, Zinside, PressureN, col=2
oplot, Zinside, PressureNs, col=4
oplot, Zinside, PressureB+PressureNs, col=6



end





