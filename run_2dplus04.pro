; PARTICLES IN A MAGNETIC RECONNECTION REGION

; SETUP: LINK TO THE SPEDAS LIBRARY SETS PATH TO MMS DATA
;setenv, 'ROOT_DATA_DIR=/Users/ree/mmsdata/'
;mms_init, local_data_dir = '/Users/ree/mmsdata'

; COMPILE RAND_E 
; #### RESET PATH #####  
;.r '/Users/ree/Bob/papers/20_IonAcceleration/idl/ree_2dplus_lib'
;.r '/Users/ree/Bob/papers/20_IonAcceleration/idl/ree_2dplus_ion_sim'
;.r '/Users/ree/Bob/papers/20_IonAcceleration/idl/ree_2dplus_diagnostics'

; FIX COLORS ###### WILL NOT WORK ON WINDOWS #######
;loadct2, 43
;tvlct, red, green, blue, /get
;red(6) = 225
;green(4) = 135
;red(5) = 200
;green(5) = 175
;tvlct, red, green, blue

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
dt     = 1d-2
qmdt   = ee/mi*dt/1000.0 ; E in mV/m
bqmdt  = ee/mi*dt ; FOR dB

; SET UP E and T
TempStart = 4000D ; eV
vth       = sqrt(TempStart*ee/mi)

; SET UP SIMULATION DOMAIN
@2dplus_setup04.txt

; ****** RESET POINT *******
it  = 0L ; TRACKS POSITION IN ELECTRIC FIELD ARRAYS
kt  = 0L ; TRACKS POSITION IN MAGNETIC FIELD ARRAYS
jt  = 0L ; TRACKS POSITION IN PARTICLE INJECTION ARRAYS
k   = 0L ; TRACKS DIAGNOSTIC ARRAYS


; SETUP PARTICLE ARRRAYS
npart = 100000L; MAXIMUM NUMBER OF ACTIVE PARTICLES
Ex   = dblarr(npart)
Ey   = dblarr(npart)
Ez   = dblarr(npart)
Xpwr = dblarr(npart)
Ypwr = dblarr(npart)
Zpwr = dblarr(npart)
dsE  = lonarr(npart)
Bx   = dblarr(npart)
By   = dblarr(npart)
Bz   = dblarr(npart)
dbx  = dblarr(npart)
dby  = dblarr(npart)
dbz  = dblarr(npart)
dsB  = lonarr(npart)
rat  = dblarr(npart) ; DUMMY ARRAY TO SAVE RUN TIME
Bmag = dblarr(npart) ; NOT NEEDED SO FAR
x    = dblarr(npart)
y    = dblarr(npart)
z    = dblarr(npart)
dx   = dblarr(npart)
px   = dblarr(npart)
py   = dblarr(npart)
pz   = dblarr(npart)
pxn  = dblarr(npart) ; DUMMY ARRAY TO SAVE RUN TIME
pyn  = dblarr(npart) ; DUMMY
pzn  = dblarr(npart) ; DUMMY
pxold= dblarr(npart) ; DUMMY
pyold= dblarr(npart) ; DUMMY
pzold= dblarr(npart) ; DUMMY
gam  = dblarr(npart)


; LOAD PARTICLES
nload = 22000L
nstart = npart-nload
load_particles3D, nload, vth, DenZ0, DenInjZ, nEtimes, nBtimes, $
  Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, $
  x=x, y=y, z=z, px=px, py=py, pz=pz, gam=gam, dsE=dsE, dsB=dsB, nstart=nstart  

; CONTROLL PARAMETERS
imax = 5000L    ; NUMBER OF TIME STEPS
nnew = 5       ; NUMBER OF INJECTED PARTICLES PER TIME STEP

; DIAGNOSTICS REFRESH (IF NEEDED)
plot_mod    = 0
plot_step   = 10
make_density = 0
make_image   = 0
make_dist    = 0

; KNOBS
Eplevel = 1.0
Bplevel = 1.0
Ey0 = 2.7D; 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 
B_str.BH0 = 16D-9

; MAIN LOOP
ree_2dplus_ion_sim, dt, it, jt, kt, npart, imax, nnew, nh, box_str, $
  B_str, Inj_str, x, y, z, px, py, pz, $
  Ex, Ey, Ez, Ey0, Epll, E1, E2, dsE, nEtimes, Eplevel, $
  Bx, By, Bz, Bmag, dBx, dBy, dBz, dB1, dB2, dB3, dsB, nBtimes, Bplevel, $
  pxn, pyn, pzn, pxold, pyold, pzold, gam, rat, plot_step, plot_mod, $
  Xpwr, Ypwr, Zpwr, psym, vth, Pdx, ind=ind

; INDICATOR OF HOW MANY ORIGONAL PARTICLES REMAIN
help, where(ind GE nstart)    
ree_print_temps, x, y, z, px, py, pz, ind, Box_str

; PLOT DISTRIBUTION
wi, 2
ree_distribution_basic, x, y, z, px, py, pz, ind, Box_str, /plot, $ 
  fpx=fpx, fpy=fpy, fpz=fpz, pax=pax, fave=fave, Wax=Wax, $
  fabspx=fabspx, fabspy=fabspy, fabspz=fabspz, pabs=pabs, /turb
Ti = 6.25d3
oplot, wax, max(fave)*exp(-Wax/Ti), col=2
wi,0

wi, 2
plot, pax, float(fpy)>0.5,/ylog, xran = [-1d7, 1d7], $
  xtit='py/m!Di!N (m/s)', ytit = 'Counts'
oplot, pax, max(fpy)*exp(-0.5D*pax*pax*mi/Ti/ee), col=2
wi,0

wi, 2
plot, pax, float(fpx)>0.5,/ylog, xran = [-1d7, 1d7], $
  xtit='px/m!Di!N (m/s)', ytit = 'Counts'
oplot, pax, max(fpx)*exp(-0.5D*pax*pax*mi/Ti/ee), col=2
wi,0

wi, 2
plot, pax, float(fpz)>0.5,/ylog, xran = [-1d7, 1d7], $
  xtit='pz/m!Di!N (m/s)', ytit = 'Counts'
oplot, pax, max(fpz)*exp(-0.5D*pax*pax*mi/Ti/ee), col=2
wi,0


W = W0*(gam - 1D)
ree_sym, 0.1
plot, y(ind)/Re, W(ind)/1d3, psym=8, xran=[ymin, ymax]/Re, ytit='Energy (keV)', $
  xtit = 'Y (R!DE!N)' 
oplot, y(ind)/Re, (y(ind)-box_str.ymin)*Ey0*1d-6, psym=3, col=6



;plot, z(ind), Bx(ind), psym=3
;plot, x(ind), Bz(ind), psym=3



end
