; PARTICLES IN A MAGNETIC RECONNECTION REGION

; SETUP: LINK TO THE SPEDAS LIBRARY SETS PATH TO MMS DATA
;setenv, 'ROOT_DATA_DIR=/Users/ree/mmsdata/'
;mms_init, local_data_dir = '/Users/ree/mmsdata'

; COMPILE RAND_E 
; #### RESET PATH #####  
;.r '/Users/ree/Bob/papers/20_IonAcceleration/idl/ree_2dplus_lib'
;.r '/Users/ree/Bob/papers/20_IonAcceleration/idl/ree_2dplus_ion_sim'

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
@2d_plus_setup00.txt

; ELECTRUC FIELD CONTROLS
plevel = 0.00D ; POWER LEVEL OF TURBULENT E
Ey0    = 0.00D ; Y ELECTRIC FIELD
Pdx    = 0.5D*Re

; BELOW PLOTS TURBULENCE WINDOWS FOR CHECK
;xax = (dindgen(1250)-625d)*1d5 
;Xpwr  = -0.5D*cos( ( ((xax+xturbmax)/Pdx) < (!dpi) ) > 0D ) - $
;        0.5D*cos( ( ((xturbmax-xax)/Pdx) < (!dpi) ) > 0D )
;plot, xax, xpwr, yran = [-0.1, 1.1]
;Zpwr  = -0.5D*cos( ( ((xax+zturbmax)/Pdx) < (!dpi) ) > 0D ) - $
;        0.5D*cos( ( ((zturbmax-xax)/Pdx) < (!dpi) ) > 0D )
;plot, xax, Zpwr, yran = [-0.1, 1.1]


; ****** RESET POINT *******
it  = 0L ; TRACKS POSITION IN ELECTRIC FIELD ARRAYS
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
ds   = lonarr(npart)
Bx   = dblarr(npart)
By   = dblarr(npart)
Bz   = dblarr(npart)
dbx  = dblarr(npart)
dby  = dblarr(npart)
dbz  = dblarr(npart)
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
load_particles3D, nload, vth, DenZ0, DenInjZ, nEtimes, $
  Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, $
  x=x, y=y, z=z, px=px, py=py, pz=pz, gam=gam, ds=ds, nstart=nstart  

; CONTROLL PARAMETERS
imax = 1000L    ; NUMBER OF TIME STEPS
nnew = 5       ; NUMBER OF INJECTED PARTICLES PER TIME STEP

; DIAGNOSTICS REFRESH (IF NEEDED)
make_density = 0
make_image   = 0
make_dist    = 0
plotXZ       = 100


; PLOT PARAMETERS
ree_sym, 1.00, /fill
psym = 8
plevel = 0.0
Ey0 = 0.0D, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8,
yplotrange = ymax

; MAIN LOOP
ree_2dplus_ion_sim, dt, it, jt, npart, imax, nnew, nh, box_str, B_str, $
  Inj_str, x, y, z, px, py, pz, Ex, Ey, Ez, Ey0, Epll, E1, E2, Bx, By, Bz, $
  Bmag, dBx, dBy, dBz, pxn, pyn, pzn, pxold, pyold, pzold, gam, ds, rat, $
  Xpwr, Ypwr, Zpwr, psym, vth, Pdx, plevel, nEtimes, ind=ind, $
  yplotrange = yplotrange  

; INDICATOR OF HOW MANY ORIGONAL PARTICLES REMAIN
help, where(ind GE nstart)    


;ree_sym, 0.1
;plot, x(ind)/Re, z(ind)/Re, psym=8
;plot, y(ind)/Re, z(ind)/Re, psym=8

ree_print_temps, x, y, z, px, py, pz, ind, Box_str

end
