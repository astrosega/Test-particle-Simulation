; PARTICLES IN A MAGNETIC RECONNECTION REGION SMALLER Y DIRECTION

; SETUP: LINK TO THE SPEDAS LIBRARY SETS PATH TO MMS DATA
;setenv, 'ROOT_DATA_DIR=/Users/ree/mmsdata/'
;mms_init, local_data_dir = '/Users/ree/mmsdata'

; COMPILE RAND_E
; #### RESET PATH #####
;.r 'ree_2dplus_lib_B'
;.r 'idl/ree_2dplus_ion_sim_B'
;.r '/Users/ree/Bob/papers/20_IonAcceleration/idl/ree_2dplus_diagnostics_B'

; FIX COLORS ###### WILL NOT WORK ON WINDOWS #######
;ree_v06_loadct

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
@2dplus_setup07.txt

; ****** RESET POINT *******
it  = 0L ; TRACKS POSITION IN ELECTRIC FIELD ARRAYS
kt  = 0L ; TRACKS POSITION IN MAGNETIC FIELD ARRAYS
jt  = 0L ; TRACKS POSITION IN PARTICLE INJECTION ARRAYS
k   = 0L ; TRACKS DIAGNOSTIC ARRAYS

; RUN-TIME TEMPERATURE PROFILE
Ntemps = 10000L
ktemp  = 0L
Buff   = 0.5D
Temp_Step = 25L
Temp_Str = {Ntemps:Ntemps, kTemp:kTemp, Temp_Step:Temp_Step, Buff:Buff, $
  t:dblarr(Ntemps), TboxPar:dblarr(Ntemps), TboxPerp:dblarr(Ntemps), $
  TboxAve:dblarr(Ntemps), TturbPar:dblarr(Ntemps), TturbPerp:dblarr(Ntemps), $
  TturbAve:dblarr(Ntemps)}

; RUN-TIME FLUX STRUCTURE
Nflux = 1000L
dW = 500D
Flux_Str = {dW:dw, Nflux:Nflux, px:dblarr(Nflux), mx:dblarr(Nflux), $
  py:dblarr(Nflux), my:dblarr(Nflux), pz:dblarr(Nflux), mz:dblarr(Nflux)}

; RUN-TIME DISTRIBUTION FUNCTIONS
Pmax  = 9.99999d7
Bin   = 1d5
npts  = ceil(2*Pmax/Bin)
Turb  = 1L
pax   = dindgen(npts)*Bin - Bin*(npts-1)/2
Wmax  = 3.9999999D6
Wbin  = 200D
Wnpts = 20000L
Dist_Step = 100L
Wax = (dindgen(Wnpts) + 0.5D)*Wbin
BuffSize = 0.5D
JetSize = 0D
Ystart  = -10D
Dist_Str = {Pmax:Pmax, Bin:Bin, npts:npts, Turb:Turb, pax:pax, $
  BuffSize:BuffSize, JetSize:JetSize, Ystart:Ystart, $
  Dist_Step:Dist_Step, Wmax:Wmax, Wbin:Wbin, Wnpts:Wnpts, Wax:Wax, $
  fpx:dblarr(npts), fpy:dblarr(npts), fpz:dblarr(npts), fw:dblarr(Wnpts)}

; RUN-TIME DENSITY STRUCTURE; NEEDS TO MATCH Box_Str
nx = 320L
nz = 60L
dx = 0.05D*Re
yminD = -Re
ymaxD = Re
DenMap = dblarr(nx,nz)
Den_Str = {dx:dx, DenMap:DenMap, ymin:yminD, ymax:ymaxD, nx:nx, nz:nz}

; PARTICLE TRACKING
ntrack = 2000L
ntimes = 5000L
Track_Step = 10L
Track_Str = {itrk:0L, ntrack:ntrack, ntimes:ntimes, Track_Step:Track_Step, $
  xtrk:dblarr(ntimes, ntrack), ytrk:dblarr(ntimes, ntrack), $
  ztrk:dblarr(ntimes, ntrack), Btrk:dblarr(ntimes, ntrack), $
  Wtrk:dblarr(ntimes, ntrack)}

; SETUP PARTICLE ARRRAYS
npart = 203100;203100L;203100;24000L;203100L;24000L;203100L;24000L; 203100L ; 10000L; MAXIMUM NUMBER OF ACTIVE PARTICLES ;34000bob ;203100L long
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
nload = 3000L; 60000L; 4800L; 12000L;
nstart = npart-nload
load_particles3D, nload, vth, DenZ0, DenInjZ, nEtimes, nBtimes, $
  Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, $
  x=x, y=y, z=z, px=px, py=py, pz=pz, gam=gam, dsE=dsE, dsB=dsB, nstart=nstart, ind=ind

; CONTROLL PARAMETERS
nnew = 2; 25; 2; 5;    ; NUMBER OF INJECTED PARTICLES PER TIME STEP

; DIAGNOSTICS REFRESH (IF NEEDED)
plot_mod     = 0
plot_step    = 25
record_temp  = 1
record_flux  = 0
record_dist  = 0
record_den   = 0
record_track = 1


; KNOBS
imax = 100000L ;100000L    ; NUMBER OF TIME STEPS
Eplevel =0;1.; 1.0;
Bplevel =0;1.; 1.0;
Ey0 = 2.7D; 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0
B_str.BH0 = 0D; 15d-9
B_str.ZH0 = Re

;WHEATHER PARTICLES ARE TACKED

tracki = 0
temps = make_array(imax+1)
radiation = make_array(npart, 7)



; MAIN LOOP
ree_2dplus_ion_sim, dt, it, jt, kt, npart, imax, nnew, nh, box_str, $
  B_str, Inj_str, x, y, z, px, py, pz, $
  Ex, Ey, Ez, Ey0, Epll, E1, E2, dsE, nEtimes, Eplevel, $
  Bx, By, Bz, Bmag, dBx, dBy, dBz, dB1, dB2, dB3, dsB, nBtimes, Bplevel, $
  pxn, pyn, pzn, pxold, pyold, pzold, gam, rat, plot_step, plot_mod, $
  Xpwr, Ypwr, Zpwr, psym, vth, Pdx, ind=ind, yplotrange = yplotrange, $
  record_flux=record_flux, Flux_str=Flux_Str, $
  record_temp=record_temp, Temp_Str=Temp_Str, $
  record_dist=record_dist, Dist_Str=Dist_Str, $
  record_den=record_den, Den_Str=Den_Str, $
   Track_Str=Track_Str $
  ,radiation=radiation ,maskzero=maskzero, zerocount=zerocount, mz3=mz3, mz4=mz4, maskleave=maskleave, wpi=wpi,$
  trajectories=trajectories, highp = highp, maskhot = maskhot, masklowp = masklowp, vinicial = vinicial, temps=temps $
  , Import_index = 1;,record_track=record_track,record_dist=record_dist



; INDICATOR OF HOW MANY ORIGONAL PARTICLES REMAIN
help, where(ind GE nstart)
ree_print_temps, x, y, z, px, py, pz, Bx, By, Bz, ind, Box_str

; PLOT TEMP PROFILE
t = Temp_Str.t(0:Temp_Str.kTemp-1)
TallT = Temp_Str.TturbAve(0:Temp_Str.kTemp-1)
TparT = Temp_Str.TturbPar(0:Temp_Str.kTemp-1)
TprpT = Temp_Str.TturbPerp(0:Temp_Str.kTemp-1)
ns = 1
plot, t, smooth(TallT,ns, /edge_t), yran = [0, 30.0]
oplot, t, smooth(Tpart,ns, /edge_t), col=2
oplot, t, smooth(TprpT,ns, /edge_t), col=254



; PLOT DISTRIBUTION
wi, 2
ree_distribution_basic, x, y, z, px, py, pz, ind, Box_str, /plot, $
  fpx=fpx, fpy=fpy, fpz=fpz, pax=pax, fave=fave, Wax=Wax, /turb, $
  fabspx=fabspx, fabspy=fabspy, fabspz=fabspz, pabs=pabs
Ti = 6.00d3
oplot, wax, max(fave)*exp(-Wax/Ti), col=2
wi,0

if n_elements(faveallt) eq 0 then faveallt = fave else faveAllT = faveAllT+fave

; PLOT FLUX
;flxAx = dindgen(nflx)*dw



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

; CHECK B
plot, x(ind), Bz(ind), psym=3
plot, z(ind), Bx(ind), psym=3
Bmag = sqrt(Bz*Bz + Bx*Bx)
plot, x(ind), Bmag(ind), psym=3
plot, z(ind), Bmag(ind), psym=3

ree_density_basic, x, y, z, px, py, pz, ind, Box_str, /plot

;plot, z(ind), Bx(ind), psym=3
;plot, x(ind), Bz(ind), psym=3

flux = dsega_energy_flux(radiation, box_str)
maskload = boolarr(npart)
maskload[where(radiation[*, 4]^2*mi/(ee*2) + radiation[*, 5]^2*mi/(ee*2)+ radiation[*, 3]^2*mi/(ee*2) gt 2*mean(radiation[*, 4]^2*mi/(ee*2) + radiation[*, 5]^2*mi/(ee*2)+ radiation[*, 3]^2*mi/(ee*2)))]=1


iyplus  = where(radiation[*, 1] GT yexitmax, nyp) ;if there's a provblem the problem is at the level of how radiatio is defined
datas = alog10(radiation[iyplus, 4]^2*mi/ee)
binsize = (3.5 * StdDev(datas)) / N_Elements(datas)^(0.3333)

maskleavey =  boolarr(npart)
maskleavey[iyplus] = 1

pyleave =radiation[iyplus, 4]^2*mi/(ee*2) + radiation[iyplus, 5]^2*mi/(ee*2)
highpy = pyleave[where(pyleave GT mean(pyleave))]
lowpy = pyleave[where(pyleave LT 2*mean(pyleave))]

;locations = linspace(0,7,100)


wi,1, wsize = [1700,1080]
;cgHistoPlot, alog10(radiation[iyplus, 4]^2*mi/(ee*2)), log=0,outline=1,binsize=binsize,mininput=0,maxinput=7
;cgHistoPlot, alog10(radiation[where(maskleavey*maskzero), 4]^2*mi/(ee*2)),/oplot,color='green',binsize=binsize,mininput=0,maxinput=7,outline=1

;cgHistoPlot, alog10(radiation[iyplus, 4]), log=0,outline=1,binsize=binsize,mininput=0,maxinput=7
;cgHistoPlot, alog10(radiation[where(maskleavey*maskzero), 4]),/oplot,color='green',binsize=binsize,mininput=0,maxinput=7,outline=1

cgHistoPlot, alog10(radiation[iyplus, 4]^2*mi/(ee*2) + radiation[iyplus, 5]^2*mi/(ee*2)+ radiation[iyplus, 3]^2*mi/(ee*2)),outline=1,binsize=binsize,mininput=0,maxinput=7,thick=3,xtitle = 'Log(Energy) [eV]',ytitle='Counts',charsize=3,color='red';,log=1
cgHistoPlot, alog10(radiation[where(maskleavey*maskzero), 4]^2*mi/(ee*2)+radiation[where(maskleavey*maskzero), 5]^2*mi/(ee*2)+radiation[where(maskleavey*maskzero), 3]^2*mi/(ee*2)),/oplot,color='green',binsize=binsize,mininput=0,maxinput=7,outline=1,thick=3;,log=1

cgoplot, alog10(findgen(1000000,start=1400)),2600*exp(-((findgen(1000000,start=1400)-4000)/4000)^2),color=cgcolor('blue'),thick=2
;cgoplot, alog10(findgen(1000000,start=3000)),2100*exp(-((findgen(1000000,start=3000)-12000)/25000)^2),color=cgcolor('blue'),thick=2


colors=['red','green','white','blue'];
items =['Particles leaving y+ boundary', 'Particles crossing z = 0','and leaving y+ boundary','4 keV Maxwellian']
linestyle=[0, 0,-99,0]
;al_legend,items,color=colors,charsize=3.0,thick=2.3,symsize=3.3,linestyle=linestyle,linsize=.1,box=0,/top,/left
al_legend,items,color=colors,charsize=3.0,thick=2.3,symsize=3.3,linestyle=linestyle,linsize=.1,box=0,position=[.002,2750]


save, trajectories, radiation,maskzero,maskload,maskleave,maskleavey,import_index,filename='trajectorieslongleavee0b0.sav'

; #########################   MAKE FLUX PLOT  #####################
faveAllB = fave
faveAllT = faveAllt

FluxT = (faveAllT>0.01)*Wax
print, max(faveAllT)
FluxT = fluxT/max(FluxT)

FluxB = (faveAllB>0.01)*Wax
print, max(faveAllB)
FluxB = fluxB/max(FluxB)

wi, 2
plot, Wax/1000, FluxT, /xlog, /ylog, xran=[10,1e4], yran=[1e-6, 1e2]
;oplot, Wax/1000, FluxB,col=2
wi, 0

wi, 2
Ti = 5.5d3
WaxHigh = Wax(72:*)
norm = 1.0/max(faveAllT)
plot, Wax/1000, norm*(float(faveAllT)>0.5), /ylog, /xlog,xran=[1,1e3], $
  xtit='Energy (keV)', ytit = 'Distribution (noramlized)'
oplot, Wax/1000, 0.97*norm*max(faveAllT)*exp(-Wax/Ti), col=2
oplot, WaxHigh/1000, 4.8e3*norm*((WaxHigh/1e5)^(-3.85)), col=4
wi,0


wi, 2
Ti = 5.25d3
plot, Wax/1000, float(faveAllB)>0.5, /ylog, /xlog,xran=[1,1e3], $
  xtit='Energy (keV)', ytit = 'Flux (noramlized)'
oplot, Wax/1000, 0.95*max(faveAllB)*exp(-Wax/Ti), col=2
oplot, Wax/1000, 3e3*((Wax/1e5)^(-4.0)), col=4
wi,0


; MAKE POSTSCRIPT
pname = 'DistI_Simulation'
dir = '/Users/ree/Bob/papers/20_IonAcceleration/SimPlots'
;popen, pname, /port, dir=dir, xsize = 3, ysize=4

!p.charsize=0.75

ree_sym, 0.25
Ti = 6.00d3
WaxHigh = Wax(74:*)
norm = 1.0/max(faveAllT)
plot, Wax/1000, norm*(float(faveAllT)>0.25), /ylog, /xlog,xran=[0.1,1e3], $
  xtit='Energy (keV)', ytit = 'Distribution (normalized)', psym=8, $
  title = 'Test-Particle Simulation'
oplot, Wax/1000, 0.97*exp(-Wax/Ti), col=85, thick=5
oplot, WaxHigh/1000, 4.8e3*norm*((WaxHigh/1e5)^(-3.80)), col=125, thick=5
oplot, Wax/1000, norm*(float(faveAllT)>0.25), psym=8

;pclose
!p.charsize=1.5




dumyyy  = -0.5D*cos( ( ((y(ind)+yturbmax)/Pdx) < (!dpi) ) > 0D ) - $
  0.5D*cos( ( ((yturbmax-y(ind))/Pdx) < (!dpi) ) > 0D ) & $

  plot, y(ind)/RE, dumyyy, psym=3

end

