; PARTICLES IN A MAGNETIC RECONNECTION REGION SMALLER Y DIRECTION
undefine,off
;off=1
;pro run_2dplus08
; SETUP: LINK TO THE SPEDAS LIBRARY SETS PATH TO MMS DATA
;setenv, 'ROOT_DATA_DIR=/Users/ree/mmsdata/'
;mms_init, local_data_dir = '/Users/ree/mmsdata'

; COMPILE RAND_E 
; #### RESET PATH #####  
;'C:\Users\super\IDLWorkspaceBob\Default
;.r 'C:\Users\super\IDLWorkspaceBob\Default'

;.r '/Users/super/Bob/papers/22_Daniel/idl/ree_2dplus_lib_B'
;.r '/Users/ree/Bob/papers/22_Daniel/idl/ree_2dplus_ion_sim_B'
;.r '/Users/ree/Bob/papers/22_Daniel/idl/ree_2dplus_diagnostics_B'
;.r '/Users/ree/Bob/papers/22_Daniel/idl/ree_2dplus_ion_sim_waves'

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
@2dplus_setup08.txt

; ****** RESET POINT *******
it  = 0L ; TRACKS POSITION IN TIME
jt  = 0L ; TRACKS POSITION IN PARTICLE INJECTION ARRAYS
k   = 0L ; TRACKS DIAGNOSTIC ARRAYS

; RUN-TIME TEMPERATURE PROFILE
Ntemps = 20000L
ktemp  = 0L
Buff   = 0.0D
Temp_Step = 25L
Temp_Str = {Ntemps:Ntemps, kTemp:kTemp, Temp_Step:Temp_Step, Buff:Buff, $
 t:dblarr(Ntemps), TboxPar:dblarr(Ntemps), TboxPerp:dblarr(Ntemps), $
 TboxAve:dblarr(Ntemps), TturbPar:dblarr(Ntemps), TturbPerp:dblarr(Ntemps), $
 TturbAve:dblarr(Ntemps)}

; RUN-TIME FLUX STRUCTURE
Nflux = 1000L
dW = 500D
Ntimes = 100000L
Nebins = 12L
Flux_Str = {dW:dw, Nflux:Nflux, px:dblarr(Nflux), mx:dblarr(Nflux), $
 py:dblarr(Nflux), my:dblarr(Nflux), pz:dblarr(Nflux), mz:dblarr(Nflux), $
 FlxVtime:dblarr(Ntimes,Nebins), Nexit:lonarr(Ntimes), it:0L}

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
BuffSize = 0.25D
Ystop   = 1.00D
Ystart  = -0.00D

Dist_Str = {Pmax:Pmax, Bin:Bin, npts:npts, Turb:Turb, pax:pax, $
  BuffSize:BuffSize, Ystart:Ystart, Ystop:Ystop, $
  Dist_Step:Dist_Step, Wmax:Wmax, Wbin:Wbin, Wnpts:Wnpts, Wax:Wax, $
  fpx:dblarr(npts), fpy:dblarr(npts), fpz:dblarr(npts), fw:dblarr(Wnpts)}

; RUN-TIME DENSITY STRUCTURE; NEEDS TO MATCH Box_Str
nx = 320L
nz = 80L
dx = 0.05D*Re
yminD = -Re
ymaxD = Re
DenMap = dblarr(nx,nz)
Den_Str = {dx:dx, DenMap:DenMap, ymin:yminD, ymax:ymaxD, nx:nx, nz:nz}

; PARTICLE TRACKING
ntrack = 100L; 2000L
ntimes = 3000L
Track_Step = 10L
Track_Str = {itrk:0L, ntrack:ntrack, ntimes:ntimes, Track_Step:Track_Step, $
             xtrk:dblarr(ntimes, ntrack), ytrk:dblarr(ntimes, ntrack), $
             ztrk:dblarr(ntimes, ntrack), Btrk:dblarr(ntimes, ntrack), $
             Wtrk:dblarr(ntimes, ntrack)}

; SETUP PARTICLE ARRRAYS
npart = 1000001L; 10000L; MAXIMUM NUMBER OF ACTIVE PARTICLES
Ex   = dblarr(npart)
Ey   = dblarr(npart)
Ez   = dblarr(npart)
Xpwr = dblarr(npart)
Ypwr = dblarr(npart)
Zpwr = dblarr(npart)
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
dsE  = lonarr(npart) ; NOT USED, BUT NEEDED
dsB  = lonarr(npart) ; NOT USED, BUT NEEDED

; LOAD PARTICLES
nload = 1L; 60000L; 4800L; 12000L;
nstart = npart-nload
DenZ0    = 1*Z0  ; NEEDED TO BALANCE B PRESSURE

load_particles3D, nload, vth, DenZ0, DenInjZ, nEtimes, nBtimes, $
  Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, $
  x=x, y=y, z=z, px=px, py=py, pz=pz, gam=gam, dsE=dsE, dsB=dsB, nstart=nstart, ind=ind

; CONTROLL PARAMETERS
nnew = 10l;     ; NUMBER OF INJECTED PARTICLES PER TIME STEP
; DIAGNOSTICS REFRESH (IF NEEDED)
plot_mod     = 0
plot_step    = 0
record_temp  = 1
record_flux  = 0
record_dist  = 0
record_den   = 0
record_track = 0


; KNOBS
imax = 100000L    ; NUMBER OF TIME STEPS
Eplevel = 1.;1.;0;1.0;  
Bplevel = 1.;1.;0;1.0; 
Ey0 = 2.7D;

if keyword_set(off) then begin
  Eplevel = 0.;1.;0;1.0;
  Bplevel = 0.;1.;0;1.0;
endif

radiation = make_array(npart, 7)

; MAIN LOOP
ree_2dplus_ion_sim_waves, dt, it, jt, imax, nnew, box_str, $
  B_str, Inj_str, x, y, z, px, py, pz, $
  Ex, Ey, Ez, Ey0, dEeStr, dEeStrP, Eplevel, $
  Bx, By, Bz, Bmag, dBx, dBy, dBz, dBgStr, Bplevel, $
  pxn, pyn, pzn, pxold, pyold, pzold, gam, rat, plot_step, plot_mod, $
  Xpwr, Ypwr, Zpwr, psym, vth, Pdx, ind=ind, yplotrange = yplotrange, $
  record_flux=record_flux, Flux_str=Flux_Str, $
  record_temp=record_temp, Temp_Str=Temp_Str, $
  record_dist=record_dist, Dist_Str=Dist_Str, $
  record_den=record_den, Den_Str=Den_Str, $
  record_track=record_track, Track_Str=Track_Str , $
  radiation=radiation,npart=npart, $
  maskzero=maskzero, zerocount=zerocount, mz3=mz3, mz4=mz4, maskleave = maskleave, wpi=wpi, $
  trajectories=trajectories, highp = highp, maskhot = maskhot, masklowp = masklowp, off=off;, import_index = 1


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


; PLOT DISTRUBUTION 
;wi, 2
;ree_plot_dist, Dist_str, Tuse=4d3, Beta=1.1, Ncal=0.04d6,xran=[5e-2, 2e3], yran=[3e3,3e8], Wfit1=[10.0, 80.0], Wfit2=[100.0, 200.0],Wfpi=fpi_en_select, Ffpi=fpi_flx_select, Wfeeps=feep_en_select, Ffeeps=feep_flx_select, Wmaxwell=[0.0,40.0]
;wi, 0



; RECORD
TTot0 = TallT
Tpar0 = TparT
Tprp0 = TprpT


;oplot, t, (TTot0 + TTot1 + TTot2 + Ttot3 + Ttot4 + Ttot5 + Ttot6 + $
;  Ttot7 + Ttot8 + Ttot9)/10, col=1, thick=5

;TtotAve = (TTot0 + TTot1 + TTot2 + Ttot3 + Ttot4 + Ttot5 + Ttot6 + $
;  Ttot7 + Ttot8 + Ttot9)/10
;TprpAve = (Tprp0 + Tprp1 + Tprp2 + Tprp3 + Tprp4 + Tprp5 + Tprp6 + $
;  Tprp7 + Tprp8 + Tprp9)/10
;TparAve = (Tpar0 + Tpar1 + Tpar2 + Tpar3 + Tpar4 + Tpar5 + Tpar6 + $
;  Tpar7 + Tpar8 + Tpar9)/10
  
  flux = dsega_energy_flux(radiation, box_str)
  maskload = boolarr(npart)
  ;maskload[nstart:-1] = 1

  iyplus  = where(radiation[*, 1] GT yexitmax, nyp) ;if there's a provblem the problem is at the level of how radiatio is defined
  datas = alog10(radiation[iyplus, 4]^2*mi/ee)
  binsize = (3.5 * StdDev(datas)) / N_Elements(datas)^(0.3333)

  maskleavey =  boolarr(npart)
  maskleavey[iyplus] = 1
  
  pyleave =radiation[iyplus, 4]^2*mi/(ee*2) + radiation[iyplus, 5]^2*mi/(ee*2)
  highpy = pyleave[where(pyleave GT mean(pyleave))]
  lowpy = pyleave[where(pyleave LT 2*mean(pyleave))]

  ;locations = linspace(0,7,100)


  wi,2, wsize = [1600,1080]
  ;cgHistoPlot, alog10(radiation[iyplus, 4]^2*mi/(ee*2)), log=0,outline=1,binsize=binsize,mininput=0,maxinput=7
  ;cgHistoPlot, alog10(radiation[where(maskleavey*maskzero), 4]^2*mi/(ee*2)),/oplot,color='green',binsize=binsize,mininput=0,maxinput=7,outline=1

  ;cgHistoPlot, alog10(radiation[iyplus, 4]), log=0,outline=1,binsize=binsize,mininput=0,maxinput=7
  ;cgHistoPlot, alog10(radiation[where(maskleavey*maskzero), 4]),/oplot,color='green',binsize=binsize,mininput=0,maxinput=7,outline=1
  
  index = where(zerocount ge 2)
  maskload[index] = 1

  cgHistoPlot, alog10(radiation[iyplus, 4]^2*mi/(ee*2) + radiation[iyplus, 5]^2*mi/(ee*2)+ radiation[iyplus, 3]^2*mi/(ee*2)), log=0,outline=1,binsize=binsize,mininput=0,maxinput=7,thick=3,xtitle = 'Log(Energy) [eV]',ytitle='Counts',charsize=3,color='red'
  ;cgHistoPlot, alog10(radiation[where(maskleavey*maskload), 4]^2*mi/(ee*2)+radiation[where(maskleavey*maskload), 5]^2*mi/(ee*2)+radiation[where(maskleavey*maskload), 3]^2*mi/(ee*2)),/oplot,color='green',binsize=binsize,mininput=0,maxinput=7,outline=1,thick=3
 
 
  cgHistoPlot, alog10(radiation[where(maskleavey*maskzero), 4]^2*mi/(ee*2)+radiation[where(maskleavey*maskzero), 5]^2*mi/(ee*2)+radiation[where(maskleavey*maskzero), 3]^2*mi/(ee*2)),/oplot,color='green',binsize=binsize,mininput=0,maxinput=7,outline=1,thick=3


  colors=['red','green'];
  items =['Particles leaving y+ boundary', 'Particles crossing z = 0','and leaving y+ boundary']
  linestyle=[0, 0,-99]
  al_legend,items,color=colors,charsize=3.0,thick=2.3,symsize=3.3,/top,/left,linestyle=linestyle,linsize=.1,box=0


  print,'particles leaving at y+', n_elements(flux.yplus)
  print,'particles leaving at y-', n_elements(flux.yminus)
  print,'particles leaving at x+', n_elements(flux.xplus)
  print,'particles leaving at x-', n_elements(flux.xminus)
  print,'particles leaving at z+', n_elements(flux.zplus)
  print,'particles leaving at z-', n_elements(flux.zminus)

if keyword_set(import_index) then begin
  
  if keyword_set(off) then save,trajectories,maskzero,maskleaveybox_str,radiation,filename='trajectoriesl0.sav'else save,trajectories,maskzero,maskleavey,box_str,radiation,filename='trajectoriesl1.sav'
endif

  import_index = where(maskleave,n)
    if keyword_set(off) then  save, import_index, filename = 'import_indexwall0.sav' else  save, import_index, filename = 'import_indexwall1.sav'
     
     maskzero2 = maskload
    if keyword_set(off) then save, maskzero2 , filename='maskindexl0.sav' else save, maskzero2 , filename='maskindexl1.sav' 
  

  
  



end







