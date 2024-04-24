; PARTICLES IN A MAGNETIC RECONNECTION REGION

; SETUP: LINK TO THE SPEDAS LIBRARY SETS PATH TO MMS DATA
;setenv, 'ROOT_DATA_DIR=/Users/ree/mmsdata/'
;mms_init, local_data_dir = '/Users/ree/mmsdata'

; COMPILE RAND_E 
; #### RESET PATH #####  
;.r '/Users/ree/Bob/papers/20_IonAcceleration/idl/ree_2dplus_lib_RevA'
;.r '/Users/ree/Bob/papers/20_IonAcceleration/idl/ree_2dplus_ion_sim_RevA'
;.r '/Users/ree/Bob/papers/20_IonAcceleration/idl/ree_2dplus_diagnostics_RevA'

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
!p.multi = [0,1,1]

; SET UP E and T
TempStart = 4000D ; eV
vth       = sqrt(TempStart*ee/mi)
pth       = sqrt(TempStart*ee*mi)

; SET UP SIMULATION DOMAIN
@2dplus_setup05A.txt

; ****** RESET POINT *******
it  = 0L ; TRACKS POSITION IN ELECTRIC FIELD ARRAYS
kt  = 0L ; TRACKS POSITION IN MAGNETIC FIELD ARRAYS
jt  = 0L ; TRACKS POSITION IN PARTICLE INJECTION ARRAYS
k   = 0L ; TRACKS DIAGNOSTIC ARRAYS


; SETUP PARTICLE ARRRAYS
npart = 2400000L; MAXIMUM NUMBER OF ACTIVE PARTICLES
;npart = 240000L; MAXIMUM NUMBER OF ACTIVE PARTICLES

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
Bmag = dblarr(npart) ; NOT NEEDED SO FARs
x    = dblarr(npart)
y    = dblarr(npart)
z    = dblarr(npart) + 900
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
nload = 1L
nstart = npart-nload
load_particles3D, nload, vth, DenZ0, DenInjZ, nEtimes, nBtimes, $
  Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, $
  x=x, y=y, z=z, px=px, py=py, pz=pz, gam=gam, dsE=dsE, dsB=dsB, nstart=nstart, ind=ind  

; CONTROLL PARAMETERS
nnew = 100l       ; NUMBER OF INJECTED PARTICLES PER TIME STEP
nnew = 10l       ; NUMBER OF INJECTED PARTICLES PER TIME STEP


; DIAGNOSTICS REFRESH (IF NEEDED)
plot_mod    = 0
plot_step   = 10
make_density = 0
make_image   = 0
make_dist    = 0

; KNOBS
imax = 500000L    ; NUMBER OF TIME STEPS
;imax = 50000L    ; NUMBER OF TIME STEPS

Eplevel = 1.
Bplevel = 1.
Ey0 = 2.7D; 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 
B_str.BH0 = 15d-9
B_str.ZH0 = Re

ind=1

;Wether particles are tracked
tracki = 0
temps = make_array(imax+1)

radiation = make_array(npart, 7)
; MAIN LOOP
ree_2dplus_ion_sim, dt, it, jt, kt, npart, imax, nnew, nh, box_str, $
  B_str, Inj_str, x, y, z, px, py, pz, $
  Ex, Ey, Ez, Ey0, Epll, E1, E2, dsE, nEtimes, Eplevel, $
  Bx, By, Bz, Bmag, dBx, dBy, dBz, dB1, dB2, dB3, dsB, nBtimes, Bplevel, $
  pxn, pyn, pzn, pxold, pyold, pzold, gam, rat, plot_step, plot_mod, $
  Xpwr, Ypwr, Zpwr, psym, vth, Pdx, ind=ind, radiation=radiation,$ 
  maskzero=maskzero, zerocount=zerocount, mz3=mz3, mz4=mz4, maskleave=maskleave, wpi=wpi,$
  trajectories=trajectories, highp = highp, maskhot = maskhot, masklowp = masklowp, vinicial = vinicial, track = track, temps=temps,npart=npart;, Import_index = 1

; INDICATOR OF HOW MANY ORIGONAL PARTICLES REMAIN
help, where(ind GE nstart)    
ree_print_temps, x, y, z, px, py, pz, Bx, By, Bz, ind, Box_str

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

; Y CURRNET
ree_y_current, z, py, gam, ind, Box_Str, /plot

; CHECK B
;plot, x(ind), Bz(ind), psym=3
;plot, z(ind), Bx(ind), psym=3
;Bmag = sqrt(Bz*Bz + Bx*Bx)
;plot, x(ind), Bmag(ind), psym=3
;plot, z(ind), Bmag(ind), psym=3

ree_density_basic, x, y, z, px, py, pz, ind, Box_str, /plot

;plot, z(ind), Bx(ind), psym=3
;plot, x(ind), Bz(ind), psym=3

flux = dsega_energy_flux(radiation, box_str)
maskload = boolarr(npart)
maskload[nstart:-1] = 1

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

cgHistoPlot, alog10(radiation[iyplus, 4]^2*mi/(ee*2) + radiation[iyplus, 5]^2*mi/(ee*2)+ radiation[iyplus, 3]^2*mi/(ee*2)), log=0,outline=1,binsize=binsize,mininput=0,maxinput=7,thick=3,xtitle = 'Log(Energy) [eV]',ytitle='Counts',charsize=3,color='red'
cgHistoPlot, alog10(radiation[where(maskleavey*maskzero), 4]^2*mi/(ee*2)+radiation[where(maskleavey*maskzero), 5]^2*mi/(ee*2)+radiation[where(maskleavey*maskzero), 3]^2*mi/(ee*2)),/oplot,color='green',binsize=binsize,mininput=0,maxinput=7,outline=1,thick=3

colors=['red','green'];
items =['Particles leaving y+ boundary', 'Particles crossing z = 0','and leaving y+ boundary']
linestyle=[0, 0,-99]

al_legend,items,color=colors,charsize=3.0,thick=2.3,symsize=3.3,/top,/left,linestyle=linestyle,linsize=.1
;cglegend, Titles = ['Particles leave', 'Particles leave & cross z=0'],color=['orange','green'],charsize=2.3,location=[.13,.8],thick=3,/box
;cgHistoPlot, radiation[where(maskleave*maskzero*mz4), 4],/oplot,McKee, Christopher F, and Ostriker, Eve C. ?Theory of Star Formation,? in  Annual Review of Astronomy and Astrophysics, vol. 45.  Palo Alto: Annual Reviews, 2007: 565?687color='cyan',binsize=binsize
;cgHistoPlot, radiation[where(maskleave*masklowp), 4],/oplot,color='violet',binsize=binsize, /log
;cgHistoPlot, alog10(highpy^2*mi/ee),/oplot,color='blue',binsize=binsize,mininput=0,maxinput=7 
;cgHistoPlot, alog10(radiation[where(maskleavey*maskload), 4]^2*mi/ee),/oplot,color='cyan',binsize=binsize,mininput=0,maxinput=7,outline=1


write_png, 'AllKE2.png',tvrd(/true)

cglegend, Titles =['Imax =' + string(imax), 'nnew =' + string(nnew), 'nload =' + string(nload)]

match, iyplus, where(zerocount GT 1, nzc), suba,subb
;cgHistoPlot, radiation[iyplus[suba], 4],/oplot,color='cyan',binsize=binsize, /log

match, highpy,  radiation[*, 4]^2+radiation[*, 5]^2+radiation[*,3]^2,suba,subb
 
maskload[*] = 0
maskload[subb] = 1











iyminus  = where(radiation[*, 1] LT yexitmin, nzp)
datas = radiation[iyminus, 4]
binsize = (3.5 * StdDev(datas)) / N_Elements(datas)^(0.3333)

maskleave1 = boolarr(npart)
maskleave1[*] = 0
maskleave1[iyminus] = 1

save,trajectories,maskzero,maskleavey,box_str,filename='trajectories.sav'

;wi, 1

;radiation[where[radiation[iymins,4] GT mean(radiation[iyminys,4]) , 4],binsize=binsize
;cgHistoPlot, radiation[iyminus, 4],binsize=binsize
;cgHistoPlot, radiation[where(maskleave1*maskzero eq 1), 4],/oplot,color='green',binsize=binsize
;cgHistoPlot, radiation[where(maskleave1*maskzero*mz4 eq 1), 4],/oplot,color='blue',binsize=binsize
;cgHistoPlot, radiation[where(zerocount GT 10), 5],/oplot,color='violet',binsize=binsize,/window
;cglegend, Titles =['Imax =' + string(imax), 'inew =' + string(nnew), 'nload =' + string(nload)]

;x=findgen(2000,increment=100000)
;cgoplot,x,2000*exp(-x/1000000)

print,'particles leaving at y+', n_elements(flux.yplus)
print,'particles leaving at y-', n_elements(flux.yminus)
print,'particles leaving at x+', n_elements(flux.xplus)
print,'particles leaving at x-', n_elements(flux.xminus)
print,'particles leaving at z+', n_elements(flux.zplus)
print,'particles leaving at z-', n_elements(flux.zminus)

if tracki then begin

goodx = where(trajectories[*, 0, 0] NE 0)
goody = where(trajectories[*, 1, 0] NE 0)
px =  trajectories[goodx, 3, *]
py =  trajectories[goodx, 4, *]
y  = trajectories[goody, 1, *]

;restore, 'import_index_e01.sav
restore, 'import_indext.sav'

fines = dindgen(n_elements(import_index))

steps = 20

x  = make_array(uint(imax/steps) + 1,n_elements(import_index), value = !VALUES.F_NAN )
y  = make_array(uint(imax/steps) + 1,n_elements(import_index), value = !VALUES.F_NAN)
z  = make_array(uint(imax/steps) + 1,n_elements(import_index), value = !VALUES.F_NAN)
px = make_array(uint(imax/steps) + 1,n_elements(import_index), value = !VALUES.F_NAN)
py = make_array(uint(imax/steps) + 1,n_elements(import_index), value = !VALUES.F_NAN)
pz = make_array(uint(imax/steps) + 1,n_elements(import_index), value = !VALUES.F_NAN)

Ex = make_array(uint(imax/steps) + 1,n_elements(import_index), value = !VALUES.F_NAN)
Ey = make_array(uint(imax/steps) + 1,n_elements(import_index), value = !VALUES.F_NAN)
Ez = make_array(uint(imax/steps) + 1,n_elements(import_index), value = !VALUES.F_NAN)
Bx = make_array(uint(imax/steps) + 1,n_elements(import_index), value = !VALUES.F_NAN)
By = make_array(uint(imax/steps) + 1,n_elements(import_index), value = !VALUES.F_NAN)
Bz = make_array(uint(imax/steps) + 1,n_elements(import_index), value = !VALUES.F_NAN)


for j = 0, n_elements(import_index)-1 do begin
  
  goodx = where(trajectories[*, 0, j] NE 0)
  goody = where(trajectories[*, 1, j] NE 0)
  goodz = where(trajectories[*, 2, j] NE 900)
  dpx = trajectories[goodx,3,j]-shift(trajectories[goodx,3,j],1)
  fin = where(dpx eq 0)
  if n_elements(fin) gt 1 then fin = fin[1]
  if fin eq -1 then fin=n_elements(dpx)-1
   
   
   
    x[0:fin,j] = trajectories[goodx[0:fin], 0, j]
    y[0:fin,j] = trajectories[goody[0:fin], 1, j]
    z[0:fin,j] = trajectories[goodz[0:fin], 2, j]
   px[0:fin,j] = trajectories[goodx[0:fin], 3, j]
   py[0:fin,j] = trajectories[goody[0:fin], 4, j]
   pz[0:fin,j] = trajectories[goodz[0:fin], 5, j]

   Ex[0:fin,j] = trajectories[goodx[0:fin], 6, j]
   Ey[0:fin,j] = trajectories[goody[0:fin], 7, j]
   Ez[0:fin,j] = trajectories[goodz[0:fin], 8, j]
   Bx[0:fin,j] = trajectories[goodx[0:fin], 9, j]
   By[0:fin,j] = trajectories[goody[0:fin], 10,j]
   Bz[0:fin,j] = trajectories[goodz[0:fin], 11,j]    

  fines[j] = fin
endfor

energy0 = mean(px[0,*]^2 + py[0:*]^2 + pz[0:*]^2,/nan)
energy1 = mean(diag_matrix(px[fines,*])^2 + diag_matrix(py[fines,*])^2 + diag_matrix(pz[fines,*])^2)

endif else begin

px1 = vinicial[0,*]^2*mi/ee
py1 = vinicial[1,*]^2*mi/ee
pz1 = vinicial[2,*]^2*mi/ee

px2 = radiation[import_index, 3]^2*mi/ee
py2 = radiation[import_index, 4]^2*mi/ee
pz2  =radiation[import_index, 5]^2*mi/ee

plot, px1 + py1 + pz1,( px2 + py2 + pz2 - (px1 + py1 + pz1))/(px1 + py1 + pz1),psym=3,yrange=[0,5]

endelse

import_index = where(maskleavey*maskzero,n)
save, import_index, filename = 'import_indexb.sav'

;import_index = where(masklowp)
;save,import_index,filename='import_indexslow.sav'

;import_index = where(maskhot)
;save,import_index,filename='import_indexfast.sav'





;izplus  = where(radiation[*, 2] GT zexitmax, nzp)
;datas = radiation[izplus, 5]
;binsize = (3.5 * StdDev(datas)) / N_Elements(datas)^(0.3333)

;maskleave[*] = 0
;maskleave[iyminus] = 1

;wi, 1

;cgHistoPlot, radiation[iyminus, 5],binsize=binsize
;cgHistoPlot, radiation[where(maskleave*maskzero eq 1), 5],/oplot,color='green',binsize=binsize
;cgHistoPlot, radiation[where(maskleave*maskzero*mz3 eq 1), 5],/oplot,color='blue',binsize=binsize
;;cgHistoPlot, radiation[where(zerocount GT 10), 5],/oplot,color='violet',binsize=binsize,/window
;cglegend, Titles =['Imax =' + string(imax), 'inew =' + string(nnew), 'nload =' + string(nload)]



end
