pro plot_yflux

restore, 'trajectorieslongleavee0b0.sav'


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
pth       = sqrt(TempStart*ee*mi)

xmin     = -8.0D*Re ; LOCATION OF INJECTION BOUNDARY (PARTICLE ENTRY)
xmax     =  8.0D*Re ;
xexitmin = -8.0D*Re ; LOCATION PARTICLE EXIT BOUNDARY
xexitmax =  8.0D*Re
xturbmin = -7.75D*Re ; LOCATION OF TURBULENCE REGION
xturbmax =  7.75D*Re

zmin     = -1.5D*Re  ;
zmax     =  1.5D*Re
zexitmin =  zmin
zexitmax =  zmax
zturbmin = -1.25D*Re  ;
zturbmax =  1.25D*Re

ymin     = -1.5*Re
ymax     =  1.5*Re
yexitmin = ymin
yexitmax = ymax
yturbmin = -1.25*Re
yturbmax =  1.25*Re

; BOX STRUCTURE CONTAINS NEEDED INFORMATION OF THE BOX
box_str = {xmin:xmin, xmax:xmax, xexitmin:xexitmin, xexitmax:xexitmax, $
  xturbmin:xturbmin, xturbmax:xturbmax, zmin:zmin, zmax:zmax, $
  zexitmin:zexitmin, zexitmax:zexitmax, zturbmin:zturbmin, $
  zturbmax:zturbmax, ymin:ymin, ymax:ymax, yexitmin:yexitmin, $
  yexitmax:yexitmax, yturbmin:yturbmin, yturbmax:yturbmax}
npart = 203100
!p.multi=[0,2,1]

wi,1, wsize = [1950,1000]


  flux = dsega_energy_flux(radiation, box_str)

  iyplus  = where(radiation[*, 1] GT yexitmax, nyp) ;if there's a provblem the problem is at the level of how radiatio is defined
  datas = alog10(radiation[iyplus, 4]^2*mi/ee)
  binsize = (3.5 * StdDev(datas)) / N_Elements(datas)^(0.3333)

  maskleavey =  boolarr(npart)
  maskleavey[iyplus] = 1

  pyleave =radiation[iyplus, 4]^2*mi/(ee*2) + radiation[iyplus, 5]^2*mi/(ee*2)
  highpy = pyleave[where(pyleave GT mean(pyleave))]
  lowpy = pyleave[where(pyleave LT 2*mean(pyleave))]
  
  cgHistoPlot, alog10(radiation[iyplus, 4]^2*mi/(ee*2) + radiation[iyplus, 5]^2*mi/(ee*2)+ radiation[iyplus, 3]^2*mi/(ee*2)),outline=1,binsize=binsize,mininput=0,maxinput=7,thick=3,xtitle = 'Log(Energy) [eV]',charsize=3.6,color='red',xmargin=[4.8,.15];, YTICKFORMAT="(A1)";,log=1
  cgHistoPlot, alog10(radiation[where(maskleavey*maskzero), 4]^2*mi/(ee*2)+radiation[where(maskleavey*maskzero), 5]^2*mi/(ee*2)+radiation[where(maskleavey*maskzero), 3]^2*mi/(ee*2)),/oplot,color='green',binsize=binsize,mininput=0,maxinput=7,outline=1,thick=3;,log=1

  cgoplot, alog10(findgen(1000000,start=1400)),2200*exp(-((findgen(1000000,start=1400)-4000)/4000)^2),color=cgcolor('blue'),thick=2
  ;cgoplot, alog10(findgen(1000000,start=3000)),2100*exp(-((findgen(1000000,start=3000)-12000)/25000)^2),color=cgcolor('blue'),thick=2

  colors=['red','green','white','blue'];
  items =['Ion exits y+ face', 'Ion crosses z=0','and exits y+ face','4 keV Maxwellian']
  linestyle=[0, 0,-99,0]
  ;al_legend,items,color=colors,charsize=3.0,thick=2.3,symsize=3.3,linestyle=linestyle,linsize=.1,box=0,/top,/left
  al_legend,items,color=colors,charsize=3.1,thick=2.3,symsize=3.3,linestyle=linestyle,linsize=.1,box=0,position=[-0.2,2300]
  
  restore, 'trajectorieslongleave.sav'

  flux = dsega_energy_flux(radiation, box_str)

  iyplus  = where(radiation[*, 1] GT yexitmax, nyp) ;if there's a provblem the problem is at the level of how radiatio is defined
  datas = alog10(radiation[iyplus, 4]^2*mi/ee)
  binsize = (3.5 * StdDev(datas)) / N_Elements(datas)^(0.3333)

  maskleavey =  boolarr(npart)
  maskleavey[iyplus] = 1

  pyleave =radiation[iyplus, 4]^2*mi/(ee*2) + radiation[iyplus, 5]^2*mi/(ee*2)
  highpy = pyleave[where(pyleave GT mean(pyleave))]
  lowpy = pyleave[where(pyleave LT 2*mean(pyleave))]

  ;locations = linspace(0,7,100)



  ;cgHistoPlot, alog10(radiation[iyplus, 4]^2*mi/(ee*2)), log=0,outline=1,binsize=binsize,mininput=0,maxinput=7
  ;cgHistoPlot, alog10(radiation[where(maskleavey*maskzero), 4]^2*mi/(ee*2)),/oplot,color='green',binsize=binsize,mininput=0,maxinput=7,outline=1

  ;cgHistoPlot, alog10(radiation[iyplus, 4]), log=0,outline=1,binsize=binsize,mininput=0,maxinput=7
  ;cgHistoPlot, alog10(radiation[where(maskleavey*maskzero), 4]),/oplot,color='green',binsize=binsize,mininput=0,maxinput=7,outline=1

  cgHistoPlot, alog10(radiation[iyplus, 4]^2*mi/(ee*2) + radiation[iyplus, 5]^2*mi/(ee*2)+ radiation[iyplus, 3]^2*mi/(ee*2)),outline=1,binsize=binsize,mininput=0,maxinput=7,thick=3,xtitle = 'Log(Energy) [eV]',ytitle='Counts',charsize=3.6,color='red',xmargin=[4.8,1];,log=1
  cgHistoPlot, alog10(radiation[where(maskleavey*maskzero), 4]^2*mi/(ee*2)+radiation[where(maskleavey*maskzero), 5]^2*mi/(ee*2)+radiation[where(maskleavey*maskzero), 3]^2*mi/(ee*2)),/oplot,color='green',binsize=binsize,mininput=0,maxinput=7,outline=1,thick=3;,log=1

  cgoplot, alog10(findgen(1000000,start=1400)),2600*exp(-((findgen(1000000,start=1400)-4000)/4000)^2),color=cgcolor('blue'),thick=2
  ;cgoplot, alog10(findgen(1000000,start=3000)),2100*exp(-((findgen(1000000,start=3000)-12000)/25000)^2),color=cgcolor('blue'),thick=2


  ; colors=['red','green','white','blue'];
  ; items =['Particles leaving y+ boundary', 'Particles crossing z = 0','and leaving y+ boundary','4 keV Maxwellian']
  ; linestyle=[0, 0,-99,0]
  ; ;al_legend,items,color=colors,charsize=3.0,thick=2.3,symsize=3.3,linestyle=linestyle,linsize=.1,box=0,/top,/left
  ; al_legend,items,color=colors,charsize=3.5,thick=2.3,symsize=3.3,linestyle=linestyle,linsize=.1,box=0,position=[.0015,2750]
  
  write_png, 'allke.png',tvrd(/true)


  end