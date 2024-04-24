pro plot_dispacement2 ,maskzero=maskzero,import_index=import_index, maskleavey=maskleavey

shift=0.25
; restore, 'import_indexwall0.sav'
  ;restore, 'trajectorieslb0e0.sav'
  restore, 'trajectoriesl0.sav'


;  restore, 'trajectorieslongalle0b0.sav'
  ;restore, 'trajectoriesBOBALL.sav
  
  flagleavey=0
  
  if flagleavey then begin
    ;if keyword_set(maskzero) then match,where(maskzero),import_index,suba,subb
    match,where(maskleavey),import_index,suba,leavey
    trajectories = temporary(trajectories[*,*,leavey])
  endif
  
 
  

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
  ey0 = 2.7d
  ;example forloopcuts, trajectories,maskzero=maskzero,import_index=import_index, maskleavey=maskleavey
  ;trajectories = trajectories[*,*,leavey]

  Re    = 6374.d3
  delvar, badind


  ; SET UP BOX
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

  ymin     = -2*Re
  ymax     =  2*Re
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

  n     = size(trajectories)
  fines = make_array(n[3])
  inis  = make_array(n[3])

  x  = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN )
  y  = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  z  = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  px = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  py = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  pz = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)

  ; Ex = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  ; Ey = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  ; Ez = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  ; Bx = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  ; By = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  ; Bz = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)


  for number=0, n[3]-1 do begin


    goodx = where(trajectories[*, 3, number] NE 0)
    goody = where(trajectories[*, 4, number] NE 0)
    goodz = where(trajectories[*, 5, number] NE 0)


  dpx = trajectories[goodx,3,number]-shift(trajectories[goodx,3,number],1)
  ddpx = dpx-shift(dpx,1)
  dddpx = ddpx-shift(ddpx,1)
  ddddpx = dddpx-shift(dddpx,1)
  dddddpx = ddddpx-shift(ddddpx,1)
  fin = where(dpx eq 0 and ddpx eq 0and dddpx eq 0and ddddpx eq 0 ,m)
  
  if m gt 1 then fin = fin[0] else fin=-1
  ini=1

    if n_elements(goodx) le 1 or fin[0] le 1 then begin
      print, 'no goodx'

      if n_elements(badind) eq 0 then badind = [number] else badind = [badind,number]
      flag=1

    endif else begin


      x[goodx[ini:fin],number] = trajectories[goodx[ini:fin], 0, number]
      y[goody[ini:fin],number] = trajectories[goody[ini:fin], 1, number]
      z[goodz[ini:fin],number] = trajectories[goodz[ini:fin], 2, number]
      px[goodx[ini:fin],number] = trajectories[goodx[ini:fin], 3, number]
      py[goody[ini:fin],number] = trajectories[goody[ini:fin], 4, number]
      pz[goodz[ini:fin],number] = trajectories[goodz[ini:fin], 5, number]

      ;    Ex[goodx[ini:fin],number] = trajectories[goodx[ini:fin], 6, number]
      ;    Ey[goody[ini:fin],number] = trajectories[goody[ini:fin], 7, number]
      ;    Ez[goodz[ini:fin],number] = trajectories[goodz[ini:fin], 8, number]
      ;   Bx[goodx[ini:fin],number] = trajectories[goodx[ini:fin], 9, number]
      ;    By[goody[ini:fin],number] = trajectories[goody[ini:fin], 10, number]
      ;    Bz[goodz[ini:fin],number] = trajectories[goodz[ini:fin], 11, number]


      fines[number] = goodx[fin]
      inis[number]  = goodx[ini]

    endelse

  endfor
  radiation=0
  trajectories=0
  ;goodx=0
  ;goody=0
  ;goodz=0
  ;undefine,trajectories
  ;undefine,goodx
  ;undefine,goody
  ;undefine,goodz
  ;undefine,maskzero,maskload



  n=size(x)
  x0e = make_array(n[2])
  xfe = make_array(n[2])
  y0e = make_array(n[2])
  yfe = make_array(n[2])
  z0e = make_array(n[2])
  zfe = make_array(n[2])
  px0e = make_array(n[2])
  pxfe = make_array(n[2])
  py0e = make_array(n[2])
  pyfe = make_array(n[2])
  pz0e = make_array(n[2])
  pzfe = make_array(n[2])

  for i=0,n[2]-1 do begin

    x0e[i] = temporary(x[inis[i],i])
    y0e[i] = temporary(y[inis[i],i])
    z0e[i] = temporary(z[inis[i],i])
    px0e[i] = temporary(px[inis[i],i])
    py0e[i] = temporary(py[inis[i],i])
    pz0e[i] = temporary(pz[inis[i],i])

    xfe[i] = temporary((x[fines[i],i]))
    yfe[i] = temporary((y[fines[i],i]))
    zfe[i] = temporary((z[fines[i],i]))
    pxfe[i] = temporary((px[fines[i],i]))
    pyfe[i] = temporary((py[fines[i],i]))
    pzfe[i] = temporary((pz[fines[i],i]))

  endfor

  remove, badind, x0e, y0e, z0e, px0e, py0e, pz0e, xfe, yfe, zfe, pxfe, pyfe, pzfe,fines,inis


  energy = py0e^2*mi/(ee*2)+px0e^2*mi/(ee*2)+pz0e^2*mi/(ee*2)
; energy=0
  energyf = pyfe^2*mi/(ee*2)+pxfe^2*mi/(ee*2)+pzfe^2*mi/(ee*2)

  deltae  =  pyfe^2*mi/(ee*2)+pxfe^2*mi/(ee*2)+pzfe^2*mi/(ee*2) - (py0e^2*mi/(ee*2)+px0e^2*mi/(ee*2)+pz0e^2*mi/(ee*2))

  points_in_bin=1
  bins = n_elements(y0e)/points_in_bin
  disp=yfe-y0e
  dispes=yfe-(y0e+shift*re)

  redispe=disp[sort(disp)]
  redispes=dispes[sort(dispes)]

  redeltae=deltae[sort(disp)]

  redispe=rebin(redispe[0:bins*points_in_bin-1],bins)
  redispes=rebin(redispes[0:bins*points_in_bin-1],bins)

  redeltae=rebin(redeltae[0:bins*points_in_bin-1],bins)
  print,mean(disp)/ymax
  print,mean(fines-inis)*0.1

  ;restore, 'trajectorieslongall.sav
  restore, 'trajectoriesl1.sav'
  ;restore, 'import_indexwall1.sav'

  
  if flagleavey then begin
    if keyword_set(maskzero) then match,where(maskzero),import_index,suba,subb
    if keyword_Set(maskleavey) then match,where(maskleavey),import_index,suba,leavey
    trajectories = trajectories[*,*,leavey]

  endif


  delvar, badind

  n     = size(trajectories)
  fines = make_array(n[3])
  inis  = make_array(n[3])

  x  = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN )
  y  = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  z  = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  px = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  py = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  pz = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)

  ;  Ex = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  ;  Ey = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  ;  Ez = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  ; Bx = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  ;  By = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  ;  Bz = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)


  for number=0, n[3]-1 do begin


    goodx = where(trajectories[*, 3, number] NE 0)
    goody = where(trajectories[*, 4, number] NE 0)
    goodz = where(trajectories[*, 5, number] NE 0)

    dpx = trajectories[goodx,3,number]-shift(trajectories[goodx,3,number],1)
    fin = where(dpx eq 0)
    ;fin = fin[1]

    fin=-1
    ini=1

    if n_elements(goodx) le 1 then begin
      print, 'no goodx'

      if n_elements(badind) eq 0 then badind = [number] else badind = [badind,number]

    endif else begin


      x[goodx[ini:fin],number] = trajectories[goodx[ini:fin], 0, number]
      y[goody[ini:fin],number] = trajectories[goody[ini:fin], 1, number]
      z[goodz[ini:fin],number] = trajectories[goodz[ini:fin], 2, number]
      px[goodx[ini:fin],number] = trajectories[goodx[ini:fin], 3, number]
      py[goody[ini:fin],number] = trajectories[goody[ini:fin], 4, number]
      pz[goodz[ini:fin],number] = trajectories[goodz[ini:fin], 5, number]

      ;  Ex[goodx[ini:fin],number] = trajectories[goodx[ini:fin], 6, number]
      ;  Ey[goody[ini:fin],number] = trajectories[goody[ini:fin], 7, number]
      ;  Ez[goodz[ini:fin],number] = trajectories[goodz[ini:fin], 8, number]
      ; Bx[goodx[ini:fin],number] = trajectories[goodx[ini:fin], 9, number]
      ;  By[goody[ini:fin],number] = trajectories[goody[ini:fin], 10, number]
      ;  Bz[goodz[ini:fin],number] = trajectories[goodz[ini:fin], 11, number]


      fines[number] = goodx[fin]
      inis[number]  = goodx[ini]

    endelse

  endfor
  undefine, trajectories
  
  m=size(x)
  x0 = make_array(m[2])
  xf = make_array(m[2])
  y0 = make_array(m[2])
  yf = make_array(m[2])
  z0 = make_array(m[2])
  zf = make_array(m[2])
  px0 = make_array(m[2])
  pxf = make_array(m[2])
  py0 = make_array(m[2])
  pyf = make_array(m[2])
  pz0 = make_array(m[2])
  pzf = make_array(m[2])
  
  for i=0,m[2]-1 do begin

    x0[i] = temporary(x[inis[i],i])
    y0[i] = temporary(y[inis[i],i])
    z0[i] = temporary(z[inis[i],i])
    px0[i] = temporary(px[inis[i],i])
    py0[i] = temporary(py[inis[i],i])
    pz0[i] = temporary(pz[inis[i],i])

    xf[i] = temporary((x[fines[i],i]))
    yf[i] = temporary((y[fines[i],i]))
    zf[i] = temporary((z[fines[i],i]))
    pxf[i] = temporary((px[fines[i],i]))
    pyf[i] = temporary((py[fines[i],i]))
    pzf[i] = temporary((pz[fines[i],i]))

  endfor

  remove, badind, x0, y0, z0, px0, py0, pz0, xf, yf, zf, pxf, pyf, pzf,fines,inis

  energy = py0^2*mi/(ee*2)+px0^2*mi/(ee*2)+pz0^2*mi/(ee*2)
;  energy=0

  energyf = pyf^2*mi/(ee*2)+pxf^2*mi/(ee*2)+pzf^2*mi/(ee*2)

  delta  =  pyf^2*mi/(ee*2)+pxf^2*mi/(ee*2)+pzf^2*mi/(ee*2) - (py0^2*mi/(ee*2)+px0^2*mi/(ee*2)+pz0^2*mi/(ee*2))

  ;points_in_bin=20
  bins = n_elements(y0)/points_in_bin
  disp=yf-y0
  disps=yf-(y0+shift*re)
  
  redisp=disp[sort(disp)]
  redisps=disps[sort(disps)]
  redelta=delta[sort(disp)]

  redisp=rebin(redisp[0:bins*points_in_bin-1],bins)
  redisps=rebin(redisps[0:bins*points_in_bin-1],bins)
  redelta=rebin(redelta[0:bins*points_in_bin-1],bins)
  
  zero=1
  final=3.1
  if flagleavey then zero = 0

  junk = scatterplot(redispe/re,redeltae/1000,xtitle='Displacement along y [$R_E$]', ytitle='$Energy gain [keV]$',sym_filled=1,sym_size=.6,xrange=[zero,final],font_size=24,dim=[1820,980],layout=[2,1,1],margin=[.2,.15,.02,.05],axis_Style=1,name=['Simulation ions (no turbulence)'],yrange=[0,130])
  junkwhite = scatterplot(redispe/re,redeltae/1000,xtitle='Displacement along y [$R_E$]', ytitle='$Energy gain [keV]$',sym_filled=1,sym_size=.6,xrange=[zero,final],font_size=24,dim=[1820,980],layout=[2,1,1],margin=[.2,.15,.02,.05],axis_Style=1,name=['Individual Ions'],sym_color='white',/current,yrange=[0,130])
  junky = scatterplot(redispe/re,redeltae/1000,xtitle='Displacement along y [$R_E$]', ytitle='$Energy gain [keV]$',sym_filled=1,sym_size=.35,xrange=[zero,final],font_size=24,dim=[1820,980],layout=[2,1,1],margin=[.2,.15,.02,.05],axis_Style=1,name=['Individual Ions'],/current,yrange=[0,130])

  junk1=  plot(redispe/re,Ey0*redispes/1000000,overplot=1,color='red',thick=4,name=['$E_{y_0}$ potential drop'])


  junk2 = scatterplot(redisp/re,redelta/1000,xtitle='Displacement along y [$R_E$]',sym_filled=1,sym_size=0.6,xrange=[zero,final],font_size=24,/current,layout=[2,1,2],margin=[.12,.15,.05,.05],axis_Style=1,sym='star',name=['Simulation Ions (turbulence)'],sym_color='blue',yrange=[0,130])
  junkwhite = scatterplot(redisp/re,redelta/1000,xtitle='Displacement along y [$R_E$]',sym_filled=1,sym_size=.7,xrange=[zero,final],font_size=24,layout=[2,1,2],margin=[.12,.15,.05,.05],axis_Style=1,name=['simulation Ions'],sym_color='white',/current,sym='star',yrange=[0,130])
  junky = scatterplot(redisp/re,redelta/1000,xtitle='Displacement along y [$R_E$]',sym_filled=1,sym_size=.35,xrange=[zero,final],font_size=24,layout=[2,1,2],margin=[.12,.15,.05,.05],axis_Style=1,name=['Simulation Ions'],/current,sym_color='blue',sym='star',yrange=[0,130])

  junk3=  plot(redisp/re,Ey0*redisps/1000000,overplot=1,color='red',thick=4,name=['$E_{y_0}$ potential drop'])

  leg = legend(target=[junk,junk3,junk2],position=[.41,.93],font_size=21,sample_width=.06)


print,mean(disp)/ymax
 print,mean(fines-inis)*0.1 ;finis is in goodx space so it doesn't count

  junk.save,'dispacementlongleavey.png'




end