pro forloopcuts, trajectories,maskzero=maskzero,import_index=import_index, maskleavey=maskleavey

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
;example forloopcuts, trajectories,maskzero=maskzero,import_index=import_index, maskleavey=maskleavey
Ey0 = 2.7D;

if keyword_set(import_index) then begin
  match,where(maskzero),import_index,suba,subb
  
 if keyword_set(maskleavey) then match,where(maskleavey),import_index,suba,leavey
endif

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

  n = size(trajectories)
  
  fines = dindgen(n[3])
  inis = dindgen(n[3])


  steps = 20
  imax = 10000L ;

  x  = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN )
  y  = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  z  = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  px = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  py = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  pz = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)

  Ex = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  Ey = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  Ez = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  Bx = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  By = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)
  Bz = make_array(n[1] + 1,n[3], value = !VALUES.F_NAN)

  
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

  Ex[goodx[ini:fin],number] = trajectories[goodx[ini:fin], 6, number]
  Ey[goody[ini:fin],number] = trajectories[goody[ini:fin], 7, number]
  Ez[goodz[ini:fin],number] = trajectories[goodz[ini:fin], 8, number]
  Bx[goodx[ini:fin],number] = trajectories[goodx[ini:fin], 9, number]
  By[goody[ini:fin],number] = trajectories[goody[ini:fin], 10, number]
  Bz[goodz[ini:fin],number] = trajectories[goodz[ini:fin], 11, number]
 
  
  fines[number] = goodx[fin]
  inis[number]  = goodx[ini]

endelse

    endfor
    
    
    x0 = diag_matrix(x[inis[*],*])
    y0 = diag_matrix(y[inis[*],*])
    z0 = diag_matrix(z[inis[*],*])
    px0 = diag_matrix(px[inis[*],*])
    py0 = diag_matrix(py[inis[*],*])
    pz0 = diag_matrix(pz[inis[*],*])
    
    xf = diag_matrix(x[fines[*],*])
    yf = diag_matrix(y[fines[*],*])
    zf = diag_matrix(z[fines[*],*])
    pxf = diag_matrix(px[fines[*],*])
    pyf = diag_matrix(py[fines[*],*])
    pzf = diag_matrix(pz[fines[*],*])
    
    remove, badind, x0, y0, z0, px0, py0, pz0, xf, yf, zf, pxf, pyf, pzf,fines,inis

    
    middle   = where(abs(z0) lt zexitmax/10.,b)
    top      = where(abs(z0) gt zexitmax/1.5 ,b)
    left     = where(y0  lt -yexitmax/2.    ,b)
    right    = where(y0  gt  yexitmax/2.    ,b)
    outto     = where(x0  gt  xexitmax/2.    ,b)
    into     = where(x0  lt  -xexitmax/2.    ,b)

    slow    = where(pyf^2*mi/(ee*2)+pxf^2*mi/(ee*2)+pzf^2*mi/(ee*2) lt mean(pyf^2*mi/(ee*2)+pxf^2*mi/(ee*2)+pzf^2*mi/(ee*2),/nan))
    fastest = where(pyf^2*mi/(ee*2)+pxf^2*mi/(ee*2)+pzf^2*mi/(ee*2) gt 2*mean(pyf^2*mi/(ee*2)+pxf^2*mi/(ee*2)+pzf^2*mi/(ee*2),/nan),nf)
    fast    = where(pyf^2*mi/(ee*2)+pxf^2*mi/(ee*2)+pzf^2*mi/(ee*2) gt mean(pyf^2*mi/(ee*2)+pxf^2*mi/(ee*2)+pzf^2*mi/(ee*2),/nan),nf)
    slowest = where(pyf^2*mi/(ee*2)+pxf^2*mi/(ee*2)+pzf^2*mi/(ee*2) lt .5*mean(pyf^2+pxf^2+pzf^2,/nan),nf)
 
   if keyword_set(maskzero) then cross   = subb

    
    energy = py0^2*mi/(ee*2)+px0^2*mi/(ee*2)+pz0^2*mi/(ee*2)
    delta  =  pyf^2*mi/(ee*2)+pxf^2*mi/(ee*2)+pzf^2*mi/(ee*2) - (py0^2*mi/(ee*2)+px0^2*mi/(ee*2)+pz0^2*mi/(ee*2))
    ;graphic1= scatterplot(y0,z0,symbol='star')
    indices=fast
    indices=indgen(n_elements(y0))
    
    print, 'middle',correlate(energy[middle],delta[middle])
    print, 'top',correlate(energy[top],delta[top])
    print, 'left',correlate(energy[left],delta[left])
    print, 'right',correlate(energy[right],delta[right])
    print, 'outto',correlate(energy[outto],delta[outto])
    print, 'into',correlate(energy[into],delta[into])
    print, 'slow',correlate(energy[slow],delta[slow])
    
    print, 'fastest',correlate(energy[fastest],delta[fastest])
    print, 'fast',correlate(energy[fast],delta[fast])
    print, 'slowest',correlate(energy[slowest],delta[slowest])
    
   if keyword_set(maskzero) then print, 'cross',correlate(energy[cross],delta[cross])
    print, 'all',correlate(energy,delta)


    
    junk = scatterplot(abs(y0-yf)/re,delta/1000.,xtitle='Displacement along y [$R_E$]',ytitle='$\Delta W [KeV]$',dimensions=[1100,800],SYM_FILLED=1,sym_size=0.3,sym_color='blue',xrange=[0,3])
    junk = plot(abs(y0-yf)/re,Ey0*abs(y0-yf)/1000000.,overplot=1,color='red',thick=4)
    junk = scatterplot(abs(y0[cross]-yf[cross])/re,delta[cross]/1000,dimensions=[1100,800],SYM_FILLED=1,sym_size=0.3,overplot=1)
    junk = plot(abs(y0[cross]-yf[cross]),Ey0*abs(y0[cross]-yf[cross])/1000000.,overplot=1,color='red',thick=4)
    junk = scatterplot(energy[leavey],delta[leavey],xtitle='$W_0$',ytitle='$\DeltaW$',dimensions=[1100,800])
    junk = scatterplot(energy,delta,xtitle='$Y_0$',ytitle='$\DeltaW [eV]$',dimensions=[1100,800],font_size=23)
    junk = scatterplot(energy[slow],delta[slow],xtitle='$W_0$',ytitle='$\DeltaW [eV]$',dimensions=[1100,800],font_size=23)
    junk = scatterplot(indices,abs(y0-yf),xtitle='time',ytitle='$\Deltay$',dimensions=[1100,800],font_size=23)



    points_in_bin = 5 ;the region has to be such that the variation in the rad can be taken to be constant
    bins     = n_elements(y0)/points_in_bin
   ; rey0 = y0[sort(y0)]
   ; reyf = yf[sort(y0)]
   ; redelta = delta[sort(y0)]
   disp = yf-y0
   print,mean(disp)
     redisp=disp[sort(disp)]
     redelta=delta[sort(disp)]
    rey0= rebin(y0[0:bins*points_in_bin-1], bins)
    reyf= rebin(yf[0:bins*points_in_bin-1], bins)
    redisp= rebin(redisp[0:bins*points_in_bin-1], bins)
    redelta= rebin(redelta[0:bins*points_in_bin-1], bins)
    
    junk = scatterplot(abs(rey0-reyf)/re,redelta/1000.,xtitle='Displacement along y [$R_E$]',ytitle='$\Delta W [KeV]$',dimensions=[1100,800],SYM_FILLED=1,sym_size=0.3,sym_color='blue',xrange=[0,3],font_size=21)
    junk = plot(abs(rey0-reyf)/re,Ey0*abs(rey0-reyf)/1000000.,overplot=1,color='red',thick=4)
    
    junk = scatterplot(redisp/re,redelta/1000.,xtitle='Displacement along y [$R_E$]',ytitle='$\Delta W [KeV]$',dimensions=[1100,800],SYM_FILLED=1,sym_size=0.3,sym_color='blue',xrange=[0,3],font_size=24)
    junk = plot(redisp/re,Ey0*redisp/1000000.,overplot=1,color='red',thick=4)




    ;plot,reradius
 ;   re = findradius(reradius, r1,r2)
 ;   reradius = reradius[re[0]:re[1]]


  
   graphic1= scatterplot(y0[indices],z0[indices],symbol='star',xrange=[yexitmin,yexitmax],yrange=[zexitmin,zexitmax],xtitle=['Y'],ytitle=['Z'])
   graphic1= scatterplot(yf[indices],zf[indices],symbol='star',xrange=[yexitmin,yexitmax],yrange=[zexitmin,zexitmax],overplot=1,sym_color='red')
   graphic1.save,'slowperticles.png'
   
   

   
    
    
    energy0 = mean(px[0,*]^2 + py[0:*]^2 + pz[0:*]^2,/nan)
    energy1 = mean(diag_matrix(px[fines,*])^2 + diag_matrix(py[fines,*])^2 + diag_matrix(pz[fines,*])^2)

end