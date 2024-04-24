pro plottrajectoriesw

  restore, 'trajectories.sav'

  flag=0
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

  Re    = 6374.d3

  ; SET UP BOX
  xmin     = -8.0D*Re ; LOCATION OF INJECTION BOUNDARY (PARTICLE ENTRY)
  xmax     =  8.0D*Re ;
  xexitmin = -8.0D*Re ; LOCATION PARTICLE EXIT BOUNDARY
  xexitmax =  8.0D*Re
  xturbmin = -7.75D*Re ; LOCATION OF TURBULENCE REGION
  xturbmax =  7.75D*Re

  zmin     = -2D*Re  ;
  zmax     =  2D*Re
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

  if flag then begin
    ixplus  = where(radiation[*, 0] GT xexitmax, nxp) ;if there's a provblem the problem is at the level of how radiatio is defined
    maskleavex =  maskleavey
    maskleavex[*]=0
    maskleavex[ixplus] = 1
    match,where(maskleavex),import_index,suba,leavex

    pxfast =where(radiation[*, 3]^2*mi/(ee*2) + radiation[*, 4]^2*mi/(ee*2)+ radiation[*, 5]^2*mi/(ee*2) gt 2*mean(radiation[*, 3]^2*mi/(ee*2) + radiation[*, 4]^2*mi/(ee*2)+ radiation[*, 5]^2*mi/(ee*2) ))
    maskfast =  maskleavey
    maskfast[*]=0
    maskfast[pxfast] = 1

    match,where(maskleavex*maskfast),import_index,suba2,fastx
    mask = where(maskleavex*maskfast)
    indices = mask[suba2]

    ;  trajectories = trajectories[*,*,fastx]
  endif

  n = size(trajectories)

  drift_towards  =  [15];[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,4564,3,53,5,2,41,4141,241,75,85,30];[5],[1142];[1065];[4503];[4392]; [4503,4350,4283];,4019,3952,3943,3336,3055,2450,1405,1208,1142,283,87] [2] -> is the crazy beamlet
  ;drift_across   = [4430, 4397, 4044];, 2842, 2144, 2139, 1541, 1659, 1631, 1563, 1535, 1447, 1116, 982]
  ;beamlets       = [65];[0];[0, 461, 609, 4392];, 4196, 4178, 4105, 4077, 3912, 3901, 3875, 3560, 2745, 1678, 998, 633]
  ;reverse_current= [4354,4308,4103,4097]
  ;other          = [182,4189,4103];,3674,2703,1655]
  jets = [939];[939];[235];,939]
  ;n3 = [drift_towards, drift_across, beamlets, reverse_current, other]
  n3 = [drift_towards, jets]

for i=0, n[3]-1 do begin
  
  junk = ddsega_plot_trajectory(trajectories, i, box_str,ddd=1)
;  junk = ddsega_plot_trajectory(trajectories, n3[1], box_str,ddd=1,overplot=1)


  cd, 'dumpw'
  write_png,'traje0b0.png',tvrd(/true)
  cd,'..'

 endfor

end