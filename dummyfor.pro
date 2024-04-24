pro dummyfor, import_index, trajectories,px,py,pz,fines,x,y,z,ex,ey,ez,bx,by,bz

  ee    = 1.6d-19
  me    = 0.911d-30
  twopi = 2D*!dpi
  mu0   = 4d-7*!dpi
  mi    = 1.67d-27

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
   ;  goody = where(trajectories[*, 1, j] NE 0)
   ;  yi[j] = trajectories[goody[2], 1, j]
    fines[j] = fin
  endfor
  
  energy0 = (px[0,*]^2 + py[0:*]^2 + pz[0:*]^2)*mi/ee
  energy1 = (diag_matrix(px[fines,*])^2 + diag_matrix(py[fines,*])^2 + diag_matrix(pz[fines,*])^2)*mi/ee
 ;   energy1 = px[fines,fines]^2 + py[fines,fines]^2 + pz[fines,fines]^2

  end