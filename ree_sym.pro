pro ree_sym, z, square=square, dash=dash, cross=cross, fill=fill

  if not keyword_set(z) then z=1
  z = z(0)

  IF keyword_set(dash) then BEGIN
    x = [-z, z]
    y = [0.0, 0.0]
    usersym, x, y
    RETURN
  ENDIF

  IF keyword_set(cross) then BEGIN
    x = [-z,    z, 0.0, 0.0, 0.0]
    y = [0.0, 0.0, 0.0,  -z,   z]
    usersym, x, y
    RETURN
  ENDIF

  IF keyword_set(square) then BEGIN
    x=[-z,z,z,-z,-z]
    y=[z, z,-z,-z,z]
    usersym, x, y, fill=fill
  ENDIF ELSE BEGIN
    x1 = findgen(21)*!dpi/10
    y1 = z*sin(x1)
    x1 = z*cos(x1)
    usersym, x1, y1, fill=fill
  ENDELSE

end