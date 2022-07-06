; SETUP (WITH SPEDAS)
;setenv, 'ROOT_DATA_DIR=/Users/ree/mmsdata/'
;mms_init, local_data_dir = '/Users/ree/mmsdata'

; FIX COLORS
;tvlct, red, green, blue, /get
;red(6) = 225
;green(4) = 135
;red(5) = 255
;green(5) = 200
;tvlct, red, green, blue
device, decomp = 0
!p.background = 255
!p.color = 0
!p.charsize = 1.75

; SETUP - PLEASE LOAD SPEDAS IN YOU PATH. IF NOT
loadct, 39; red = 254, blue = 50, green = 130

; CONSTANTS
ee = 1.6D-19
me = 0.911D-30
twopi = 2D*!dpi
n = 0.05D6
mu0 = 4d-7*!dpi
mi = 1.67D-27
RE = 6374D3
ds = 1e-2

; DEFINE BOX
zmax = 2.5*Re
zmin = -2.5*Re
xmax = 10*Re
xmin = -10*Re

; SET UP LOBE AND Z FIELD STRENGTHS
Bx0   = 20D-9           ; Bx STRENGTH IN LOBES
z0    = 0.25D * RE    ; THICKNESS OF CURRENT SHEET
Bz0   = 2.5d-9          ; Bz STRENGTH AT BOUNDARIES
X0    = 2.5D * RE   ; USE 10:1 RATIO FOR MAGNETIC RECONNECTION

; LINE PLOTS OF B
; MAKE AXIS
nax   = 1001L
nhalf = (nax-1)/2

zax   = (dindgen(nax)-nhalf)*(zmax-zmin)/(nax-1)
Bx  = Bx0*tanh(zax/Z0)  ; HARRIS LIKE
;plot, zax/re, Bx/Bx0

xax = (dindgen(1001)-500D)*(xmax-xmin)/1000D
Bz  = Bz0*tanh(xax/X0)
;plot, xax/re, Bz/Bz0



; BELOW IS MAKES AN IMAGE OF B - NOT GENERALIZED. USE CONTOUR???
; MAKE A 2D ARRAY OF MAGNITUDE OF B

nx = 401 ; 20 RE
nz = 101 ; 5 RE
cal = 1D/20D

Bmag = fltarr(nx, nz) ; RATIO IS IN KEEPING WITH XMAX AND ZMAX
for j = 0, nx-1 do $
  for k=0, nz-1 do $
  Bmag(j,k) = sqrt( (Bx0*tanh( (float(k)*cal*RE + Zmin)/z0 ) )^2.0 + $
  (Bz0*tanh( (float(j)*cal*RE + Xmin)/x0 ) )^2.0 )
;SET UP PLOT
loadct, 1
!p.charsize=1.4

; MAKE SURE THAT RRE_SYM and IMLPLOT ARE IN PATH
ree_sym, 0.01
imlplot, Bmag*1d9, mag=2.5, colmin=160, colmax=240, xran=[Xmin/RE, Xmax/RE], $
  yran = [Zmin/Re,Zmax/RE], xtit = 'X (R!DE!N)', zran = [0.0, 20.0], $
  ytit = 'Z (R!DE!N)', ztit = '|B|', /reverse_color

; NEXT SECTION ADDS LINES - BRUTE FOR TRACING
FOR dum=0, 11 do BEGIN & $
  if (dum EQ 0) then sz0 = [0.6D]
  BxB = Bx0*tanh(sz0*Re/z0)
  sz0 = sz0 + Bx0/abs(BxB)*0.25 & $
  sz = sz0 & $
  sgn = (sz/(abs(sz)>1d-19) < 1D) > (-1D) & $
  sx = [-10D] & $
  FOR k=0, 2500 do BEGIN & $
  Bx = Bx0*tanh(sz*Re/Z0) & $
  Bz = Bz0*tanh(sx*Re/X0) & $
  B = sqrt(Bx*Bx + Bz* Bz) & $
  sx = sx + sgn*Bx/B*ds & $
  sz = sz + sgn*Bz/B*ds & $
  if (sx GE -10D) AND (sx LE 10D) AND (sz GE -2.5D) AND (sz LE 2.5D) then $
  oplot, sx, sz, psym=3 & $
  if (sx GE -10D) AND (sx LE 10D) AND (sz GE -2.5D) AND (sz LE 2.5D) then $
  oplot, -sx, sz, psym=3 & $
ENDFOR & $
ENDFOR

FOR dum=0, 11 do BEGIN & $
  BxB = Bx0*tanh(sz0*Re/z0) & $
  if (dum EQ 0) then sz0 = [-0.6D] else $
  sz0 = sz0 - Bx0/abs(BxB)*0.25 & $
  sz = sz0 & $
  sgn = (sz/(abs(sz)>1d-19) < 1D) > (-1D) & $
  sx = [-10D] & $
  FOR k=0, 2500 do BEGIN & $
  Bx = Bx0*tanh(sz*Re/Z0) & $
  Bz = Bz0*tanh(sx*Re/X0) & $
  B = sqrt(Bx*Bx + Bz* Bz) & $
  sx = sx + sgn*Bx/B*ds & $
  sz = sz + sgn*Bz/B*ds & $
  if (sx GE -10D) AND (sx LE 10D) AND (sz GE -2.5D) AND (sz LE 2.5D) then $
  oplot, sx, sz, psym=3 & $
  if (sx GE -10D) AND (sx LE 10D) AND (sz GE -2.5D) AND (sz LE 2.5D) then $
  oplot, -sx, sz, psym=3 & $
ENDFOR & $
ENDFOR

END