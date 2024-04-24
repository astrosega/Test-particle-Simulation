; ### PRO TO PLOT IMAGE IN BOX ###
pro imlplot, image, trans=trans, xs=xs, ys=ys, mag=mag, $
  xtitle=xtitle, xrange=xrange, ytitle =ytitle, yrange=yrange, $
  ztitle=ztitle, zrange=zrange, zlog=zlog, title=title, $
  charsize=charsize, vals=vals, colmax=colmax, colmin=colmin, $
  reverse_color=reverse_color, xmag=xmag, ymag=ymag, nwindow=nwindow

  ; NOTE XS AND YS ARE IN SCREEN PIXELS

  ; DEFAULT MARGINS
  xm = 50 + !P.charsize*25
  xp = 140 + !P.charsize*10; LEAVE ROOM OF COLORBAR
  ym = 50 + !p.charsize*10
  yp = 30 + !p.charsize*10
  colbargap   = 25
  colbarwidth = 15
  ztickwidth  = 5
  zticklable  = 50
  zlablegap   = 125

  ; CHECK KEYWORDS
  if keyword_set(trans) then iml = transpose(image) else iml=image
  if (keyword_set(xs) OR keyword_set(ys)) then regrid=1
  if (keyword_set(xs) AND keyword_set(ys)) then mag = 0
  if not keyword_set(colmax) then colmax = 254
  if not keyword_set(colmin) then colmin = 0
  if not keyword_set(nwindow) then nwindow = 0


  ; DETERMINE IMAGE SIZE
  temp    = size(iml)                        ; IMAGE SIZE
  if not keyword_set(xs) then xs = temp[1]   ; XSIZE
  if not keyword_set(ys) then ys = temp[2]   ; XSIZE

  ; CHECK FOR MAGNIFICATION
  if keyword_set(xmag) then xs = xs*xmag else $
    if keyword_set(mag) then xs = xs*mag
  if keyword_set(ymag) then ys = ys*ymag else $
    if keyword_set(mag) then ys = ys*mag
  if keyword_set(mag) OR keyword_set(xmag) OR keyword_set(ymag) then regrid=1

  ; REGRID
  if keyword_set(regrid) then iml = congrid(iml, xs, ys, /interp)

  ; CHECK ZRANGE
  if keyword_set(zlog) then iml = alog10(iml)
  if (keyword_set(zlog) AND keyword_set(zrange)) then zrange=alog10(zrange)
  if not keyword_set(zrange) then zrange = [min(iml), max(iml)]
  if zrange(0) EQ zrange(1) then zrange(1) = zrange(0) + 1

  ; NORMALIZE VALUES
  iml = (colmax-colmin)*(iml-zrange(0))/(zrange(1)-zrange(0)) + colmin
  if keyword_set(reverse_color) then iml = colmax-(iml-colmin)
  iml = iml > colmin
  iml = iml < colmax

  ; CHECK XRANGE and YRANGE
  if not keyword_set(xrange) then xrange = [0, xs]
  if not keyword_set(yrange) then yrange = [0, ys]

  ; MAKE WINDOW OF PROPER SIZE
  if (!d.name EQ 'X') OR (!d.name EQ 'WIN') then $
    window, Nwindow, xsize=xs+xm+xp, ysize=ys+ym+yp  ; FIX WINDOW SIZE
  if (!d.name EQ 'WIN') then erase, 0

  ; MAKE COLORBAR
  colbar = bindgen(1,colmax-colmin) + colmin
  if keyword_set(reverse_color) then colbar = colmax-(colbar-colmin)
  colbar = congrid(colbar, colbarwidth, ys)

  ; POSTSCRIPT FIX
  IF (!d.name EQ 'PS') THEN BEGIN
    ColorbarRat = (xm > (0.125*xs) ) /xm
    xm = xm*ColorbarRat
    xp = xp*ColorbarRat
    ym = ym*ColorbarRat
    yp = yp*ColorbarRat
    scale = (float(!d.x_size)/(xs+xm+xp)) < (float(!d.y_size)/(ys+ym+yp))
    xm = xm*scale
    xp = xp*scale
    ym = ym*scale
    yp = yp*scale
    colbargap   = colbargap*scale*ColorbarRat
    colbarwidth = colbarwidth*scale*ColorbarRat
    xs = xs*scale
    ys = ys*scale
    zlablegap   = zlablegap*scale*ColorbarRat
    ztickwidth  = ztickwidth*scale*ColorbarRat
    zticklable  = zticklable*scale*ColorbarRat
  ENDIF

  tickstart = colbargap+colbarwidth
  tickend   = tickstart+ztickwidth

  ; PLOT
  tv,iml,xm,ym, xsize=xs, ysize=ys
  tv, colbar, xm+xs+colbargap, ym, xsize=colbarwidth, ysize =ys
  
  backColor = cgColor("White")
  
  plot, [xrange(0),xrange(1)], [yrange(0),yrange(1)], /nodata, xstyle=1, ystyle=1, $
    position=[xm,ym,xs+xm,ys+ym], /device, /noerase, xrange=xrange, yrange=yrange, $
    xtitle=xtitle, ytitle =ytitle, title=title, charsize=charsize, background=backcolor

  ; MAKE COLORBAR BOX
  xrat = float(xrange(1) - xrange(0))/xs
  xr1  = xrange(1)
  x = [xr1+colbargap*xrat, xr1+colbargap*xrat, $
    xr1+(colbargap+colbarwidth)*xrat, xr1+(colbargap+colbarwidth)*xrat, $
    xr1+colbargap*xrat]
  y = [yrange(0), yrange(1), yrange(1), yrange(0), yrange(0)]
  oplot, x, y, /noclip

  ; TICKS AND ZLABLES
  ff_ylim, '', zrange, vals=vals, names=ytickname, yticks=yticks
  xx = [xr1+tickstart*xrat, xr1+tickend*xrat]
  vals = yrange(0) + (yrange(1)-yrange(0))*(vals-zrange(0))/(zrange(1)-zrange(0))
  for i =0, n_elements(vals)-1 DO oplot, xx, [vals(i), vals(i)], /noclip

  xx = xr1+zticklable*xrat
  for i =0, n_elements(vals)-1 DO $
    xyouts, xx, vals(i), strcompress(ytickname(i), /rem), charsize=charsize

  xx = xr1+zlablegap*xrat
  yy = total(yrange)/2.0
  if keyword_set(ztitle) then $
    xyouts, xx, yy, ztitle, align=0.5, orien = 90.0, charsize=charsize
end

; ### END