; MMS_R06_FITER
; FILTERS DATA TO A FREQUENCY RANGE
;
; FREQ MUST BE A VECTOR

pro mms_r06_filter_l2, t, dat, freq, iterate=iterate, tr=tr, tf=tf, datf=datf

; ESTABLISH NYQUIST FREQUENCY
npts = n_elements(t)
dt = median( t(1:(100<(npts-1))) - t(0:(99<(npts-2))) )
fnq = 1.0/(2.0*dt)

; ESTABLISH FHIGH AND FLOW AS FUNCTION OF NYQUIST
flow = (freq(0)/fnq) < 0.9999
fhigh = (freq(1)/fnq) < 1.0

; DO CASE IF NO REDUCE
IF not keyword_set(iterate) then BEGIN
  ; BRUTE FORCE
  if (flow EQ 0.0) then nterms = round(!pi/fhigh) else $
                        nterms = round(!pi/flow) < (npts-2)/2
  cof  = digital_filter(Flow, Fhigh, 50.0, Nterms)
  if (flow EQ 0.0) then cof  = cof/total(cof)
  datf = convol(dat, cof, /edge_t, /nan) 
  tf   = t
  if keyword_set(tr) then datf = interpol(datf, tf, tr)
  if keyword_set(tr) then tf = tr 
ENDIF ELSE BEGIN
  
  ; SET UP ITERATE
  kerl = [-8.0, 0.0, 72.0, 128.0, 72.0, 0.0, -8.0]/256.0
  kers = [0.25, 0.5, 0.25]
  datr = dat
  tr   = t
  nptr = npts

  ihigh = (round(alog(fnq)/alog(2.0)) - ceil(alog(freq(1))/alog(2.0))) < $
          (floor(alog(npts)/alog(2.0))-3) ; CORRECTED REE 17-10-17
  FOR i=0,ihigh-1 do BEGIN & $
    datr = convol(datr, kerl, /edge_t, /nan) & $
    datr = convol(datr, kers, /edge_t, /nan) & $
    ind = lindgen((nptr-1)/2)*2 + 1 & $
    datr=datr(ind) & $
    tr = tr(ind) & $
    nptr = n_elements(ind) & $
  ENDFOR

  IF (freq(0) GT 0) then BEGIN
    ilow = (round(alog(fnq)/alog(2.0)) - floor(alog(freq(0))/alog(2.0))) < $
          (floor(alog(npts)/alog(2.0))-3) ; CORRECTED REE 17-10-17
    datl = dat
    tl   = t
    nptl = npts

    FOR i=0,ilow-1 do BEGIN & $
      datl = convol(datl, kerl, /edge_t, /nan) & $
      datl = convol(datl, kers, /edge_t, /nan) & $
      ind = lindgen((nptl-1)/2)*2 + 1 & $
      datl=datl(ind) & $
      tl = tl(ind) & $
      nptl = n_elements(ind) & $
    ENDFOR

    datl = interpol(datl, tl, tr, /spline)
    datr = datr - datl
  ENDIF
 
  datf = datr
  tf = tr
ENDELSE

return
end







