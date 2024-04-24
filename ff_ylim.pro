;+
; PROCEEDURE: FF_YLIM, pname, ylim, log=log, vals=vals, clear=clear, dsp=dsp
;
;
; PURPOSE: A utility to set up y-ticks and lables for fields plots.
;          ONLY WORKS FOR INTEGER RANGES! IF YOU DON'T LIKE WHAT YOU GET,
;          DO IT YOUR SELF (6 commands).
; INPUT:
;       pname -       REQUIRED. A tplot name.
;       ylim -        REQUIRED. A vetor yrange.
;
; KEYWORDS:
;       log -         If set, yrange is log otherwise linear.
;       vals -        If given, they are the values!
;       dsp -         If set, a special set of lables comes up.
;
; CALLING: ff_ylim,'SFA_OMNI',[0,600]
;
; OUTPUT: Stored in tplot.
;
; INITIAL VERSION: REE 97-11-18
; MODIFICATION HISTORY:
; Space Scienes Lab, UCBerkeley
;
;-
;

pro ff_ylim, pname, ylim, log=log, vals=vals, names=names, yticks=yticks, six=six

  IF n_elements(ylim) NE 2 then BEGIN
    if keyword_set(vals) then ylim=[min(vals),max(vals)] $
    ELSE BEGIN
      print, "FF_YLIM: STOPPED! Must give vector yrange."
      return
    ENDELSE
  ENDIF

  ; DO LOG CASE
  IF keyword_set(log) then BEGIN

    ; HANDLE CASE WHERE VALS ARE GIVEN
    IF keyword_set(vals) then BEGIN
      FOR i=0, n_elements(vals)-1 DO BEGIN
        x = round(alog10(vals(i)))
        IF not keyword_set(names) then BEGIN
          name = strcompress('10!A' + string(x) + '!N', /rem)
          if i EQ 0 then names = ' ' + name else names=[names,' ' + name]
        ENDIF
      ENDFOR
    ENDIF ELSE BEGIN

      ; VALS NOT SET
      ; DETERMINE HOW MANY AND WHERE TICKS GO.
      max = floor(alog10(ylim(1)) + 0.000001)
      min = ceil (alog10(ylim(0)) - 0.000001)

      nticks = max - min + 1
      IF nticks LT 1 then BEGIN
        print, "FF_YLIM: STOPPED! Invalid range."
        return
      ENDIF

      FOR i=min,max DO BEGIN
        name = strcompress('10!A' + string(i) + '!N', /rem)
        if i EQ min then names = ' ' + name else names=[names,' ' + name]
        if i EQ min then vals  = 10.0^i else vals=[vals,10.0^i]
      ENDFOR
    ENDELSE

    ;SET OPTIONS
    IF pname NE '' then BEGIN
      options,pname,'yticks',n_elements(vals)-1
      options,pname,'ytickname',names
      options,pname,'ytickv',vals
      options,pname,'ystyle',1
      options,pname,'ylog',1
      options,pname,'yrange',ylim
    ENDIF

    yticks = n_elements(vals) -1

  ENDIF ELSE BEGIN ; LINEAR CASE

    ; HANDLE CASE WHERE VALS ARE GIVEN
    IF keyword_set(vals) then BEGIN
      FOR i=0, n_elements(vals)-1 DO BEGIN
        x = round(vals(i))
        IF not keyword_set(names) then BEGIN
          name = strcompress(string(x), /rem)
          if i EQ 0 then names = ' ' + name else names=[names,' ' + name]
        ENDIF
      ENDFOR
    ENDIF ELSE BEGIN

      ; DETERMINE HOW MANY AND WHERE TICKS GO.
      range = ylim(1) - ylim(0)
      IF range LE 0 then BEGIN
        print, "FF_YLIM: STOPPED! Invalid range."
        return
      ENDIF

      mag = floor(alog10(range))
      mag = 10.0^mag
      ran = range/mag

      delt = 0.5
      if ran GT 1.2 then delt = 0.5
      if ran GT 3.0 then delt = 1.0
      if ran GT 5.5 then delt = 2.0

      ; SPECIAL CASE (SIX)
      if (keyword_set(six) AND ran EQ 1.2) THEN delt = 0.2

      ; SPECIAL CASE [-5,5], [-6,6], [0,50]
      if ran EQ 1.2 AND ylim(0)*10/mag EQ -6 AND ylim(1) GT 1 then delt = 0.30
      if ran EQ 1.0 AND ylim(0)*10/mag EQ -5 AND ylim(1) GT 10 then delt = 0.25
      if ran EQ 5.0 AND ylim(0) EQ 0 AND ylim(1) GT 10 then delt = 2.5

      min = ceil ( ylim(0)/(delt*mag) - 0.000001 )
      max = floor( ylim(1)/(delt*mag) + 0.000001 )

      FOR i=min,max DO BEGIN
        name = strcompress(string(floor(i*delt*mag + 0.01)), /rem)
        IF (mag LE 1.1) THEN BEGIN
          name = string(i*delt*mag,format='(F5.1)')
        ENDIF
        IF (mag LE 0.11)THEN BEGIN
          name = string(i*delt*mag,format='(F5.2)')
        ENDIF
        IF (mag LE 0.011)THEN BEGIN
          name = string(i*delt*mag,format='(F6.3)')
        ENDIF
        if i EQ min then names = ' ' + name else names=[names,' ' + name]
        if i EQ min then vals  = i*delt*mag else vals=[vals,i*delt*mag]
      ENDFOR
    ENDELSE

    ylim(0) = ylim(0) < vals(0)
    ylim(1) = ylim(1) > max(vals)

    ;SET OPTIONS
    IF pname NE '' then BEGIN
      options,pname,'yticks',n_elements(vals)-1
      options,pname,'ytickname',names
      options,pname,'ytickv',vals
      options,pname,'ystyle',1
      options,pname,'ylog',0
      options,pname,'yrange',ylim
    ENDIF

    yticks = n_elements(vals) -1

  ENDELSE ; Linear case.

  return
end