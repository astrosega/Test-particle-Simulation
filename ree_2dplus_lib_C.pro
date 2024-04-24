; PROS FOR 2DPlus ION SIMULATION


; %%%%%%%%%%%%%%%%%%%%%%%%% IONS 3D %%%%%%%%%%%%%%%%%%%%%%%%%%

; %%%%%%%%%%%%%%%%%%%%%%%%% IONS VIA SPECTRA %%%%%%%%%%%%%%%%%%%%%%%%%%

; CREATE POINT TIME SERIES

pro ion_e_set_up, fmin, fmax, nfreqs=nfreqs, Erms=Erms, $
  alpha=alpha, ntimes=ntimes, nit=nit, EIonT=EIonT

twopi   = 2D*!dpi
; CHECK KEYWORDS
if not keyword_set(fmin) then fmin = 0.01D
if not keyword_set(fmax) then fmax = 50D
if not keyword_set(nfreqs) then nfreqs=5000L
if not keyword_set(Erms) then Erms=10D
if not keyword_set(alpha) then alpha=1.1D
if not keyword_set(ntimes) then ntimes=1200000L
if not keyword_set(nit) then nit=10L
;print, fmin, fmax, nfreqs, Erms, alpha, ntimes, nit


; GENERATE FREQUENCIES, AMPLITUDES AND PHASES
df = (fmax-fmin)/nfreqs
ws     = twopi*((dindgen(nfreqs)+0.5D)*df + randomn(seed, nfreqs, /double)/8D*df) 
dum    = ((dindgen(nfreqs)+0.5D)*df + randomn(seed, nfreqs, /double)/32D*df)
phis   = randomu(seed, nfreqs, /double)*twopi
As     = 1D/(dum^(alpha/2D))
ind = where(dum LT 0.25D)
As(ind) = As(ind)*(1.325D - 1.3D*dum(ind))
ind = where(dum GT 0.25D)
As(ind) = As(ind)*(1.0D + 0.165D*alog(dum(ind)/0.25D))
Scale  = 2D*Erms/sqrt(total(As*As))


; GENERATE EionT
dum = 0D
t = dindgen(ntimes)*0.01
for jj = 0, nfreqs-1 do dum= dum + scale*As(jj)*cos(ws(jj)*t+phis(jj))
kern = [0.18D, 0.64D, 0.18D]
EIonT = convol(dum, kern)

RETURN
end


; %%%%%%%%%%%%%%%%%%%%%%%%% IONS VIA SPECTRA %%%%%%%%%%%%%%%%%%%%%%%%%%

; CREATE POINT TIME SERIES

pro ion_db_set_up, Brms=Brms, alpha=alpha, ntimes=ntimes, dbIonT=dbIonT

; CHECK KEYWORDS
if not keyword_set(Brms) then Brms=2D
if not keyword_set(alpha) then alpha = sqrt(2D^(-2.65D))
if not keyword_set(ntimes) then ntimes=1000000L

db0 = Brms*20D ; IMPERICALLY DERIVED

; GENERATE FLAT SIGNAL
dbbase = randomn(seed,ntimes,/double)*db0
t = dindgen(ntimes)*0.01D

; FILTER IN SETS
mms_r06_filter_l2, t, dbbase, [0.0, 0.2499], /iterate, tf=tf, datf=datf
db = interpol(datf, tf, t, /spline)
filt = [0.25001D, 0.49999D]
FOR i = 0, 7 do BEGIN & $
  mms_r06_filter_l2, t, dbbase, filt*(2D^i), /iterate, tf=tf, datf=datf & $
  db = db + interpol(datf, tf, t, /spline) * (alpha^(i+1)) & $
ENDFOR

kern = [0.18D, 0.64D, 0.18D]
dbIonT = convol(db, kern)

RETURN
end


; %%%%%%%%%%%%%%%%%%%%%%%%%  INJECT PARTICLES 3D SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%
pro inject_particles3D_setup, ninjarr, vth, DenZ0, DenInjZ, nEtimes, nBtimes, $
   Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Inj_str=Inj_str

; SET UP FOR INJECTION 
; NOTE THAT ZMIN MUST EQUAL -ZMAX FOR INJECTION TO WORK
DenExtra = ( DenInjZ - 1D/(cosh(Zmax/DenZ0)^2) ) > 0D
ZdistEff = 2D*tanh(Zmax/DenZ0)*DenZ0 ; INTEGRATION OF SECH^2

; CALCULATE AREAS; NOTE THAT X AND Y AREA ACCOUNTS FOR OUT-OF-BOUNDS PARTICLES
AreaX = (Ymax-Ymin)*(ZdistEff*(1D - DenExtra) +  DenExtra*(Zmax-Zmin))
AreaY = (Xmax-Xmin)*(ZdistEff*(1D - DenExtra) +  DenExtra*(Zmax-Zmin))
AreaZ = (Xmax-Xmin)*(Ymax-Ymin)*DenInjZ ; NORMALIZES Z INJECTION
AreaT = AreaX + AreaY + AreaZ

AxNorm = AreaX/AreaT
AyNorm = AreaY/AreaT
AzNorm = AreaZ/AreaT

; FILL ARRAYS WITH DEFAULTS
xinjarr = randomu(seed, ninjarr, /double)*(Xmax-Xmin) + Xmin
yinjarr = randomu(seed, ninjarr, /double)*(Ymax-Ymin) + Ymin

; Z IS COMPLICATED - MAKE SECH^2 ARRAY
dum     = randomu(seed, ninjarr,/double)
ztemp   = 0.5D * alog((1D + dum)/(1D - dum))
sgn     = 2L*round(randomu(seed, ninjarr,/double)) - 1L
zinjarr = ztemp*sgn*DenZ0 ; SECH^2 HISTROGRAM

; RANDOMLY INTERMIX FLAT ARRAY
thresh = DenExtra*(Zmax-Zmin)/(ZdistEff*(1D - DenExtra) +  DenExtra*(Zmax-Zmin))
indmix = where(randomu(seed, ninjarr, /double) LT thresh, nindmix)
if nindmix GT 0 then $
  zinjarr(indmix) = randomu(seed, nindmix, /double)*(Zmax-Zmin) + Zmin

; CREATE THERMAL ARRAYS
dum1     = randomn(seed, ninjarr, /double)*vth
dum2     = randomn(seed, ninjarr, /double)*vth
vth2arr  = sqrt(dum1*dum1+dum2*dum2)

; SET DEFAULT MOMENTUM
pxtharr  = randomn(seed, ninjarr, /double)*vth
pytharr  = randomn(seed, ninjarr, /double)*vth
pztharr  = randomn(seed, ninjarr, /double)*vth

; SELECT THE FACE
dum     =  randomu(seed, ninjarr, /double)

; NEGATIVE X FACE
ind = where(dum LT (AxNorm/2D))
xinjarr(ind) = Xmin+1d-3
pxtharr(ind) = vth2arr(ind)

; POSITIVE X FACE
ind = where( (dum GE (AxNorm/2D)) AND (dum LT AxNorm) )
xinjarr(ind) = Xmax-1d-3
pxtharr(ind) = -vth2arr(ind)

; NEGATIVE Y FACE
ind = where( (dum GE AxNorm) AND (dum LT (AxNorm+Aynorm/2D)) )
yinjarr(ind) = Ymin+1d-3
pytharr(ind) = vth2arr(ind)

; POSITIVE Y FACE
ind = where( (dum GE (AxNorm+Aynorm/2D)) AND (dum LT (AxNorm+Aynorm)) )
yinjarr(ind) = Ymax-1d-3
pytharr(ind) = -vth2arr(ind)

; NEGATIVE Z FACE
ind = where( (dum GE (AxNorm+Aynorm)) AND (dum LT (AxNorm+Aynorm+Aznorm/2D)) )
zinjarr(ind) = Zmin+1d-3
pztharr(ind) = vth2arr(ind)

; POSITIVE Z FACE
ind = where( (dum GE (AxNorm+Aynorm+Aznorm/2D)) )
zinjarr(ind) = Zmax-1d-3
pztharr(ind) = -vth2arr(ind)

; GAMMA
c = 3d8
c2 = c*c
gaminj = sqrt(pxtharr*pxtharr+pytharr*pytharr + pztharr*pztharr+c2)/c

; DS - POSITION IN ELECTRIC FIELD TURBULENCE
injds = long(randomu(seed, ninjarr)*(nEtimes-1))
inkds = long(randomu(seed, ninjarr)*(nBtimes-1))

Inj_str = {x:xinjarr, y:yinjarr, z:zinjarr, px:pxtharr, py:pytharr, $
           pz:pztharr, gam:gaminj, ninjarr:ninjarr, dsE:injds, dsB:inkds}

RETURN
end


; %%%%%%%%%%%%%%%%%%%%%%%%% LOAD PARTICLES 3D %%%%%%%%%%%%%%%%%%%%%%%%%%
pro load_particles3D, nload, vth, DenZ0, DenInjZ, nEtimes, nBtimes, $
  Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, $
  x=x, y=y, z=z, px=px, py=py, pz=pz, gam=gam, dsE=dsE, dsB=dsB, nstart=nstart, ind=ind

; CHECK KEYWORDS
if not keyword_set(nstart) then nstart = 0L
nend = nstart+nload-1

; FILL ARRAYS WITH DEFAULTS
x(nstart:nend) = randomu(seed, nload, /double)*(Xmax-Xmin) + Xmin
y(nstart:nend) = randomu(seed, nload, /double)*(Ymax-Ymin) + Ymin

; Z IS COMPLICATED - MAKE SECH^2 ARRAY
dum     = randomu(seed, nload,/double)
ztemp   = 0.5D * alog((1D + dum)/(1D - dum))
sgn     = 2L*round(randomu(seed, nload,/double)) - 1L
z(nstart:nend) = ztemp*sgn*DenZ0 ; SECH^2 HISTROGRAM

; RANDOMLY INTERMIX FLAT ARRAY
DenExtra = ( DenInjZ - 1D/(cosh(Zmax/DenZ0)^2) ) > 0D
ZdistEff = 2D*tanh(Zmax/DenZ0)*DenZ0 ; INTEGRATION OF SECH^2
thresh = DenExtra*(Zmax-Zmin)/(ZdistEff*(1D - DenExtra) +  DenExtra*(Zmax-Zmin))
indmix = where(randomu(seed, nload, /double) LT thresh, nindmix)
if nindmix GT 0 then $
  z(nstart+indmix) = randomu(seed, nindmix, /double)*(Zmax-Zmin) + Zmin

px(nstart:nend) = randomn(seed, nload, /double)*vth
py(nstart:nend)  = randomn(seed, nload, /double)*vth
pz(nstart:nend)  = randomn(seed, nload, /double)*vth

c = 3d8
c2 = c*c
gam(nstart:nend) = sqrt(px(nstart:nend)*px(nstart:nend) + $
  py(nstart:nend)*py(nstart:nend) + pz(nstart:nend)*pz(nstart:nend) +c2)/c
dsE(nstart:nend) = long(randomu(seed, nload)*(nEtimes-1))
dsB(nstart:nend) = long(randomu(seed, nload)*(nBtimes-1))

ind = where(gam NE 0D, nind)


RETURN
end


; %%%%%%%%%%%%%%%%%%%%%%%%% INJECT PARTICLES 3D %%%%%%%%%%%%%%%%%%%%%%%%%%
pro inject_particles3D, nnew, jt=jt, x=x, y=y, z=z, px=px, py=py, pz=pz, $
  gam=gam, dsE=dsE, dsB=dsB, Inj_str=Inj_str

  ; GENERATE INDOPEN
indopen = where(gam EQ 0)
indnew = indopen(0:nnew-1)
   
Ijmax = jt+nnew-1

; ROLL OVER IF jt EXCEEDS ARRAY SIZE
IF (Ijmax GE Inj_str.ninjarr) then BEGIN
  jt = 0L
  Ijmax = jt+nnew-1
ENDIF

x(indnew)   = Inj_str.x(jt:Ijmax)
y(indnew)   = Inj_str.y(jt:Ijmax)
z(indnew)   = Inj_str.z(jt:Ijmax)
px(indnew)  = Inj_str.px(jt:Ijmax)
py(indnew)  = Inj_str.py(jt:Ijmax)
pz(indnew)  = Inj_str.pz(jt:Ijmax)
gam(indnew) = Inj_str.gam(jt:Ijmax)
dsE(indnew)  = Inj_str.dsE(jt:Ijmax)
dsB(indnew)  = Inj_str.dsB(jt:Ijmax)

jt = jt + nnew

RETURN
end


; %%%%%%%%%%%%%%%%%%%% MAGNETIC FIELD %%%%%%%%%%%%%%%%%%%%%%%%%

pro bcalc, x, z, ind, B_str, Bx=Bx, Bz=Bz, Bmag=Bmag, $
  sechx=sechx, sechz=sechz, tanhx=tanhx, atan_sinhz=atan_sinhz

IF keyword_set(B_str.BH0) then BEGIN 
  sechx(ind) = 1D/cosh(x(ind)/B_str.xHw)
  sechz(ind) = 1D/cosh((z(ind)-B_str.zH0)/B_str.zHw)
  tanhx(ind) = tanh(x(ind)/B_str.xHw)
  atan_sinhz(ind) = atan(sinh((z(ind)-B_str.zH0)/B_str.zHw))
  Bx(ind) = (-B_str.BH0)*sechx(ind)*sechz(ind)
  Bz(ind) = (-B_str.BH0)*(B_str.zHw/B_str.xHw)* $
            tanhx(ind)*sechx(ind)*atan_sinhz(ind)
  Bx(ind) = Bx(ind) + B_str.Bx0*tanh(z(ind)/B_str.Z0)
  Bz(ind) = Bz(ind) + B_str.Bz0*tanh(x(ind)/B_str.X0)
ENDIF ELSE BEGIN
  Bx(ind) = B_str.Bx0*tanh(z(ind)/B_str.Z0)
  Bz(ind) = B_str.Bz0*tanh(x(ind)/B_str.X0)
ENDELSE

Bmag(ind) = sqrt(Bx(ind)*Bx(ind) + Bz(ind)*Bz(ind)) > 1d-12

RETURN


end

; %%%%%%%%%%%%%%%%%%%%%%%%% BORIS 1/2 STEP ADVANCE %%%%%%%%%%%%%%%%%%%%%%%%%%

pro boris_half_step, px, py, pz, gam, Bx, By, Bz, Ex, Ey, Ez, ind, bqmdt, qmdt, $
  pxn, pyn, pzn, pxold, pyold, pzold

; DO NOMINAL ADVANCE
pxn(ind) = px(ind) + Ex(ind)*qmdt + $
          (py(ind)*Bz(ind)-pz(ind)*By(ind))*bqmdt/gam(ind) & $
pyn(ind) = py(ind) + Ey(ind)*qmdt + $
          (pz(ind)*Bx(ind)-px(ind)*Bz(ind))*bqmdt/gam(ind) & $
pzn(ind) = pz(ind) + Ez(ind)*qmdt + $
          (px(ind)*By(ind)-py(ind)*Bx(ind))*bqmdt/gam(ind) & $

; RECORD OLD VALUES
pxold(ind) = px(ind) & $
pyold(ind) = py(ind) & $
pzold(ind) = pz(ind) & $

; HALF-STEP ADVANCE
  px(ind) =  px(ind) + Ex(ind)*qmdt + ( (pyold(ind)+pyn(ind))*Bz(ind)- $
                       (pzold(ind)+pzn(ind))*By(ind))*bqmdt/gam(ind)/2D & $
  py(ind) =  py(ind) + Ey(ind)*qmdt + ( (pzold(ind)+pzn(ind))*Bx(ind)- $
                       (pxold(ind)+pxn(ind))*Bz(ind))*bqmdt/gam(ind)/2D & $
  pz(ind) =  pz(ind) + Ez(ind)*qmdt + ( (pxold(ind)+pxn(ind))*By(ind)- $
                       (pyold(ind)+pyn(ind))*Bx(ind))*bqmdt/gam(ind)/2D & $

RETURN
end

; %%%%%%%%%%%%%%%%%%%%%%%%% BORIS ADVANCE %%%%%%%%%%%%%%%%%%%%%%%%%%

pro boris_advance, px, py, pz, gam, Bx, By, Bz, Ex, Ey, Ez, ind, bqmdt, qmdt, $
  Ax, Ay, Az, Bgx, Bgy, Bgz, det
  
; CALCULATE VELOCITIES
Bgx(ind) = Bx(ind)*bqmdt/gam(ind)/2D
Bgy(ind) = By(ind)*bqmdt/gam(ind)/2D
Bgz(ind) = Bz(ind)*bqmdt/gam(ind)/2D

; ADVANCE FIRST 1/2 STEP
Ax(ind) = px(ind) + Ex(ind)*qmdt + (py(ind)*Bgz(ind)-pz(ind)*Bgy(ind))
Ay(ind) = py(ind) + Ey(ind)*qmdt + (pz(ind)*Bgx(ind)-px(ind)*Bgz(ind))
Az(ind) = pz(ind) + Ez(ind)*qmdt + (px(ind)*Bgy(ind)-py(ind)*Bgx(ind))

; CALCULATE DETERMINATE
det(ind) = 1D + Bgx(ind)*Bgx(ind) + Bgy(ind)*Bgy(ind) + Bgz(ind)*Bgz(ind)

; ADVANCE
px(ind) = ( Ax(ind)*(1D + Bgx(ind)*Bgx(ind)) + $
            Ay(ind)*(Bgy(ind)*Bgx(ind) + Bgz(ind)) + $
            Az(ind)*(Bgz(ind)*Bgx(ind) - Bgy(ind)) ) / det(ind)
py(ind) = ( Ax(ind)*(Bgx(ind)*Bgy(ind) - Bgz(ind)) + $
            Ay(ind)*(1D + Bgy(ind)*Bgy(ind)) + $
            Az(ind)*(Bgz(ind)*Bgy(ind) + Bgx(ind)) ) / det(ind)
pz(ind) = ( Ax(ind)*(Bgz(ind)*Bgx(ind) + Bgy(ind)) + $
            Ay(ind)*(Bgz(ind)*Bgy(ind) - Bgx(ind)) + $
            Az(ind)*(1D + Bgz(ind)*Bgz(ind)) ) / det(ind)
RETURN
end


