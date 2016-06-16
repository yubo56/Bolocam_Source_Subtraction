function subtractmax_multi, sig, specDens, sigm, binwidth, emissivity, tdust, normband=normband, bbody=bbody, eq_sigm=eq_sigm, real_x=real_x, real_y=real_y, real_amp=real_amp
compile_opt idl2, HIDDEN
    ; Input
    ;   sig         - struct containing
    ;       signal  - input multi-band map from which we subtract tallest peak
    ;       freqs   - frequencies of bands
    ;   specdens    - list of spectral densities of each band of noise
    ;   sigm        - sigma of profile to be subtracted (correspending to first band)
    ;   binwidth    - binwidth of signal map, in arccminutes
    ;   /bbody      - use bbody instead of R-J power law
    ;   /EQ_SIGM    - 0 [default] - 1/f dependence
    ;               - 1 - all equal
    ;   /emissivity - emissivity of body (default = 1.5)
    ;   real_pos    - known positions of beam (debugging purposes mostly)

    ; Output - Struct containing
    ;   xparam      - estimated x-coordinate of peak
    ;   yparam      - estimated y-coordinate of peak
    ;   Aest        - estimated height of peaks in first band
    ;   map         - signal minus peak
    ;   sigm_x0     - theoretical deviation of xparam
    ;   sigm_a      - theoretical deviation of Aest
    ;   dropmask    - which bands dropped in calculation
    ;
    ; subtracts maximum peak of width sigm from map signal containing noise specdens

freqs = sig.freqs
signal = sig.signal
temp = max(freqs, max_freq); get max frequency
if ~keyword_set(normband) then normband = 0 ; normalize to leading band by default

; parameters
fitRange = 3        ; fit an area of +- fitRange * sigm around max to find peak
range = double(sqrt(n_elements(signal[*,*,0])))
num_bands = n_elements(freqs)
; set up default values
if emissivity eq !NULL then emissivity = 1.5D else emissivity = double(emissivity)
if tdust eq !NULL then tdust = 40D else tdust = double(tdust)

; mask out DC bin
mask = make_array(range, range, value=1D)
mask[(range - 1) / 2, (range - 1) / 2] = 0 ; mask out DC bin

; generate amps, sigm
if keyword_set(emissivity) then begin
    if keyword_set(tdust) then begin
        if ~ keyword_set(bbody) then amps = amps_multi(1, freqs, emissivity, tdust, normband=normband) else amps = amps_multi(1, freqs, emissivity, tdust, bbody=bbody, normband=normband) ; whatever we get out of the estimator is just multiplied by amps for return
    endif else begin
        if ~ keyword_set(bbody) then amps = amps_multi(1, freqs, emissivity, normband=normband) else amps = amps_multi(1, freqs, emissivity, bbody=bbody, normband=normband) ; whatever we get out of the estimator is just multiplied by amps for return
    endelse
endif else begin
    if keyword_set(tdust) then begin
        if ~ keyword_set(bbody) then amps = amps_multi(1, freqs, tdust, normband=normband) else amps = amps_multi(1, freqs, tdust, bbody=bbody, normband=normband) ; whatever we get out of the estimator is just multiplied by amps for return
    endif else begin
        if ~ keyword_set(bbody) then amps = amps_multi(1, freqs, normband=normband) else amps = amps_multi(1, freqs, bbody=bbody, normband=normband) ; whatever we get out of the estimator is just multiplied by amps for return
    endelse
endelse
if ~ keyword_set(eq_sigm) then sigms = sigm_multi(sigm, freqs) else sigms = sigm_multi(sigm, freqs, eq_sigm=eq_sigm)

; first, convolve
temp_sigm_a = [] ; sigm_a per band, up to binwidth etc. factors
convSig = [] ; list of convolutions
weights = dblarr(num_bands) ; compute the weights for later part in here
for i=0, num_bands - 1 do begin
    kSignal = fft_shift(fft(signal[*,*,i]))
    paddedKernel = addgauss(amps[i], sigms[i], 0, 0, dblarr(range, range))
    kkernel = fft_shift(fft(paddedKernel))

    ; convolve and figure out position
    CSTemp = real_part(fft(fft_shift(mask * conj(kkernel) * kSignal / specDens[*,*,i], /REVERSE), /INVERSE))
    CSNorm = max(real_part(fft(fft_shift(mask * conj(kkernel) * kkernel / specDens[*,*,i], /REVERSE), /INVERSE)))
    convSig = [[[convSig]], [[CSTemp/CSNorm]]]

    temp_sigm_a = [temp_sigm_a, CSNorm]

    ; weights (1 / sigm^2)
    weights[i] = REAL_PART(TOTAL( (fft_shift(dist(range, 1)^2) # replicate(1.0 / (binwidth * range)^2, range)) * conj(kkernel) * kkernel / specDens[*,*,i] )) * (2 * !PI *  (binwidth * range))^2
endfor

dropMask = make_array(num_bands, value=1)
; if both positions given, use
if keyword_set(real_x) and keyword_set(real_y) then begin
    xparam = real_x
    yparam = real_y
    weights = dblarr(num_bands) ; compute the weights for later part
    for i=0, num_bands - 1 do begin
        paddedKernel = addgauss(amps[i], sigms[i], 0, 0, dblarr(range, range))
        kkernel = fft_shift(fft(paddedKernel))
        ; weights (1 / sigm^2)
        weights[i] = REAL_PART(TOTAL( (fft_shift(dist(range, 1)^2) # replicate(1.0 / (binwidth * range)^2, range)) * conj(kkernel) * kkernel / specDens[*,*,i] )) * (2 * !PI *  (binwidth * range))^2
    endfor
endif else begin ; if not given position then compute it
    ; rough heuristic for A
    aest_temp = 0
    for i=0,5 do begin
        aest_temp += max(convSig[*,*,i]) * temp_sigm_a[i]
    endfor
    aest_temp /= total(temp_sigm_a)

    ; fit total convolved map (highest SNR)
    convSigTot = dblarr(range, range)
    for i=0,5 do convSigTot += convSig[*,*,i]
    temp = max(convSigTot, loc)
    loc = to2d(range, loc)
    xcent = (loc[0] + range) MOD range ; forces positive index
    ycent = (loc[1] + range) MOD range

    ; compute weighted average/sigm_x0
    xparams = dblarr(num_bands)
    yparams = dblarr(num_bands)
    for i=0, num_bands - 1 do begin
        ; ; fit peak on convolved map
        ; ; sigm_fit = FIX(1/(sqrt(weights[i]) * aest_temp * binwidth)) + 1 ; sigm_pixels for this band rounded up, roughly
        ; minx = max([xcent - fitRange, 0]) ; bounds checking
        ; maxx = min([xcent + fitRange, range - 1])
        ; miny = max([ycent - fitRange, 0])
        ; maxy = min([ycent + fitRange, range - 1])
        ; ; minx = max([xcent - fitRange * sigm_fit, 0]) ; bounds checking
        ; ; maxx = min([xcent + fitRange * sigm_fit, range - 1])
        ; ; miny = max([ycent - fitRange * sigm_fit, 0])
        ; ; maxy = min([ycent + fitRange * sigm_fit, range - 1])
        ; submap = convSig[minx:maxx, miny:maxy, i]
        ; params = fitquad(alog(submap))

        ; ; extract centers; backwards from expected because IDL x,y axis are flipped
        ; xtemp = params[3] / (2 * params[1]) + fix(xcent - fitRange)
        ; ytemp = params[2] / (2 * params[1]) + fix(ycent - fitRange)
        ;     ; last term compensates for non-centered kernel
        ; if xtemp eq xtemp and ytemp eq ytemp then begin
        ;     xparams[i] = (xtemp + range) mod range
        ;     yparams[i] = (ytemp + range) mod range
        ; endif else begin
        ;     xparams[i] = 0
        ;     yparams[i] = 0
        ; endelse
        ret = subtractmax(sig.signal[*,*,i], specDens[*,*,i], sigm, binwidth)
        xparams[i] = ret.xparam
        yparams[i] = ret.yparam
    end

    ; make estimates for center, avoiding xparams = 0 weights
    xparam = TOTAL(weights * xparams) / TOTAL(weights[where(xparams)])
    yparam = TOTAL(weights * yparams) / TOTAL(weights[where(xparams)])

    ; 05/16/16 -- Change: re-estimate until deviations from xparam, yparam within 5sigma
    chi2_pos = TOTAL(weights * (xparams - xcent)^2 + weights * (yparams - ycent)^2) * binwidth^2 * aest_temp^2
    stop
    while chi2_pos gt chisqr_cvf(0.01, 10) do begin
        temp = max([[xparams - xcent]^2 * dropMask * weights, [yparams - ycent]^2 * dropMask  * weights], loc)
        dropMask[loc] = 0
        chi2_pos = TOTAL(weights * dropMask * (xparams - xcent)^2 + $
            weights * dropMask * (yparams - ycent)^2) * $
            binwidth^2 * aest_temp^2
    endwhile

    ; if only one param is set, use it
    if keyword_set(real_x) then xparam=real_x
    if keyword_set(real_y) then yparam=real_y
endelse


; compute uncertainties and amplitude estimator
Aesttop = 0
Aestbottom = 0
inv_sigm_a = 0

; amplitude estimator + uncertainty
kfilters = []
for i=0, num_bands - 1 do begin
    ; redo filter for each band
    filter = addgauss(amps[i], sigms[i], xparam, yparam, dblarr(range, range)) ; start with gaussian in correct place
    kFilter = fft_shift(fft(filter))
    kSignal = fft_shift(fft(signal[*,*,i]))

    ; compute top, bottom of aest
    Aesttop += real_part(TOTAL(conj(kFilter) * kSignal / specDens[*,*,i] * mask))
    Aestbottom += real_part(TOTAL(conj(kFilter) * kFilter / specDens[*,*,i] * mask))

    ; sigm_a
    inv_sigm_a += REAL_PART( (binwidth * range)^2 * TOTAL(kfilter * CONJ(kfilter) / specDens[*,*,i] * mask))
    kfilters = [[[kfilters]], [[kfilter]]]
ENDFOR

; if amplitude estimator is set, use it, else use computed one
; we still need the above loop anyways to compute a, beta uncertainties
if keyword_set(real_amp) then begin
    aest = real_amp
endif else begin
    aest = aesttop / aestbottom
endelse

if ~keyword_set(bbody) then retknown = subtractknown_multi(sig, specDens, sigm, binwidth, emissivity, real_pos=[xparam, yparam], real_amp=Aest, normband=normband)$
    else retknown = subtractknown_multi(sig, specDens, sigm, binwidth, emissivity, tdust, real_pos=[xparam, yparam], real_amp=Aest, /bbody, normband=normband)
signal = retknown.sig.signal

; residual/chi2/derivative dbeta
dbeta = 0 ; derivative dbeta
dt_dust = 0 ; derivative dt_dust
inv_sigm_T = 0
inv_covarbt = 0
inv_covarat = 0
inv_sigm_beta = 0
inv_covar = 0
for i=0, num_bands - 1 do begin
    ; compute beta derivatives
    inv_sigm_beta += REAL_PART( (binwidth * range)^2 * aest^2 * TOTAL(kfilters[*,*,i] * CONJ(kfilters[*,*,i]) / specDens[*,*,i] * mask)) * alog(freqs[i] / 1500)^2
    inv_covar += REAL_PART( (binwidth * range)^2 * aest * TOTAL(kfilters[*,*,i] * CONJ(kfilters[*,*,i]) / specDens[*,*,i] * mask)) * alog(freqs[i] / 1500)
    dbeta += REAL_PART((binwidth * range )^2 * TOTAL((-conj(fft_shift(fft(signal[*,*,i]))) * aest * kfilters[*,*,i]) / specdens[*,*,i] * mask) * 2 * alog(freqs[i] / 1500))

    ; compute tdust derivatives
    if keyword_set(bbody) then begin
        ; recall 0.04799 is h/k_B
        dt_dust += 2 * (binwidth * range )^2 * TOTAL(- (conj(fft_shift(fft(signal[*,*,i]))) * aest * kfilters[*,*,i]) /$
            specdens[*,*,i] * mask) * 0.04799 * exp(0.04799 * freqs[i] / tdust) * freqs[i] / ( tdust^2 * (exp(0.04799 * freqs[i] / tdust) - 1))
        inv_sigm_T += (binwidth * range )^2 * aest^2 * TOTAL(abs(kfilters[*,*,i])^2 / specdens[*,*,i] * mask) * (0.04799 * exp(0.04799 * freqs[i] / tdust) * freqs[i] / ( tdust^2 * (exp(0.04799 * freqs[i] / tdust) - 1)))^2
        inv_covarbt += REAL_PART( (binwidth * range)^2 * aest^2 * TOTAL(kfilters[*,*,i] * CONJ(kfilters[*,*,i]) / specDens[*,*,i] * mask)) * alog(freqs[i] / 1500) * (0.04799 * exp(0.04799 * freqs[i] / tdust) * freqs[i] / ( tdust^2 * (exp(0.04799 * freqs[i] / tdust) - 1)))
        inv_covarat += REAL_PART( (binwidth * range)^2 * aest * TOTAL(kfilters[*,*,i] * CONJ(kfilters[*,*,i]) / specDens[*,*,i] * mask)) * (0.04799 * exp(0.04799 * freqs[i] / tdust) * freqs[i] / ( tdust^2 * (exp(0.04799 * freqs[i] / tdust) - 1)))
    endif
endfor

; to actually get covariances, must invert Hessian
Hess = [[ inv_sigm_a, inv_covar, inv_covarat], [inv_covar, inv_sigm_beta, inv_covarbt], [inv_covarat, inv_covarbt, inv_sigm_T]]
; Hess = [[inv_sigm_beta, inv_covarbt], [inv_covarbt, inv_sigm_T]]
covar = invert(Hess)

return, {xparam:xparam,$
    yparam:yparam,$
    Aest: Aest,$
    sig:retknown.sig,$
    chi2:retknown.chi2,$
    sigm_x0:sqrt(1/TOTAL(weights * dropMask)) / (aest * binwidth),$
    ; sigm_a: sqrt(covar[0,0]),$
    sigm_a: sqrt(1 / inv_sigm_a),$
    sigm_beta:sqrt(covar[1,1]),$
    sigm_tdust:sqrt(covar[2,2]),$
    covar_betaamp:0,$
    covar_bt:covar[0,1],$
    covar_at:covar[0,2],$
    emissivity:emissivity,$
    tdust:tdust,$
    dbeta:dbeta,$
    dt_dust:real_part(dt_dust),$
    convSig:convSig,$
    dropMask:dropMask $
    }
    ; uncertainty in pixels, not arcmins
end
