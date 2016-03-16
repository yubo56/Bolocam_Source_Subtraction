function subtractmax_multi, sig, specDens, sigm, binwidth, emissivity, tdust, bbody=bbody, eq_sigm=eq_sigm, real_x=real_x, real_y=real_y, real_amp=real_amp
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
    ;   real_pos    - known positions of beam (debugging purposes mostly

    ; Output - Struct containing
    ;   xparam      - estimated x-coordinate of peak
    ;   yparam      - estimated y-coordinate of peak
    ;   Aest        - estimated height of peaks in first band
    ;   map         - signal minus peak
    ;   sigm_x0     - theoretical deviation of xparam
    ;   sigm_a      - theoretical deviation of Aest
    ;
    ; subtracts maximum peak of width sigm from map signal containing noise specdens

freqs = sig.freqs
signal = sig.signal
temp = max(freqs, max_freq); get max frequency

; parameters
fitRange = 3        ; fit an area of +- fitRange around max to find peak
range = double(sqrt(n_elements(signal[*,*,0])))
num_bands = n_elements(freqs)
; set up default values
if emissivity eq !NULL then emissivity = 1.5D else emissivity = double(emissivity)
if tdust eq !NULL then tdust = 40D else tdust = double(tdust)

; generate amps, sigm
if keyword_set(emissivity) then begin
    if keyword_set(tdust) then begin
        if ~ keyword_set(bbody) then amps = amps_multi(1, freqs, emissivity, tdust) else amps = amps_multi(1, freqs, emissivity, tdust, bbody=bbody) ; whatever we get out of the estimator is just multiplied by amps for return
    endif else begin
        if ~ keyword_set(bbody) then amps = amps_multi(1, freqs, emissivity) else amps = amps_multi(1, freqs, emissivity, bbody=bbody) ; whatever we get out of the estimator is just multiplied by amps for return
    endelse
endif else begin
    if keyword_set(tdust) then begin
        if ~ keyword_set(bbody) then amps = amps_multi(1, freqs, tdust) else amps = amps_multi(1, freqs, tdust, bbody=bbody) ; whatever we get out of the estimator is just multiplied by amps for return
    endif else begin
        if ~ keyword_set(bbody) then amps = amps_multi(1, freqs) else amps = amps_multi(1, freqs, bbody=bbody) ; whatever we get out of the estimator is just multiplied by amps for return
    endelse
endelse
if ~ keyword_set(eq_sigm) then sigms = sigm_multi(sigm, freqs) else sigms = sigm_multi(sigm, freqs, eq_sigm=eq_sigm)

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
    ; convolve maps
    convSig = [] ; list of convolutions
    weights = dblarr(num_bands) ; compute the weights for later part in here
    for i=0, num_bands - 1 do begin
        kSignal = fft_shift(fft(signal[*,*,i]))
        paddedKernel = addgauss(amps[i], sigms[i], 0, 0, dblarr(range, range))
        kkernel = fft_shift(fft(paddedKernel))

        ; convolve and figure out position
        convSig = [[[convSig]], [[real_part(fft(fft_shift(conj(kkernel) * kSignal / specDens[i], /REVERSE), /INVERSE))]]]

        ; weights (1 / sigm^2)
        weights[i] = REAL_PART(TOTAL( (fft_shift(dist(range, 1)^2) # replicate(1.0 / (binwidth * range)^2, range)) * conj(kkernel) * kkernel / specDens[*,*,i] )) * (2 * !PI *  (binwidth * range))^2
    endfor

    ; fit max-frequency map (brightest source)
    temp = max(convSig[*, *, max_freq], loc)
    loc = to2d(range, loc)
    xmax = (loc[0] + range) MOD range ; forces positive index
    ymax = (loc[1] + range) MOD range

    ; fit peak on convolved map
    minx = max([xmax - fitRange, 0]) ; bounds checking
    maxx = min([xmax + fitRange, range - 1])
    miny = max([ymax - fitRange, 0])
    maxy = min([ymax + fitRange, range - 1])
    submap = convSig[minx:maxx, miny:maxy, max_freq]
    submap[where(submap le 0, /null)] = 1e-8 ; sets minimum to be 1e-8, since we only use this for position determination anyways
    params = fitquad(alog(submap))

    ; extract centers; backwards from expected because IDL x,y axis are flipped
    xcent = params[3] / (2 * params[1]) + xmax - fitRange
    xcent = (xcent + range) mod range
        ; last term compensates for non-centered kernel
    ycent = params[2] / (2 * params[1]) + ymax - fitRange
    ycent = (ycent + range) mod range

    ; compute weighted average/sigm_x0
    xparams = dblarr(num_bands)
    yparams = dblarr(num_bands)
    for i=0, num_bands - 1 do begin
        if i eq max_freq then begin; if it is the max_map, we already have its poesition estimate
            xparams[i] = xcent
            yparams[i] = ycent
        endif else begin
            ; fit peak on convolved map
            minx = max([xcent - fitRange, 0]) ; bounds checking
            maxx = min([xcent + fitRange, range - 1])
            miny = max([ycent - fitRange, 0])
            maxy = min([ycent + fitRange, range - 1])
            submap = convSig[minx:maxx, miny:maxy, i]
            params = fitquad(alog(submap))

            ; extract centers; backwards from expected because IDL x,y axis are flipped
            xtemp = params[3] / (2 * params[1]) + fix(xcent - fitRange)
            xparams[i] = (xtemp + range) mod range
                ; last term compensates for non-centered kernel
            ytemp = params[2] / (2 * params[1]) + fix(ycent - fitRange)
            yparams[i] = (ytemp + range) mod range
        endelse
    end

    ; make estimates for center
    xparam = TOTAL(weights * xparams) / TOTAL(weights)
    yparam = TOTAL(weights * yparams) / TOTAL(weights)

    ; if only one param is set, use it
    if keyword_set(real_x) then xparam=real_x
    if keyword_set(real_y) then yparam=real_y
endelse


; mask out DC bin
mask = make_array(range, range, value=1D)
mask[(range - 1) / 2, (range - 1) / 2] = 0 ; mask out DC bin

; compute uncertainties and amplitude estimator
Aesttop = 0
Aestbottom = 0
inv_sigm_a = 0
inv_sigm_beta = 0
inv_covar = 0

; aest/sigm_a/sigm_beta
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
    inv_sigm_a += REAL_PART( (range)^2 * TOTAL(kfilter * CONJ(kfilter) / specDens[*,*,i] * mask))
    kfilters = [[[kfilters]], [[kfilter]]]

    ; compute uncertaintes for both black body and power law case
    inv_sigm_beta += REAL_PART( (range)^2 * TOTAL(kfilter * CONJ(kfilter) / specDens[*,*,i] * mask)) * alog(freqs[i] / 1500)^2
    inv_covar += REAL_PART( (range)^2 * TOTAL(kfilter * CONJ(kfilter) / specDens[*,*,i] * mask)) * alog(freqs[i] / 1500)
ENDFOR

; if amplitude estimator is set, use it, else use computed one
; we still need the above loop anyways to compute a, beta uncertainties
if keyword_set(real_amp) then begin
    aest = real_amp
endif else begin
    aest = aesttop / aestbottom
endelse

; residual/chi2/derivative dbeta
dbeta = 0 ; derivative dbeta
dt_dust = 0 ; derivative dt_dust
inv_sigm_T = 0
inv_covarbt = 0
for i=0, num_bands - 1 do begin
    ; compute residual
    ; compute derivatives for both black body and power law case
    ; note that signal contains residual now
    if keyword_set(bbody) then begin
        ; recall 0.04799 is h/k_B
        dt_dust += (range )^2 * TOTAL((conj(fft_shift(fft(signal[*,*,i]))) * kfilters[*,*,i] * $
            0.04799 * exp(0.04799 * freqs[i] / tdust) * freqs[i] / ( tdust^2 * (exp(0.04799 * freqs[i] / tdust) - 1)^2)) / specdens[*,*,i] * mask)
        inv_sigm_T += (range )^2 * TOTAL(abs(kfilters[*,*,i])^2 * (0.04799 * exp(0.04799 * freqs[i] / tdust) * freqs[i] / ( tdust^2 * (exp(0.04799 * freqs[i] / tdust) - 1)))^2 / specdens[*,*,i] * mask)
        inv_covarbt += REAL_PART( (range)^2 * TOTAL(kfilters[*,*,i] * CONJ(kfilters[*,*,i]) / specDens[*,*,i] * mask)) * alog(freqs[i] / 1500) * (0.04799 * exp(0.04799 * freqs[i] / tdust) * freqs[i] / ( tdust^2 * (exp(0.04799 * freqs[i] / tdust) - 1)))
    endif
    dbeta += (range )^2 * TOTAL((-conj(fft_shift(fft(signal[*,*,i]))) * aest * kfilters[*,*,i]) / specdens[*,*,i] * mask) * 2 * alog(freqs[i] / 1500)
endfor

; to actually get covariances, must invert non-diagonal subspace of covariance matrix
covar = [[ inv_sigm_a, aest * inv_covar], [aest * inv_covar, aest^2 * inv_sigm_beta]]
covar = invert(covar)
covar_bt = [[ inv_sigm_beta, inv_covarbt], [inv_covarbt, inv_sigm_T]]
covar_bt = invert(covar_bt)

if ~keyword_set(bbody) then retknown = subtractknown_multi(sig, specDens, sigm, binwidth, emissivity, real_pos=[xparam, yparam], real_amp=Aest)$
    else retknown = subtractknown_multi(sig, specDens, sigm, binwidth, emissivity, tdust, real_pos=[xparam, yparam], real_amp=Aest, /bbody)

return, {xparam:xparam,$
    yparam:yparam,$
    Aest: Aest,$
    sig:retknown.sig,$
    sigm_x0:sqrt(1/TOTAL(weights)) / (aest * binwidth),$
    sigm_a: sqrt(covar[0,0]),$
    chi2:retknown.chi2,$
    sigm_beta:sqrt(covar[1,1]),$
    sigm_tdust:sqrt(1/covar_bt[1,1]) / aest,$
    covar_betaamp:mean([covar[0,1], covar[1,0]]),$
    covar_bt:covar_bt[0,1] / aest,$
    emissivity:emissivity,$
    tdust:tdust,$
    dbeta:real_part(dbeta),$
    dt_dust:real_part(dt_dust)}
    ; uncertainty in pixels, not arcmins
end
