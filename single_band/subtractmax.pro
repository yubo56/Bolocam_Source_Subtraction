function subtractmax, signal, specDens, sigm, binwidth
compile_opt idl2, HIDDEN
    ; Input
    ;   signal      - input map from which we subtract tallest peak
    ;   specdens    - spectral density of noise
    ;   sigm        - sigma of profile to be subtracted
    ;   binwidth    - binwidth of signal map, in arccminutes

    ; Output - Struct containing
    ;   xparam      - estimated x-coordinate of peak
    ;   yparam      - estimated y-coordinate of peak
    ;   Aest        - estimated height of peak
    ;   map         - signal minus peak
    ;   sigm_x0     - theoretical deviation of xparam
    ;   sigm_y0     - theoretical deviation of yparam
    ;   sigm_a      - theoretical deviation of Aest

; subtracts maximum peak of width sigm from map signal containing noise specdens

fitRange = 3        ; fit an area of +- fitRange around max
range = double(sqrt(n_elements(signal)))

; convolve map
kSignal = fft_shift(fft(signal))
    ; no need to manually insert units, since taken care of with L^2 correspondence
paddedKernel = addgauss(1, sigm, 0, 0, dblarr(range, range)) ; kernel is normalized
    ; no units
kkernel = fft_shift(fft(paddedKernel))
; convolve and figure out position
convSig = fft(fft_shift(conj(kkernel) * kSignal / specDens, /REVERSE), /INVERSE)
     ; optimal filter = kernel/specDens

; fit map
temp = max(convSig, loc)
loc = to2d(range, loc)
xmax = (loc[0] + range) MOD range ; forces positive index
ymax = (loc[1] + range) MOD range

; fit peak on convolved map
minx = max([xmax - fitRange, 0]) ; bounds checking
maxx = min([xmax + fitRange, range - 1])
miny = max([ymax - fitRange, 0])
maxy = min([ymax + fitRange, range - 1])
; xparam = xmax
; yparam = ymax
submap = convSig[minx:maxx, miny:maxy]
params = fitquad(alog(submap))

; extract centers; backwards from expected because IDL x,y axis are flipped
xparam = params[3] / (2 * params[1]) + xmax - fitRange
xparam = (xparam + range) mod range
    ; last term compensates for non-centered kernel
yparam = params[2] / (2 * params[1]) + ymax - fitRange
yparam = (yparam + range) mod range
filter = addgauss(1, sigm, xparam, yparam, dblarr(range, range)) ; start with gaussian in correct place
kFilter = fft_shift(fft(filter))
    ; no units to insert; filter is unitless, and /arcmin^-2 comes with L substitution
Aest = real_part(TOTAL(conj(kFilter) * kSignal / specDens) / TOTAL(conj(kFilter) * kFilter / specDens))

; compute expected sigmas
sigm_x0 = SQRT( REAL_PART(1 / (TOTAL( (fft_shift(dist(range, 1)^2) # replicate(1.0 / (binwidth * range)^2, range)) * conj(kfilter) * kfilter / specDens )))) / (2 * !PI * Aest * (binwidth * range))
sigm_y0 = SQRT( REAL_PART(1 / (TOTAL( (replicate(1.0 / (binwidth * range)^2, range) # fft_shift(dist(range, 1)^2)) * conj(kfilter) * kfilter / specDens )))) / (2 * !PI * Aest * (binwidth * range))
    ; replicate contains a 1/(binwidth * range)^2 b/c nu_x needs units of 1/(binwidth * range), 1/arcmin
sigm_A = REAL_PART(SQRT(1  / ( (binwidth * range)^2 * TOTAL(kfilter * CONJ(kfilter) / specDens))))

return, {xparam:xparam, yparam:yparam, Aest: Aest, signal: (signal - Aest * filter), sigm_x0:sigm_x0 / binwidth, sigm_y0:sigm_y0 / binwidth, sigm_a:sigm_a} ; divide by binwidth so that uncertainties are in pixels, not arcmins
end
