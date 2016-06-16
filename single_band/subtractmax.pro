function subtractmax, signal, specDens, sigm, binwidth, real_pos=real_pos, range=range
compile_opt idl2, HIDDEN
    ; Input
    ;   signal      - input map from which we subtract tallest peak
    ;   specdens    - spectral density of noise
    ;   sigm        - sigma of profile to be subtracted
    ;   binwidth    - binwidth of signal map, in arccminutes
    ;   real_pos    - actual position of source (opt)

    ; Output - Struct containing
    ;   xparam      - estimated x-coordinate of peak
    ;   yparam      - estimated y-coordinate of peak
    ;   Aest        - estimated height of peak
    ;   map         - signal minus peak
    ;   sigm_x0     - theoretical deviation of xparam
    ;   sigm_y0     - theoretical deviation of yparam
    ;   sigm_a      - theoretical deviation of Aest

; subtracts maximum peak of width sigm from map signal containing noise specdens
kSignal = fft_shift(fft(signal))
    ; no need to manually insert units, since taken care of with L^2 correspondence
range = double(sqrt(n_elements(signal)))

; mask out DC bin
mask = make_array(range, range, value=1D)
mask[(range - 1) / 2, (range - 1) / 2] = 0 ; mask out DC bin

if ~keyword_set(real_pos) then begin
    fitRange = 3        ; fit an area of +- fitRange around max

    ; convolve map
    paddedKernel = addgauss(1, sigm, 0, 0, dblarr(range, range)) ; kernel is normalized
        ; no units
    kkernel = fft_shift(fft(paddedKernel))
    ; convolve and figure out position
    convSig = fft(fft_shift(mask * conj(kkernel) * kSignal / specDens, /REVERSE), /INVERSE)
         ; optimal filter = kernel/specDens

    ; fit map
    if ~keyword_set(range) then temp = max(convSig, loc) else$
        temp = max(convSig[range[0]:range[1], range[2]:range[3]], loc)
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
endif else begin
    xparam = real_pos[0]
    yparam = real_pos[1]
end

filter = addgauss(1, sigm, xparam, yparam, dblarr(range, range)) ; start with gaussian in correct place
kFilter = fft_shift(fft(filter))
    ; no units to insert; filter is unitless, and /arcmin^-2 comes with L substitution
Aest = real_part(TOTAL(mask * conj(kFilter) * kSignal / specDens) / TOTAL(mask * conj(kFilter) * kFilter / specDens))

; compute expected sigmas
sigm_x0 = SQRT( REAL_PART(1 / (TOTAL( (fft_shift(dist(range, 1)^2) # replicate(1.0 / (binwidth * range)^2, range)) * mask * conj(kfilter) * kfilter / specDens )))) / (2 * !PI * Aest * (binwidth * range))
sigm_y0 = SQRT( REAL_PART(1 / (TOTAL( (replicate(1.0 / (binwidth * range)^2, range) # fft_shift(dist(range, 1)^2)) * mask * conj(kfilter) * kfilter / specDens )))) / (2 * !PI * Aest * (binwidth * range))
    ; replicate contains a 1/(binwidth * range)^2 b/c nu_x needs units of 1/(binwidth * range), 1/arcmin
sigm_A = REAL_PART(SQRT(1  / ( (binwidth * range)^2 * TOTAL(kfilter * CONJ(kfilter) / specDens))))

return, {xparam:xparam,$
    yparam:yparam,$
    Aest: Aest,$
    signal: (signal - Aest * filter),$
    sigm_x0:sigm_x0 / binwidth,$
    sigm_y0:sigm_y0 / binwidth,$
    sigm_a:sigm_a} ; divide by binwidth so that uncertainties are in pixels, not arcmins
end
