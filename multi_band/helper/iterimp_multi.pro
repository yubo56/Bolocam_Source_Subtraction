function iterimp_multi, orig, specdens, sigm_gauss, binwidth, SIGM=sigm, NUM_PULSES=num_pulses, no_opt=no_opt
compile_opt idl2, HIDDEN
    ; Input
    ;   orig        - input map from which we subtract tallest peak
    ;   specdens    - spectral density of noise
    ;   sigm_gauss  - sigma of profile to be subtracted
    ;   binwidth    - binwidth of signal map, in arccminutes
    ;   /SIGM       - whether or not to store uncertainty values

    ; Output - Struct containing
    ;   xs          - estimated x-coordinates of peaks
    ;   ys          - estimated y-coordinates of peaks
    ;   As          - estimated heights of peaks
    ;
    ;   (if /SIGM is set, struct also contains)
    ;   sigm_x0     - theoretical deviations of xs
    ;   sigm_a      - theoretical deviations of As

    ; Iteratively improves subtractmax's estimates until position estimates converge to within sigma/factor
    ; Also identifies all pulses above factor * noise without knowledge of how many pulses

SNR             = 4.5 ; considers a peak found when above SNR * specDens(DC bin) - USEFUL ONLY WHEN NUM_PULSES NOT SET

; relevant data are stored here
xs = []
ys = []
As = []

; store sigmas only if SIGM keyword is set
if keyword_set(SIGM) then begin
    sigm_x0 = []
    sigm_a = []
endif

range = long(sqrt(n_elements(orig.signal[*,*,0])))
signal = orig ; make a local copy

; initial estimates; do differently if know how many pulses
if keyword_set(num_pulses) then begin
    for i = 1, num_pulses do begin
        ; subtract out maximum pulse
        retVal = subtractmax_multi(signal, specdens, sigm_gauss, binwidth)

        ; when know how many pulses, don't cutoff at SNR
        xs = [xs, retVal.xparam]
        ys = [ys, retVal.yparam]
        As = [As, retval.Aest]
        ; store sigmas only if SIGM keyword is set
        if keyword_set(SIGM) then begin
            sigm_x0 = [sigm_x0, retval.sigm_x0]
            sigm_a = [sigm_a, retval.sigm_a]
        endif
        signal = retval.signal
    endfor
endif else begin
    while 1 do begin
        ; subtract out maximum pulse
        retVal = subtractmax_multi(signal, specdens, sigm_gauss, binwidth)
        if retVal.aest ge SNR * retval.sigm_a then begin
            xs = [xs, retVal.xparam]
            ys = [ys, retVal.yparam]
            As = [As, retval.Aest]
            ; store sigmas only if SIGM keyword is set
            if keyword_set(SIGM) then begin
                sigm_x0 = [sigm_x0, retval.sigm_x0]
                sigm_a = [sigm_a, retval.sigm_a]
            endif
            signal = retval.sig
        endif else begin
            break
        endelse
    endwhile
endelse

; call optimization function only if there are amplitudes AND optimization requested
if keyword_set(no_opt) then begin
    if keyword_set(SIGM) then return, {xs: xs, ys:ys, As:As, sigm_x0:mean(sigm_x0), sigm_a:mean(sigm_a), signal:signal} else return, {xs: xs, ys:ys, As:As, signal:signal}
endif else begin
    if (n_elements(As) gt 0) then optval = estopt_multi(signal, specdens, As, xs, ys, sigm_gauss, binwidth) else return, {As:[-1], signal:orig}
    if keyword_set(SIGM) then return, {xs: optval.xs, ys:optval.ys, As:optval.As, sigm_x0:mean(sigm_x0), sigm_a:mean(sigm_a), signal:optval.signal} else return, {xs: optval.xs, ys:optval.ys, As:optval.As, signal:optval.signal}
endelse
end
