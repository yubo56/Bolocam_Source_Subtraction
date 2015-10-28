function submatwrap_multi, in_sig, specdens, sigm_gauss, binwidth, NUM_PULSES=num_pulses, no_opt=no_opt
compile_opt idl2, HIDDEN
    ; simple wrapper to feed iterimp output to subtractmat


; iterimp
if keyword_set(no_opt) then retval = iterimp_multi(in_sig, specdens, sigm_gauss, binwidth, /SIGM, NUM_PULSES=num_pulses, /no_opt) else retval = iterimp_multi(in_sig, specdens, sigm_gauss, binwidth, /SIGM, NUM_PULSES=num_pulses)

; matrix compute
if retval.as ne -1 then begin
    retvalmat = subtractmat_multi(in_sig, specdens, retval.xs, retval.ys, sigm_gauss, binwidth)

    ; subtract
    range = long(sqrt(n_elements(in_sig.signal[*,*,0])))
    num_bands = n_elements(in_sig.freqs)
    subtract = addgauss_multi(retvalmat.amps, replicate(sigm_gauss, n_elements(retvalmat.amps)), retval.xs, retval.ys, {signal:dblarr(range, range, num_bands), freqs:in_sig.freqs})
    res = in_sig.signal - subtract.signal

    ; compute chi^2, DOF
;   DOF = range^2 - 3 * n_elements(retvalmat.amps) - 1 ; minus 1 b/c specdens[DC] = 0
;   mask = make_array(range, range, value=1D)
;   mask[(range - 1) / 2, (range - 1) / 2] = 0 ; mask out DC bin
;   chi2 = (range * binwidth)^2 * TOTAL(abs(fft_shift(fft(res)))^2 / (specdens * binwidth) * mask)

    ; return
    return, {xs: retval.xs,$
        ys:retval.ys,$
        As:retvalmat.amps,$
        sigm_x0:retval.sigm_x0,$
        sigm_a:retval.sigm_a,$
;       DOF:DOF,$
;       chi2:chi2,$
        signal: {signal:res, freqs:in_sig.freqs}}
endif else return, {numIters: -1, converged:1, signal:in_sig}

end
