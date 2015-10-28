function submatwrap, orig, specdens, sigm_gauss, binwidth, NUM_PULSES=num_pulses, no_opt=no_opt
compile_opt idl2, HIDDEN
    ; simple wrapper to feed iterimp output to subtractmat


; iterimp
if keyword_set(no_opt) then retval = iterimp(orig, specdens, sigm_gauss, binwidth, /SIGM, NUM_PULSES=num_pulses, /no_opt) else retval = iterimp(orig, specdens, sigm_gauss, binwidth, /SIGM, NUM_PULSES=num_pulses)

; matrix compute
if retval.numiters ne -1 then begin
    retvalmat = subtractmat(orig, specdens, retval.xs, retval.ys, sigm_gauss, binwidth)

    ; subtract
    range = long(sqrt(n_elements(orig)))
    subtract = addgauss(retvalmat.amps, replicate(sigm_gauss, n_elements(retvalmat.amps)), retval.xs, retval.ys, dblarr(range, range))
    res = orig - subtract

    ; compute chi^2, DOF
    DOF = range^2 - 3 * n_elements(retvalmat.amps) - 1 ; minus 1 b/c specdens[DC] = 0
    mask = make_array(range, range, value=1D)
    mask[(range - 1) / 2, (range - 1) / 2] = 0 ; mask out DC bin
    chi2 = (range * binwidth)^2 * TOTAL(abs(fft_shift(fft(res)))^2 / (specdens) * mask)

    ; return
    return, {xs: retval.xs,$
        ys:retval.ys,$
        As:retvalmat.amps,$
        sigm_x0:retval.sigm_x0,$
        sigm_y0:retval.sigm_y0,$
        sigm_a:retval.sigm_a,$
        numIters: retval.numIters,$
        converged:retval.converged,$
        DOF:DOF,$
        chi2:chi2,$
        signal: res}
endif else return, {numIters: -1, converged:1, signal:orig}

end
