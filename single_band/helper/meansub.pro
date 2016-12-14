function meansub, signal, specdens, sigm_gauss, binwidth
compile_opt idl2, HIDDEN
    ; Input -
    ;   signal - what from which we subtract
    ;   specdens - spectral density
    ;   sigm_gauss - width of Gaussians
    ;   binwidth - width of pixels
    ;
    ; Output -
    ;   signal, As, xs, ys - residual map; amplitudes, positions of beams
    ;
    ; Wrapper around submatwrap that also goes back to look at estimates to check whether double incidence


range = double(sqrt(n_elements(signal)))
fitwidth = 3 * sigm_gauss ; width around which to search for coincident peaks

; get threshold
threshold = 10 * specdens[(range - 1)/2, (range - 1)/2] / (fitwidth)^2 ; prolly eventually want this in terms of the residual RMS in terms of sigmas, derived in 2D.pdf
minheight = 0.3 ; minimum height before reject

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;routine;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
newxs = []
newys = []
newAs = []

; initial estimates
retval = submatwrap(signal, specdens, sigm_gauss, binwidth)
if retval.as[0] eq -1 then begin
    return, {signal:signal, numIters:-1}
    ; whaaaaaaaaaaaaaat why would this ever happen
endif
orig = retval.signal
; fix 'em up!
for i=0, (n_elements(retval.xs) - 1) do begin
    num_pulses = 1

    ; set up
    sig = orig[retval.xs[i] - fitwidth:retval.xs[i] + fitwidth, retval.ys[i] - fitwidth:retval.ys[i] + fitwidth]
    xs = [retval.xs[i] - FIX(retval.xs[i] - fitwidth)]
    ys = [retval.ys[i] - FIX(retval.ys[i] - fitwidth)]
    As = [retval.As[i]]
    ; if mean too large
    while mean(sig) ge threshold do begin
        ; add signal back
        temp =addgauss(As, replicate(sigm_gauss, num_pulses), xs, ys, sig)

        ; resubtract
        num_pulses++
        retvaltemp = submatwrap(temp, specdens, sigm_gauss, binwidth, num_pulses = num_pulses)

        ; check whether we beat minimum height
        if min(retvaltemp.as) le minheight then break ; break without overwriting As.

        ; update values
        sig = retvaltemp.signal
        xs = retvaltemp.xs
        ys = retvaltemp.ys
        As = retvaltemp.As
    endwhile

    ; store changes, update map
    newxs = [newxs, xs + FIX(retval.xs[i] - fitwidth)]
    newys = [newys, ys + FIX(retval.ys[i] - fitwidth)]
    newAs = [newAs, As]
    orig[retval.xs[i] - fitwidth:retval.xs[i] + fitwidth, retval.ys[i] - fitwidth:retval.ys[i] + fitwidth] = sig
endfor

; reoptimize now that everything is fleshed out
optret = estopt(orig, specdens, newAs, newxs, newys, sigm_gauss, binwidth)
orig = optret.signal
newxs = optret.xs
newys = optret.ys
newAs = optret.As
retvalmat = subtractmat(orig, specdens, newxs, newys, sigm_gauss, binwidth)
subtract = addgauss(retvalmat.amps, replicate(sigm_gauss, n_elements(retvalmat.amps)), newxs, newys, dblarr(range, range))


return, {signal:orig - subtract, xs:newxs, ys:newys, As:newAs, numIters:1}

end
