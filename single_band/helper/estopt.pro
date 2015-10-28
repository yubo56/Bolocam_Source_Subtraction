function estopt, signal, specdens, As, xs, ys, sigm_gauss, binwidth
compile_opt idl2, HIDDEN
    ; Input
    ;   signal  - map after subtracting As, xs, ys
    ;   As      - height of pulses
    ;   xs      - x position of pulses
    ;   ys      - y position of pulses
    ;   binwidth- bin width as usual
    ;
    ; Return - Struct containing
    ;   signal  - res after re-optimising
    ;   As, xs, ys - beams after re-optimising
    ;   numIters - number of iterations for convergence
    ;   converged - whether converged
    ;
    ; Heavy lifter for iterimp, re-optimises a set of estimates
    ; also called after scrutinizing for double coincidences

cnvrgnceFact    = 20 ; considers a peak's position determined when iteration changes position by less than sigma_x/cnvrgnceFact
numConsec       = 10 ; number of consecutive subtractions that must change by less than sigm_x/cnvrgnceFact
maxcnv          = 500 ; maximum number of iterations before giving up

; set up flags
numPeaks = n_elements(As)
flags = replicate(numConsec, numPeaks)

numIters = 1
converged = 1
; amps = As
; exes = Xs
; whys = Ys
; iteratively improve estimates
while TOTAL(flags) gt 0 do begin ; while at least one flag needs to be examined
    for i=0,(N_ELEMENTS(As) - 1) do begin
        ; only if flagged for examination
        if flags[i] ne 0 then begin
            ; perform subtraction
            signal = addgauss(As[i], sigm_gauss, xs[i], ys[i], signal) ; add signal back
            retVal = subtractmax(signal, specdens, sigm_gauss, binwidth)
            signal = retval.signal

            ; flag appropriately
            if (abs(retVal.xparam - xs[i]) le retval.sigm_x0 / cnvrgnceFact) and (abs(retVal.yparam - ys[i]) le retval.sigm_y0 / cnvrgnceFact) then begin
                flags[i] -- ; if within convergence factor for positions, unflag
            endif else begin
                flags[i] = numConsec
            endelse

            ; always update values
            Xs[i] = retval.xparam
            Ys[i] = retval.yparam
            As[i] = retval.aest
            if keyword_set(SIGM) then begin
                sigm_x0[i] = retval.sigm_x0
                sigm_y0[i] = retval.sigm_y0
                sigm_a[i] = retval.sigm_a
            endif
        endif
    endfor
    ; amps = [[amps], [As]]
    ; exes = [[exes], [Xs]]
    ; whys = [[whys], [Ys]]
    if numIters++ ge maxcnv then begin
        ; print, 'Maximum convergence of' + strcompress(string(maxcnv)) + ' iterations reached, still' + strcompress(string(TOTAL(n_elements(flags)))) + ' peaks non-convergent.'
        converged = 0
        break
    endif
endwhile

; print, As, xs, ys
;       HANDY THING IS TO plot, exes[0,*], whys[0,*],xrange=[min(exes),max(exes)],yrange=[min(whys),max(whys)]
;                         oplot, exes[1,*], whys[1,*]
;                         etc.
; stop
return, {signal:signal, As:As, xs:xs, ys:ys, numIters:numIters, converged:converged}
end
