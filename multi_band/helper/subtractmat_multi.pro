function subtractmat_multi, in_sig, specdens, foundX, foundY, sigm, binwidth, emissivity, tdust, bbody=bbody, eq_sigm=eq_sigm
compile_opt idl2, HIDDEN
    ; Arguments
    ;       in_sig - real signal in multi bands, usual struct containing .signal, .freqs components
    ;       specdens - spectral density of noisE
    ;       foundX - x-locations of pulses to examine
    ;       foundY - y-locations of pulses to examine
    ;       sigm - sigma of kernel
    ;       binwidth - binwidth
    ;       eq_amp, eq_sigm - see amps_multi
    ;       emissivity - emissivity of body (default = 3.0)
    ;       tdust - temperature of interstellar dust

    ; Return
    ;       Struct containing
    ;           struct.amps - amplitudes of beams in first band
    ;           struct.x, y - same as input foundX, foundY (for convenience)
    ;
    ; matrix subtraction for multiple pulses, given foundX, foundY

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; CONSTANTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; get k-space forms of stuff
freqs = in_sig.freqs
signal = in_sig.signal
num_bands = n_elements(freqs)
range = double(sqrt(N_ELEMENTS(signal[*,*,0])))
numBeams = N_ELEMENTS(foundX)

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

; allocate arrays
M = dcomplexarr(numBeams, numBeams)
S = dcomplexarr(numBeams)
fouriers = []

; let's do stupidest implementation for now

; get Fourier Transforms
for j=0, (numBeams - 1) do begin
    sub_fouriers = []
    for i=0, num_bands - 1 do begin
        sub_fouriers = [[[sub_fouriers]], [reform(fft_shift(fft(addgauss(amps[i], sigms[i], foundX[j], foundY[j], dblarr(range, range)))), range * range)]]
    endfor
    fouriers = [[[fouriers]], [[sub_fouriers]]]
    ; access with reform(Fouriers[*,i,j], [range, range])
    ;   i is the i-th band, j is the j-th beam
endfor

for i=0, num_bands - 1 do begin
    for j=0, (numBeams - 1) do begin
        for k=0, (numBeams - 1) do begin
            M[j,k] += TOTAL(conj(reform(Fouriers[*, i, j], [range, range])) * reform(Fouriers[*, i, k], [range, range]) / specDens[*,*,i])
        endfor
        ksignal = fft_shift(fft(signal[*,*,i]))
        S[j] += TOTAL(conj(reform(Fouriers[*,i,j], [range, range])) * ksignal / specDens[*,*,i])
    endfor
endfor

print, M, S
; calculate amps!
amps = invert(M) # S

return, {amps:real_part(amps), x:foundX, y:foundY}

END
