function subtractmat, signal, specdens, foundX, foundY, sigm, binwidth
compile_opt idl2, HIDDEN
    ; Arguments
    ;       signal - real signal from which to subtract kernel once, at max convolution value
    ;       specdens - spectral density of noise
    ;       foundX - x-locations of pulses to examine
    ;       foundY - y-locations of pulses to examine
    ;       sigm - sigma of kernel
    ;       binwidth - binwidth

    ; Return
    ;       Struct containing
    ;           struct.amps - amplitudes of pulses
    ;           struct.sigm - theoretical sigma of pulses
    ;
    ; matrix subtraction for multiple pulses, given foundX, foundY

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; CONSTANTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; get k-space forms of stuff
range = LONG(sqrt(N_ELEMENTS(signal)))
numBeams = N_ELEMENTS(foundX)
ksignal = fft_shift(fft(signal))

; allocate arrays
M = dcomplexarr(numBeams, numBeams)
S = dcomplexarr(numBeams)
Fouriers = []

; let's do stupidest implementation for now

; get Fourier Transforms
for i=0, (numBeams - 1) do begin
    Fouriers = [[Fouriers], [reform(fft_shift(fft(addgauss(1, sigm, foundX[i], foundY[i], dblarr(range, range)))), range * range)]]
    ; access with reform(Fouriers[*,i], [range, range])
endfor

for i=0, (numBeams - 1) do begin
    for j=0, (numBeams - 1) do begin
        M[i,j] = TOTAL(conj(reform(Fouriers[*,i], [range, range])) * reform(Fouriers[*,j], [range, range]) / specDens)
    endfor
    S[i] = TOTAL(conj(reform(Fouriers[*,i], [range, range])) * ksignal / specDens)
endfor

; calculate amps!
amps = invert(M) # S

return, {amps:real_part(amps), x:foundX, y:foundY}

END
