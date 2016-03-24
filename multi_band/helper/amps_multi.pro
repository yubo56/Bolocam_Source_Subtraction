function amps_multi, in_a, in_freqs, emissivity, tdust, bbody=bbody
compile_opt idl2, HIDDEN
; Input
;   in_a - amplitude in first band
;   in_freqs - frequencies of bands
;   /bbody - if set, computes full black body spectrum, else computes
;          - f^(2 + epsilon) power law approx to black body
;   /emissivity - emissivity of body (default = 3.0)
;   tdust - temperature of interstellar dust
; Output
;   amps - computed amplitudes corresponding to the input bands
;
; given leading amplitude and frequencies, return full set of ones

if emissivity eq !NULL then emissivity = 1.5 ; default emissivity
if tdust eq !NULL then tdust = 15 ; default tdust
freqs = double(in_freqs)
num_bands = n_elements(freqs)
amps = dblarr(num_bands)

if keyword_set(bbody) then begin
    ; black body
    amps = (freqs / freqs[0])^(3.0 + emissivity) / (exp(0.04799 * freqs / tdust) - 1) ; 4.799e-11 is h/k_B, but freqs are in GHz
    amps *= in_a ; normalize to first peak
endif else begin
    ; f^2
    amps = (freqs / freqs[0])^(2.0 + emissivity) * in_a ; power law means 2 + epsilon b/c 1/(exp - 1) cancels with one power of nu
endelse
return, amps
end
