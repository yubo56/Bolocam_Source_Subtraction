function getsedpts, sig, specdens, sigm, binwidth, real_pos
; gets the estimates of the fluxes in each frequency band and their uncertainties
; inputs:
;   sig, specdens, sigm, binwidth - the usual multi-band
;   real_pos - [xpos, ypos] position of the mean (using multi-band convolution method)
; outputs:
;   amps: amps in each band
;   sigms: sigms in each band

freqs = sig.freqs
signal = sig.signal
range = double(sqrt(n_elements(signal[*,*,0])))
num_bands = n_elements(freqs)
sigms = sigm_multi(sigm, freqs)

; store results
amps = dblarr(num_bands)
sigm_a = dblarr(num_bands)

for i=0, num_bands - 1 do begin
    ret = subtractmax(signal[*,*,i], specdens[*,*,i], sigms[i], binwidth, real_pos=real_pos)
    amps[i] = ret.aest
    sigm_a[i] = ret.sigm_a
endfor

return, {amps:amps, sigms:sigm_a}
end
