function subtractknown_multi, emissivity, dp, signal=sig, specdens=specdens, amp=amp, sigm=sigm, xparam=xparam, yparam=yparam
; return chi2 for known amp, sigm, xparam, yparam, emissivity
; power law usage only (for now)

freqs = sig.freqs
signal = sig.signal
range = double(sqrt(n_elements(signal[*,*,0])))
num_bands = n_elements(freqs)

amps = amps_multi(amp, freqs, emissivity=emissivity)
sigms = sigm_multi(sigm, freqs)

chi2=0
dbeta=0
for i=0, num_bands - 1 do begin
    ; compute residual
    kernel = addgauss(amps[i], sigms[i], xparam, yparam, dblarr(range, range))
    signal[*,*,i] = signal[*,*,i] - kernel

    ; mask out DC bin
    mask = make_array(range, range, value=1D)
    mask[(range - 1) / 2, (range - 1) / 2] = 0 ; mask out DC bin
    chi2 += (range )^2 * TOTAL(abs(fft_shift(fft(signal[*,*,i])))^2 / (specdens[*,*,i]) * mask)

    dbeta += real_part((range )^2 * TOTAL((-conj(fft_shift(fft(signal[*,*,i]))) * fft_shift(fft(kernel))) / specdens[*,*,i] * mask) * 2 * alog(freqs[i] / freqs[0]))
endfor

dp = dbeta
print, dp, chi2

return, chi2
end
