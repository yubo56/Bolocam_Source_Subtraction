function calc_vars, emissivity, tdust, amp
; computes variances
; currently only does beta-tdust
; in_f is flux in highest band

; definitions
range = 256D
binwidth = 480 / range ; units of arcsec, total width 480 arcsec
freqs = [400, 352.94, 272.73, 230.77, 150] ; frequencies in GHz
num_bands = n_elements(freqs)
sigm = 15D / (2 * sqrt(2 * alog(2))) ; sigm in arcsec, in terms of FWHM
sigm /= binwidth ; sigm in bins
xparam = 127.7
yparam = 127.3
emis = 1.5

; setup
amps = amps_multi(1, freqs, emissivity, tdust, /bbody)
sigms = sigm_multi(sigm, freqs)
specdens = spdensgen_multi(range, $
    [0.181, 0.137, 0.112, 0.0947, 0.049] / 2, $
    replicate(0,num_bands), replicate(8.0/3, num_bands), 1) ; PSD in units of mJy^2

inv_sigm_beta = 0
inv_covarbt = 0
inv_sigm_T = 0
inv_sigm_a = 0
inv_covar = 0
inv_covarAT = 0

for i=0, num_bands - 1 do begin
    mask = make_array(range, range, value=1D)
    mask[(range - 1) / 2, (range - 1) / 2] = 0 ; mask out DC bin

    ; compute filter, = nu/nu0^beta * nu^3 / (exp(...) - 1) * s
    filter = addgauss(amps[i], sigms[i], xparam, yparam, dblarr(range, range)) 
    kFilter = fft_shift(fft(filter))

    inv_sigm_beta += REAL_PART( (range)^2 * TOTAL(kfilter * CONJ(kfilter) / specDens[*,*,i] * mask)) * alog(freqs[i] / 1500)^2
    inv_sigm_T += real_part((range )^2 * TOTAL(kfilter * CONJ(kfilter) / specdens[*,*,i] * mask)) * (0.04799 * exp(0.04799 * freqs[i] / tdust) * freqs[i] / ( tdust^2 * (exp(0.04799 * freqs[i] / tdust) - 1)))^2

    ; 0.04799 = hbar/k_B
    inv_covarbt += REAL_PART( (range)^2 * TOTAL(kfilter * CONJ(kfilter) / specDens[*,*,i] * mask)) * alog(freqs[i] / 1500) * (0.04799 * exp(0.04799 * freqs[i] / tdust) * freqs[i] / ( tdust^2 * (exp(0.04799 * freqs[i] / tdust) - 1)))
    inv_covarAT += REAL_PART( (range)^2 * TOTAL(kfilter * CONJ(kfilter) / specDens[*,*,i] * mask)) * (0.04799 * exp(0.04799 * freqs[i] / tdust) * freqs[i] / ( tdust^2 * (exp(0.04799 * freqs[i] / tdust) - 1)))

    inv_sigm_a += REAL_PART( (range)^2 * TOTAL(kfilter * CONJ(kfilter) / specDens[*,*,i] * mask))
    inv_covar += REAL_PART( (range)^2 * TOTAL(kfilter * CONJ(kfilter) / specDens[*,*,i] * mask)) * alog(freqs[i] / 1500)
endfor

covar_btmat = [[inv_sigm_beta * amp^2, inv_covarbt * amp^2], [inv_covarbt * amp^2, inv_sigm_T * amp^2]]
covar_btmat = invert(covar_btmat)

covar_mat = [[inv_sigm_beta * amp^2, inv_covarbt * amp^2, inv_covar * amp],$
      [inv_covarbt * amp^2, inv_sigm_T * amp^2, inv_covarAT * amp],$
          [inv_covar * amp, inv_covarAT * amp, inv_sigm_a]]
covar_mat = invert(covar_mat)

return, [covar_btmat[0,0], covar_btmat[1,1], covar_btmat[1,0]]
; return, [covar_mat[0,0], covar_mat[1,1], covar_mat[1,0]]
end
