function sigm_multi, in_sigm, freqs, eq_sigm=eq_sigm
compile_opt idl2, HIDDEN
; Input
;   in_sigm - sigm in first band
;   freqs - frequencies of bands
;   /EQ_AMP - 0 [default] - f^(2 + epsilon) power law approx to black body
;           - 1 - all equal
;           - 2 - full black body spectrum    [ TODO implement ]
; Output
;   sigms - computed amplitudes corresponding to the input bands
;
; given leading sigm and frequencies, return rest of sigms

num_bands = n_elements(freqs)
sigms = dblarr(num_bands)
if ~ keyword_set(eq_sigm) || eq_sigm eq 0 then begin
    ; 1/f
    for i=0, num_bands - 1 do begin
        sigms[i] = in_sigm * freqs[0] / freqs[i]
    endfor
endif else begin ; constant sigm
    sigms = replicate(in_sigm, num_bands)
endelse

return, sigms
end
