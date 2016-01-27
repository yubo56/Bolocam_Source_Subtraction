function addgauss_multi, in_a, in_sigm, x0, y0, in_sig, emissivity, tdust, bbody=bbody, EQ_SIGM = eq_sigm
compile_opt idl2, HIDDEN
; Inputs
;   in_a - amplitude of source in first-frequency band (first in the in_signal order)
;   in_sigm - width of source in first-frequency band
;   x0 - location of source
;   y0 - location of source
;   in_sig - struct containing
;       signal - [N, N, num_bands] array containing sources
;       freqs - [num_bands] array containing frequencies; should correspond to ordering of in_signal
;   /EQ_SIGM    - 0 [default] - 1/f dependence
;               - 1 - all equal
;   /emissivity - emissivity (default = 3.0)
;   /tdust      - t_dust (default = 40K)
;   /bbody      - use black body spectrum instead of power law
; Outputs
;   signal - [N, N, num_bands] array containing inserted sources
;
; multi-band wrapper around addgauss.pro


; for safety
freqs = double(in_sig.freqs)
signal = in_sig.signal
size_amps = size(in_a)
num_amps = n_elements(in_a)

; get num_bands
num_bands = n_elements(freqs)

if size_amps[0] eq 0 then begin ; if just single beam
    ; generate amps, sigms
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
    ; sigms
    if ~ keyword_set(eq_sigm) then sigms = sigm_multi(in_sigm, freqs) else sigms = sigm_multi(in_sigm, freqs, eq_sigm=eq_sigm)


    ; add in each band
    for i=0, num_bands - 1 do begin
        signal[*,*,i] = addgauss(amps[i], sigms[i], x0, y0, signal[*, *, i])
    endfor
endif else begin ; loop over input parameters
    for i=0, num_amps - 1 do begin
        ; generate amps, sigms
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
        ; sigms
        if ~ keyword_set(eq_sigm) then sigms = sigm_multi(in_sigm[i], freqs) else sigms = sigm_multi(in_sigm[i], freqs, eq_sigm=eq_sigm)


        ; add in each band
        for j=0, num_bands - 1 do begin
            signal[*,*,j] = addgauss(amps[j], sigms[j], x0[i], y0[i], signal[*, *, j])
        endfor
    endfor
endelse


return, {signal:signal, freqs:freqs}
end
