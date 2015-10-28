FUNCTION genNoise_multi, specDens, binwidth, in_freqs, DUPLICATE=duplicate
compile_opt idl2, HIDDEN
; argument
;       specDens - list of spectral densities, [N, N, num_bands]
;       binwidth - arcmin/pixel edge
;       /DUPLICATE - take first specdens and copy it DUPLICATE times
; return
;       noise - Noise is ensured to be fully real, and we only return real part
;
; generates list of noises from a list of specdenses

; get number of bands
sizeSD = size(specdens)
if sizeSD[0] lt 3 then num_bands = 0 else num_bands = sizeSD[3] 

; store to return noise

; if duplicate is set, two cases
if keyword_set(duplicate) then begin
    if num_bands eq 0 then begin
        ; only passed a 2D specdens, so duplicate the 2D specdens
        noise = gennoise(specdens, binwidth)
        for i=1, duplicate - 1 do begin
            noise = [[[noise]], [[gennoise(specdens, binwidth)]]]
        endfor
    endif else begin
        ; passed a 3D specdens, take first one
        noise = gennoise(specdens[*,*,0], binwidth)
        for i=1, duplicate - 1 do begin
            noise = [[[noise]], [[gennoise(specdens[*,*,0], binwidth)]]]
        endfor
    endelse
endif else begin
    ; duplicate not set
    noise = gennoise(specdens[*,*,0], binwidth)
    for i=1, num_bands - 1 do begin
        noise = [[[noise]], [[gennoise(specdens[*,*,i], binwidth)]]]
    endfor
endelse

return, {signal:noise, freqs:in_freqs}
end
