function spdensgen_multi, range, whites, elbows, pows, binwidth
compile_opt idl2, HIDDEN
; Inputs
;   whites - list of white noise components of specdens, in Jy (will square here)
;   elbows - TODO - for now, just constant multiple * power term in specdens
;   pows - power laws for various specdens bands
;   binwidth - usual binwidth
; Output
;   specdens - [range, range, m] array, with m = n_elements(whites) = n_elements(elbows) = n_elements(pows)
;
; Given parameters, makes a multi-band specdens

num_bands = n_elements(whites)

; generate specdens
specdens = []
for i=0, (num_bands - 1) do begin
    specdens = [[[specdens]], [[spdensgen(range, whites[i], elbows[i], pows[i], binwidth)]]]
endfor

return, specdens
end
