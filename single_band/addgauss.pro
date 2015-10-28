function addgauss, a, sigm, x0, y0, signal
compile_opt idl2, HIDDEN
    ; adds a Gaussian height a, sigma sigma at arbitrary offset $x0, y0$ into signal
    ; only process around 3sigma of Gaussian
    
range = LONG(sqrt(N_ELEMENTS(signal)))
sizeSig = size(a)
numSigs = N_ELEMENTS(a) ; assume N_ELEMENTS(A) = N_ELEMENTS(sigm, x0, y0)
a = double(a)
sigm = double(sigm)
x0 = double(x0)
y0 = double(y0)

retVal = signal ; store return

; if not list
if sizeSig[0] eq 0 then begin
    for x=round(x0 - 3 * sigm), round(x0 + 3 * sigm) do begin
        for y=round(y0 - 3 * sigm), round(y0 + 3 * sigm) do begin
            retVal[x mod range,y mod range] += a * EXP(-((x - x0)^2 + (y - y0)^2) / (2 * sigm^2))
        endfor
    endfor
endif else begin
    ; add signal
    for i=0, numSigs - 1 do begin
        for x=round(x0[i] - 3 * sigm[i]), round(x0[i] + 3 * sigm[i]) do begin
            for y=round(y0[i] - 3 * sigm[i]), round(y0[i] + 3 * sigm[i]) do begin
                retVal[x mod range,y mod range] += a[i] * EXP(-((x - x0[i])^2 + (y - y0[i])^2) / (2 * sigm[i]^2))
            endfor
        endfor
    endfor
endelse

return, retVal
end
