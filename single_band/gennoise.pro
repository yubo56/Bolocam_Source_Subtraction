FUNCTION genNoise, specDens, binwidth
compile_opt idl2, HIDDEN
; argument
;       specDens - Spectral density of noise
; return
;       noise - Noise is ensured to be fully real, and we only return real part
; generates real noise given a spectral density
;
; assumes specDens is passed in the sane way, [-N/2 + 1, N/2]
; for real noise, re(knoise) has to have the following characteristics:
;   1) Inside square [-N/2+1:N/2-1]^2 has to be symmetric, reverse(reverse()) leaves invariant
;   2) Outside bins [N/2, -N/2+1: N/2-1] and the corresponding have to be symmetric
;   3) High bin [N/2, N/2] is irrelevant
;
; and for imaginary, just antisymmetric, 3) must vanish



range = long(sqrt(N_ELEMENTS(specDens)))
sigm = sqrt(specDens * (binwidth * range)^2) / sqrt( 2 ) ; specdens = Jy^2/arcmin^-2, sigma = Jy/arcmin^-2
    ; recall J = v^2/L^2

temp = sigm * RANDOMN(SEED, range, range)

if (range MOD 2 eq 0) then begin ; even cases are different b/c set unpaired to zero

    knoseRe = dblarr(range, range) ; we can only reflect part of it, not including highest bin
    knoseRe[0,0] = (temp[0:range - 2,0:range - 2] + reverse(reverse(temp[0:range - 2, 0:range - 2], 1), 2)) / sqrt(2)
    knoseRe[range - 1, 0:range - 2] = (temp[range - 1, 0:range - 2] + reverse(temp[range - 1, 0:range - 2], 2)) / sqrt(2)
    knoseRe[0:range - 2, range - 1] = (temp[0:range - 2, range - 1] + reverse(temp[0:range - 2, range - 1], 1)) / sqrt(2)
    knoseRe[range - 1, range - 1] = temp[range - 1, range - 1]; copy high bin
    ; we drop the outermost elements b/c otherwise reverse(reverse) goes bonkers

    knoseIm = dblarr(range, range)
    knoseIm[0,0] = ((temp[0: range - 2, 0:range - 2] - reverse(reverse(temp[0: range-2, 0:range-2], 1), 2)) / sqrt(2))
    knoseIm[range - 1, 0:range - 2] = (temp[range - 1, 0:range - 2] - reverse(temp[range - 1, 0:range - 2], 2)) / sqrt(2)
    knoseIm[0:range - 2, range - 1] = (temp[0:range - 2, range - 1] - reverse(temp[0:range - 2, range - 1], 1)) / sqrt(2)

    ; guarantees Re is symmetric, Im is antisymmetric about -vec-k (reverse(reverse()) is just -vec-k)
    ; such tricks are why I am not a Caltech professor yet
endif else begin
    knoseRe = (temp + reverse(reverse(temp),2)) / sqrt(2)
    knoseIm = (temp - reverse(reverse(temp),2)) / sqrt(2)
endelse

knose = COMPLEX(knoseRe, knoseIm); generate full knoise
    ; units Jy/arcmin^-2 still

; perform the FFT
noise = FFT(fft_shift(knose,/REVERSE), /INVERSE) / (binwidth * range)^2 ; units Jy
RETURN, real_part(noise)
END
