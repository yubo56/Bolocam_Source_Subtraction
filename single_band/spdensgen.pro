function spdensgen, range, const, constpow, power, binwidth
compile_opt idl2, HIDDEN
    ; takes various parameters and returns a spectral density in 2D

range = long(range) ; just in case, often issues since FIX is freaking tiny

k = fft_shift(dist(range))
k[(range - 1) / 2, (range - 1) / 2 ] = 1
specDens = make_array(range, range, VALUE=const^2) + (constPow * const^2 / k ^ power)
    ; use white noise for now
    ; units of Jy^2
specDens *= (binwidth * range)^2 ; units of Jy^2 / arcmin^-2
specDens /= range^2 ; missing 1/N^2 division that is never performed since we generate specdens directly rather than from fft

return, specDens
end
