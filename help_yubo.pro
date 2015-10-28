pro help_yubo
compile_opt idl2, HIDDEN
print
print, 'Most function headers for Yubo functions written SURF 2015'
print
print, 'addfunc, noise, gaussHeight, sigm, numGaussians, MINSEP=minsep'
print, '    Wrapper around addgauss for random Gaussians'
print, 'addgauss, a, sigm, x0, y0, signal'
print, '    Adds gaussians to map'
print, 'gennoise, specDens, binwidth'
print, '    Generates noise by specdens'
print, 'spdensgen, range, const, constpow, power, binwidth'
print, '    Generates specdens from parameters'
print, 'submatwrap, orig, specdens, sigm_gauss, binwidth, NUM_PULSES=num_pulses'
print, '    Wrapper that feeds iterimp amplitude estimates to subtractmat'
print, 'subtractmax, signal, specDens, sigm, binwidth'
print, '    Subtracts tallest pulse from map'
print
print, 'multi-band functions'
print
print, 'addfunc_multi, noise, gaussHeight, sigm, numGaussians, MINSEP=minsep, emissivity=emissivity'
print, '    Wrapper around addgauss_multi for random Gaussians with optional minimum separation (default 0)'
print, 'addgauss_multi, a, sigm, x0, y0, signal, emissivity=emissivity'
print, '    Adds gaussians with optional emissivity (default 1.5)'
print, 'gennoise_multi, specDens, binwidth, freqs'
print, '    Generates noise by specdens'
print, 'spdensgen_multi, range, whites, elbows, pows, binwidth'
print, '    Generates multi-band specdens from parameters'
print, 'emissivity_multi, sig, specdens0, sigm, binwidth, init=init'
print, '    Computes emissivity of maximum brightness source, with optional initial guess (default 1.5)'
print, 'subtractmax_multi, sig, specDens, sigm, binwidth, emissivity=emissivity'
print, '    Subtracts tallest pulse from multi-band map, with optional emissivity (default 1.5)'
print, 'submatwrap_multti, orig, specdens, sigm_gauss, binwidth, NUM_PULSES=num_pulses, no_opt=no_opt'
print, '    Wrapper that feeds iterimp amplitude estimates to subtractmat; optional num_pulses number of beams, /no_opt sets so no reoptimization of guesses'
print
end
