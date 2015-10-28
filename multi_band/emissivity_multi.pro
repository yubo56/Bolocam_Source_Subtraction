function emissivity_multi, sig, specdens0, sigm, binwidth, init=init, real_pos=real_pos
compile_opt idl2, HIDDEN
; usual inputs
; computes emissivity by calling tnmin

if ~ keyword_set(init) then init = 1.5D else init = double(init) ; init value for tnmin search
    ; idk why but ensuring init is double is necessary
limits = [0,3] ; limits on emissivity; 1.5 +- 1.5

num_iters = 0
if keyword_set(real_pos) then begin
    emissivity = tnmin('tnmin_wrap', init, functargs={signal:sig, specdens:specdens0, sigm:sigm, binwidth:binwidth, real_pos:real_pos}, autoderivative=0, parinfo={limited:[1,1], limits:limits}, niter=num_iters, /quiet)
    ret = subtractmax_multi(sig, specdens0, sigm, binwidth, emissivity = emissivity, real_pos=real_pos)
    return, {xparam:ret.xparam, yparam:ret.yparam, Aest:ret.Aest, sig:ret.sig, sigm_x0:ret.sigm_x0, sigm_a:ret.sigm_a, chi2:ret.chi2, emissivity:emissivity, num_iters:num_iters, sigm_beta:ret.sigm_beta}
endif else begin
    emissivity = tnmin('tnmin_wrap', init, functargs={signal:sig, specdens:specdens0, sigm:sigm, binwidth:binwidth}, autoderivative=0, parinfo={limited:[1,1], limits:limits}, niter=num_iters, /quiet)
    ret = subtractmax_multi(sig, specdens0, sigm, binwidth, emissivity = emissivity)
    return, {xparam:ret.xparam, yparam:ret.yparam, Aest:ret.Aest, sig:ret.sig, sigm_x0:ret.sigm_x0, sigm_a:ret.sigm_a, chi2:ret.chi2, emissivity:emissivity, num_iters:num_iters, sigm_beta:ret.sigm_beta}
endelse

end
