; add ourselves to path
!PATH = expand_path('+' + '~/SURF2015/routines') + path_sep(/search_path) + !PATH
!PATH = expand_path('+' + '~/bolocam_svn') + path_sep(/search_path) + !PATH

; various bolocam_svn needs
astrolib
def_user_common
def_elec_common
def_bolo_common
def_phys_common
def_astro_common
def_cosmo_common

; compile routines (maintain this manually)
; .compile addfunc addgauss estopt fitquad gennoise meansub submatwrap subtractmat subtractmax to2d spdensgen iterimp help_yubo hist_yubo

; multi band versions
; .compile amps_multi sigm_multi gennoise_multi addgauss_multi addfunc_multi subtractmax_multi subtractmat_multi spdensgen_multi
; .compile estopt_multi iterimp_multi submatwrap_multi

; laziness ^.^
; single
; specdens = spdensgen(256, 0.1, 3, 1, 1)

; multiple
range = 256
freqs = [400, 352.94, 272.73, 230.77, 150] ; frequencies in GHz
num_bands = n_elements(freqs)

specdens0 = spdensgen_multi(range, $
    [0.181, 0.137, 0.112, 0.0947, 0.049] / 2, $
    replicate(3,num_bands), replicate(8.0/3, num_bands), 1) ; PSD in units of mJy^2
noise = gennoise_multi(specdens0, 1, freqs)
sig = addgauss_multi(1, 5, range/2 + 0.5, range/2 + 0.5, noise)
sig_bbody = addgauss_multi(1, 5, range/2 + 0.5, range/2 + 0.5, noise, /bbody)
; ret = subtractmax_multi(sig, specdens0, 5, 0.3)

; sig1 = addfunc_multi(noise, 1, 3, 5, minsep = 6) ; noise, height, sigm, num_gaussians
; ret1 = submatwrap_multi(sig1.sig, specdens0, 5, 0.3)
