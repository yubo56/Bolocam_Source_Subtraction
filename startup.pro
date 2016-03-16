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

; laziness ^.^
; single
; specdens = spdensgen(256, 0.1, 3, 1, 1)

; multiple
range = 256D
binwidth = 480 / range
freqs = [400, 352.94, 272.73, 230.77, 150] ; frequencies in GHz
num_bands = n_elements(freqs)

specdens0 = spdensgen_multi(range, $
    [0.181, 0.137, 0.112, 0.0947, 0.049] / 2, $
    replicate(0,num_bands), replicate(8.0/3, num_bands), binwidth) ; PSD in units of 0.01 mJy^2
noise = gennoise_multi(specdens0, binwidth, freqs)
sigm = 15D / (2 * sqrt(2 * alog(2))) ; sigm in arcsec, in terms of FWHM
sigm /= binwidth ; sigm in bins
amp = sqrt((0.181 / 2) / (!PI * sigm^2)) * 5 ; ref 2D.tex for why this is 5SNR
sig_bbody = addgauss_multi(amp, sigm, range/2 + 0.5, range/2 + 0.5, noise, 1.7, 13, /bbody)
; ret = subtractmax_multi(sig, specdens0, sigm, binwidth)

; sig1 = addfunc_multi(noise, 1, 3, 5, minsep = 6) ; noise, height, sigm, num_gaussians
; ret1 = submatwrap_multi(sig1.sig, specdens0, 5, 0.3)
