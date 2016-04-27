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
; sd = spdensgen(256, 0.1, 3, 1, 1)

; multiple
range = 256D
bw = 480 / range
emis = 1.7
td = 13
freqs = [400, 352.94, 272.73, 230.77, 150, 100] ; frequencies in GHz
num_bands = n_elements(freqs)

convrms = [0.181, 0.137, 0.112, 0.0947, 0.049, 0.029] / 2
s = 15D / (2 * sqrt(2 * alog(2))) ; sigm in arcsec, in terms of FWHM
s /= bw ; sigm in bins
a = convrms[0] * 5 / sqrt(!PI * s^2); SNR 5
sd = spdensgen_multi(range, $
    convrms, $ ; convert noise/beam -> noise/pixel
    replicate(0,num_bands), replicate(8.0/3, num_bands), bw) ; PSD in units of mJy^2
noise = gennoise_multi(sd, bw, freqs)
sig = addgauss_multi(a, s, range/2 + 0.5, range/2 + 0.5, noise, emis, td, /bbody, normband=0)
sed = getsedpts(sig, sd, s, bw, [127.5, 127.5])  
ret = subtractmax_multi(sig, sd, s, bw, /bbody, normband=0)

; sig1 = addfunc_multi(noise, 1, 3, 5, minsep = 6) ; noise, height, s, num_gaussians
; ret1 = submatwrap_multi(sig1.sig, sd, 5, 0.3)
