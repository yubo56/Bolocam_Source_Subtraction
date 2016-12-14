# Yubo's Source Subtraction IDL routines

Goal is to be able to both generate sky maps and to source subtract them across
multiple frequency bands. Included code does just that. Run help\_yubo for help
with arguments of various functions.

Documentation/coding style will probably be revised as merging into general kSZ
simulation

## Organization

### Overall flow

The flow goes roughly as follows:
- Generate a PSD using `spdensgen`
- Generate a noise realisation using `gennoise`
- Compute a fit using `submatwrap` in the single band case and either `submatwrap_multi` in the multi band fixed-SED case or `emissivity_multi`/`bbody_multi` in the beta-T\_dust case (I never implemented `submatwrap` for either of those two since I focused on trying to break the degenerate chi2 first)

Analysis techniques include:
- Run best-estimator routine multiple times and create multi-dimensional histogram, take a Gaussian fit for covariance matrix
- For a fixed noise realization, compute the chi2 surface as a function of the various parameters, interpret as a likelihood function to obtain a covariance matrix
- Analytically compute, given a signal map, spectral density etc., what the expected covariance matrix should look like.

The covariance matricies for all three of these should agree assuming a quadratic chi2 surface, which does not hold for beta-T\_dust simultaneously unconstrained

Some examples of these routines in action can be found under my [Bitbucket](https://bitbucket.org/yubo56/)

### File list
- `help_yubo.pro`
- `helper`
    - `to2d.pro`
    - `sfit_p.pro`
    - `invert_matuncert.pro`
    - `getcovars.pro`
    - `hist_wrap.pro`
    - `contour_wrap.pro`
- `single_band`
    - `helper`
        - `subtractmat.pro`
        - `meansub.pro`
        - `estopt.pro`
        - `iterimp.pro`
    - `subtractmax.pro`
    - `addfunc.pro`
    - `submatwrap.pro`
    - `spdensgen.pro`
    - `fitquad.pro`
    - `addgauss.pro`
    - `gennoise.pro`
- `multi_band`
    - `helper`
         - `iterimp_multi.pro`
         - `estopt_multi.pro`
         - `tnmin_wrap.pro`
         - `subtractmat_multi.pro`
         - `amps_multi.pro`
         - `sigm_multi.pro`
    - `spdensgen_multi.pro`
    - `addfunc_multi.pro`
    - `getsedpts.pro`
    - `subtractknown_multi.pro`
    - `bbody_multi.pro`
    - `subtractmax_multi.pro`
    - `gennoise_multi.pro`
    - `emissivity_multi.pro`
    - `addgauss_multi.pro`
    - `submatwrap_multi.pro`
    - `plot_bbody_snu.pro`

### Quick introduction to each file
- `help_yubo.pro` (procedure)
    - Prints out a short documentation on most of the routines that I use. Slightly outdated now, this readme might be more useful.
- `helper`
    - `to2d.pro`
        - Short little routine to convert a 1D array index to a 2D array index
    - `sfit_p.pro`
        - Canned IDL `sfit.pro` routine. Copied to `sfit_p` so would not override built-in `sfit`, but do not remember why I made a copy. I think it was an `IDL_PATH` conflict
    - `invert_matuncert.pro`
        - Inverts a 2x2 matrix that has uncertainties in its elements. The return also has uncertainties attached.
    - `getcovars.pro`
        - For a fit to the chi2 surface, computes the covariance matrix obtained by treating the chi2 as a likelihood function
    - `hist_wrap.pro`
        - A wrapper around the build in idl `hist` routine that has some nice fonting and colors and a useful options concerning what values to write where. Automatically computes placement of text for both logarithmic and non-logarithmic scales.
    - `contour_wrap.pro`
        - A wrapper around the build in idl `contour` routine that has some nice fonting and colors and a useful options concerning where to draw the contours.
- `single_band` Routines for the single frequency band case. Most of these are described in 2D.pdf.
    - `helper`
        - `subtractmat.pro`
            - Given the X, Y positions of some Gaussians, find the set of amplitudes of each Gaussian that minimize the chi2. Also known as simultaneously fitting for Gaussian amplitudes given fixed positions.
        - `meansub.pro`
            - This is an artifact of an old hypothesis. With `MINSEP=0`, it is possible for two nearly-coincident sources to look like a single source. `meansub` computes a list of sources by running `submatwrap` (the best subtraction routine I have, discussed later), then examines each individual source to see whether the mean of the residual map decreases if we add another source.
            - This is a faulty hypothesis, as is discussed in my thesis, since the mean of the residual map is a meaningless number.
        - `estopt.pro`
            - This is in charge of the "iterative improvement" procedure described above. Given a set of flux estimates, it iteratively adds and re-subtracts them until positions converge
        - `iterimp.pro`
            - Calls `estopt` but also computes initial flux estimates upon which to improve
    - `subtractmax.pro`
        - For a given 1D signal, spectral density and binwidth, find the maximum signal for a Gaussian with width `sigm`. Can constrain `real_pos` if the position of the Gaussian is to be fixed (useful for certain exercises)
    - `addfunc.pro`
        - Add a bunch of Gaussians of a fixed height to a `noise` array. Can specify a `MINSEP` to make sure the Gaussians do not overlap too much
    - `submatwrap.pro`
        - Feeds iterimp (which iterates until position estimators converge) into subtractmat (which computes optimal flux estimators given positions.
    - `spdensgen.pro`
        - Given some parameters, computes a 1D spectral power density
    - `fitquad.pro`
        - Fits a 4x4 patch to a quadratic. Fits the log of a Gaussian.
    - `addgauss.pro`
        - Add a gaussian of width `sigm` and height `a` to a given signal array
    - `gennoise.pro`
        - Given a spectral density and a binwidth, generate a random noise realisation
- `multi_band` multi-band versions of the above routines. I will only elaborate on those routines that don't exist in the single-band versions. The formalism is developed mostly only in my thesis, as by this time I had grown very bad at maintaining 2D.pdf
    - `helper`
         - `iterimp_multi.pro`
         - `estopt_multi.pro`
         - `tnmin_wrap.pro`
            - Since with multiple frequency bands we must also fit for beta, T\_dust, we must use a gradient descent fitter. IDL has a built in `tnmin`, and `tnmin_wrap` is a wrapper that specializes `tnmin` for our specific application. Supports both beta-only fitting or beta+T\_dust (`/bbody` flag)
         - `subtractmat_multi.pro`
         - `amps_multi.pro`
            - Given a beta or a beta+T\_dust, compute amplitudes relative to lowest frequency band.
         - `sigm_multi.pro`
            - Compute beam widening relative to lowest frequency band due to diffraction and aperture size.
    - `spdensgen_multi.pro`
    - `addfunc_multi.pro`
    - `getsedpts.pro`
        - computes the single-band estimators + uncertainties for a multi-band signal. Used to compare the multi-band flux estimators to the single-band estimators.
    - `subtractknown_multi.pro`
        - Given a known SED, subtract it from the map. Used to subtract the true parameters and compare the chi2 to that obtained from subtracting the best fit parameters
    - `bbody_multi.pro`
        - calls `subtractmax_multi` and performs gradient descent to find the best-fit beta/T\_dust. Can fix either parameter or let both float (bad degeneracies...)
    - `subtractmax_multi.pro`
        - Should be noted this only works for a fixed beta/T\_dust
    - `gennoise_multi.pro`
    - `emissivity_multi.pro`
        - Same as bbody_multi but only performs gradient descent on beta
    - `addgauss_multi.pro`
    - `submatwrap_multi.pro`
    - `plot_bbody_snu.pro`
        - Plots SEDs for a list of emissitivities and T\_dusts. Useful to compare two SEDs that have similar chi2 values

### Random notes
- `knose` b/c used to generate noise in k-space, k-noise
- `sigm` b/c `sigma` is a reserved word in IDL
