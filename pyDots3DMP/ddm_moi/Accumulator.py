import numpy as np
from scipy.stats import multivariate_normal as mvn, _mvn

import matplotlib.pyplot as plt
from matplotlib import animation

from dataclasses import dataclass, field
from typing import Optional, Union, Sequence

import concurrent.futures
from functools import partial

# from codetiming import Timer
    
class Accumulator:
    """
    2-D accumulator, calculated via method of images.
    
    Instantiate an object of the class with a list of drift rates, bound, and grid settings (time vector and grid vector).
    Drift rates can be single values (constant rate, or time series matching the length of tvec).
    Bound can be a single value (same for both accumulators, or separate for each one).
    The class then has a bunch of method calls associated for calculating the model predictions:
    - cdf returns a proportion of choices for the "positive" accumulator for each drift rate
        and the distribution of bound crossing times for each drift rate (independent of which accumulator hits first)
    - pdf returns the pdf of the accumulator model at each timepoint
        if full_pdf is False (default), it will return separate pdfs for each marginal (i.e. correct and errors).
        if full_pdf is True, it will return a single 3-D array of the square grid, with a 2-D pdf for each timepoint
    - dist runs the cdf method, and optionally the pdf method too (if return_pdf is true), with full_pdf set to False
    - log_posterior_odds uses the losing accumulator pdfs given correct and errors to calculate log odds of correct choice

    """

    def __init__(self, tvec: np.ndarray, grid_vec: np.ndarray, 
                 drift_rates: Optional[list] = None, 
                 bound: Optional[Union[float, Sequence, np.ndarray]] = 1.0,
                 num_images: Optional[int]=7,
                 wager_theta: Optional[float]=1.0):
        self.tvec = tvec
        self.grid_vec = grid_vec
        self.bound = bound
        self.drift_rates = drift_rates if drift_rates is not None else []
        self.num_images = num_images
        self.wager_theta = wager_theta

        if self.drift_rates:
            self.drift_labels = self.drift_rates.copy()
        
    @property
    def bound(self):
        """return symmetric bound as 2-element array"""
        b = self._bound
        if isinstance(b, (int, float)):
            b = [b, b]
        self._bound = np.array(b)
        return self._bound

    @bound.setter
    def bound(self, bound):
        self._bound = bound

    def apply_drifts(self, drifts, labels: Optional[list] = None, 
                     sensitivity: float = 1.0, urgency: Optional[Union[np.ndarray,float]] = None):
        """
        Set accumulator drift rates. Optionally add label for each drift.
        This also adds a mirrored drift rate for the anti-correlated accumulator, and 
        updates drift rates based on sensitivity and urgency parameters."""

        if isinstance(drifts, np.ndarray):
            drifts = np.split(drifts, drifts.shape[1], axis=1)
        
        # add corresponding negated value for anti-correlated accumulator
        # also update drift rates based on sensitivity and urgency, if provided

        if labels is not None:
            assert len(drifts) == len(labels), "drift rates and provided labels must match in length"
            self.drift_labels = labels
        
        for d, drift in enumerate(drifts):
            drift = drift * np.array([1, -1])
            drifts_posneg = _urgency_scaling(drift * sensitivity, self.tvec, urgency)
            self.drift_rates.append(drifts_posneg)
            

    def cdf(self):
        """calculate cdf at boundaries, and corresponding choices and times"""
        # faster to preallocate array rather than appending to list
        p_corr = np.zeros(len(self.drift_rates))
        rt_dist = np.zeros((len(self.drift_rates), len(self.tvec)))
        
        for d, drift in enumerate(self.drift_rates):
            p_corr[d], rt_dist[d, :], flux1, flux2 = _moi_cdf(
                self.tvec, drift, self.bound, 0.025, self.num_images)

        self.p_corr_ = p_corr
        self.rt_dist_ = rt_dist

    def pdf(self, full_pdf=False):

        # TODO allow flexible specification of grid_vec, to use mgrid
        if full_pdf:
            xmesh, ymesh = np.meshgrid(self.grid_vec, self.grid_vec)
        else:
            xmesh1, ymesh1 = np.meshgrid(self.grid_vec, self.grid_vec[-1])
            xmesh2, ymesh2 = np.meshgrid(self.grid_vec[-1], self.grid_vec)

        pdfs, marg_up, marg_lo = [], [], []

        for drift in self.drift_rates:

            if full_pdf:
                pdf_3d = _moi_pdf(xmesh, ymesh, self.tvec, drift,
                                  self.bound, self.num_images)
                pdfs.append(pdf_3d)

                pdf_up = pdf_3d[:, :, -1]  # right bound
                pdf_lo = pdf_3d[:, -1, :]  # top bound

            else:
                # only need to calculate pdf at the boundaries!

                # vectorized implementation is about 10x faster!
                pdf_lo = _moi_pdf_vec(xmesh1, ymesh1, self.tvec, drift,
                                      self.bound, self.num_images)
                pdf_up = _moi_pdf_vec(xmesh2, ymesh2, self.tvec, drift,
                                      self.bound, self.num_images)

            # distribution of losing accumulator, GIVEN winner has hit bound
            marg_up.append(pdf_up)  # right bound
            marg_lo.append(pdf_lo)  # top bound

        if full_pdf:
            self.pdf3D_ = np.stack(pdfs, axis=0)

        self.up_lose_pdf_ = np.stack(marg_up, axis=0)
        self.lo_lose_pdf_ = np.stack(marg_lo, axis=0)

        return self.up_lose_pdf_, self.lo_lose_pdf_

    def log_posterior_odds(self):
        """Return the log posterior odds given pdfs"""
        self.log_odds_ = log_odds(self.up_lose_pdf_, self.lo_lose_pdf_)
        return self.log_odds_

    @property
    def wager_map(self):
        return self.log_odds_ >= self.wager_theta

    def dv(self, drift, sigma):
        """Return accumulated DV for given drift rate and diffusion noise."""
        return _moi_dv(mu=drift*self.tvec.reshape(-1, 1),
                       s=sigma, num_images=self.num_images)

    def compute_distrs(self, return_pdf=False):
        """Calculate cdf and pdf for accumulator object. Returns self for chaining commands"""

        self.cdf()
        if return_pdf:
            self.pdf()
        return self

    def plot(self, d_ind: int = -1):
        """
        Plot summary of accumulator results.

        Parameters
        ----------
        d_ind : INT, optional
            index of which drift rate to plot. The default is the last one.

        Returns
        -------
        fig_cdf & fig_pdf: figure handles
        """
        if not self.is_fitted:
            return

        fig_cdf, axc = plt.subplots(2, 1, figsize=(4, 5))
        axc[0].plot(self.drift_labels, self.p_corr_)
        axc[0].set_xlabel('drift')
        axc[0].set_xticks(self.drift_labels)
        axc[0].set_ylabel('prob. correct choice')

        axc[1].plot(self.tvec, self.rt_dist_.T)
        # axc[1].legend(self.drift_labels, frameon=False)
        axc[1].set_xlabel('Time (s)')
        axc[1].set_title('RT distribution (no NDT)')
        fig_cdf.tight_layout()

        fig_pdf = None
        has_log_odds = hasattr(self, 'log_odds_')
        n = 3 if has_log_odds else 2
        if hasattr(self, 'up_lose_pdf_'):
            fig_pdf, axp = plt.subplots(n, 1, figsize=(5, 6))
            contour = axp[0].contourf(self.tvec, self.grid_vec,
                                      log_pmap(np.squeeze(self.up_lose_pdf_[d_ind, :, :])).T,
                                      levels=100)
            axp[1].contourf(self.tvec, self.grid_vec,
                            log_pmap(np.squeeze(self.lo_lose_pdf_[d_ind, :, :])).T,
                            levels=100)
            axp[0].set_title(f"Losing accumulator | Correct, drift rate {self.drift_labels[d_ind]}")
            axp[1].set_title(f"Losing accumulator | Error, drift rate {self.drift_labels[d_ind]}")
            cbar = fig_cdf.colorbar(contour, ax=axp[0])
            cbar = fig_cdf.colorbar(contour, ax=axp[1])

        if has_log_odds:
            vmin, vmax = 0, 3
            contour = axp[2].contourf(self.tvec, self.grid_vec,
                                        self.log_odds_.T, vmin=vmin, vmax=vmax,
                                        levels=100)
            axp[2].set_title("Log Odds of Correct Choice given Losing Accumulator")
            cbar = fig_pdf.colorbar(contour, ax=axp[2])
        fig_pdf.tight_layout()

        return fig_cdf, fig_pdf


    def plot_3d(self, d_ind=-1):
        raise NotImplementedError("3D plot is extremely slow or getting stuck somehow, need to improve")

        def animate_wrap(i):
            z = log_pmap(self.pdf3D[d_ind, i, :, :])
            cont = plt.contourf(self.grid_vec, self.grid_vec, z, levels=100)
            plt.title(f"Frame: {i + 1} - {self.tvec[i]:.2f}")
            return cont

        # def init():
        #     cont = plt.contourf(self.grid_vec, self.grid_vec, log_pmap(self.pdf3D[d_ind, 0, :, :]),
        #                         levels=100)
        #     return cont

        fig = plt.figure()

        anim = animation.FuncAnimation(fig, animate_wrap, frames=len(self.tvec))
        writervideo = animation.PillowWriter(fps=10)
        anim.save(f'pdf_animation_{self.drift_labels[d_ind]}.gif', writer=writervideo)


## ----------------------------------------------------------------
## % Private functions

# Python does not have explicit private/public functions, but by convention, private functions
# are prefaced with an underscore, meaning they are not advised to be called directly from outside
# the module, but exist only for internal use.

## ----------------------------------------------------------------

# TODO make these all static methods and move them inside the class namespace?

# @staticmethod
#@np_cache_minimal
def _sj_rot(j, s0, k):
    """
    Image rotation formalism.

    :param j: jth image
    :param s0: starting_point, length 2 array
    :param k: 2*k-1 is the number of images
    :return:
    """
    alpha = (k - 1)/k * np.pi
    sin_alpha = np.sin(j * alpha)
    sin_alpha_plus_k = np.sin(j * alpha + np.pi / k)
    sin_alpha_minus_k = np.sin(j * alpha - np.pi / k)

    if j % 2 == 0:
        s = np.array([[sin_alpha_plus_k, sin_alpha], [-sin_alpha, -sin_alpha_minus_k]])
    else:
        s = np.array([[sin_alpha, sin_alpha_minus_k], [-sin_alpha_plus_k, -sin_alpha]])

    return (1 / np.sin(np.pi / k)) * (s @ s0.T)

def _weightj(j, mu, sigma, sj, s0):
    """weight of the jth image"""
    return (-1) ** j * np.exp(mu @ np.linalg.inv(sigma) @ (sj - s0).T)

#@lru_cache(maxsize=32)
def _corr_num_images(num_images: int) -> tuple[np.ndarray, int]:
    """returns 2-D accumulator correlation given a number of images"""

    k = int(np.ceil(num_images / 2))
    rho = -np.cos(np.pi / k)
    sigma = np.array([[1, rho], [rho, 1]])

    return sigma, k

# TODO this has the same input and outputs as _moi_pdf, so should be a single function with a flag
def _moi_pdf_vec(
    xmesh: np.ndarray, 
    ymesh: np.ndarray, 
    tvec: np.ndarray,
    mu: np.ndarray,
    bound=np.array([1, 1]),
    num_images: int = 7):
    
    """
    Calculate 2-D pdf according to method of images (vectorized implementation).

    :param xmesh: x-values for pdf computation
    :param ymesh: y-valuues, should match shape of xmesh
    :param tvec: 1-D array containing times to evaluate pdf
    :param mu: drift rate 2xlen(tvec) array (to incorporate any urgency signal)
    :param bound: bound, length 2 array
    :param num_images: number of images for MOI, default is 7
    :return: 2-D probability density function evaluated at points in xmesh-ymesh
    """
    
    sigma, k = _corr_num_images(num_images)

    nx, ny = xmesh.shape
    pdf_result = np.zeros((len(tvec), nx, ny)).squeeze()
    # pdf_result2 = np.zeros_like(pdf_result)

    xy_mesh = np.dstack((xmesh, ymesh))
    new_mesh = xy_mesh.reshape(-1, 2)

    # xy_mesh is going to be an X*Y*2 mesh.
    # for vectorized pdf calculation, reshape it to N*2
    # TODO this might require more testing for the full_pdf

    s0 = -bound

    covs = tvec[:, None, None] * sigma
    mu_t = tvec[:, None] * mu

    pdf_result += np.exp(_multiple_logpdfs_vec_input(new_mesh, s0 + mu_t, covs))

    for j in range(1, k*2):
        sj = _sj_rot(j, s0, k)

        aj_all = np.zeros_like(tvec).reshape(-1, 1)

        # skip the first sample (t starts at 1)
        for t in range(1, len(tvec)):
            a_j = _weightj(j, mu[t, :].T, sigma, sj, s0)
            aj_all[t] = a_j

            # pdf_result[t, ...] += a_j * mvn(mean=sj + mu[t, :], cov=sigma*tvec[t]).pdf(xy_mesh)

        # use vectorized implementation!
        # TODO unit tests with np.allclose vs scipy mvn result
        pdf_result += (aj_all * np.exp(_multiple_logpdfs_vec_input(new_mesh, sj + mu_t, covs)))

    return pdf_result

def _moi_pdf(
    xmesh: np.ndarray, 
    ymesh: np.ndarray,
    tvec: np.ndarray,
    mu: np.ndarray,
    bound=np.array([1, 1]), 
    num_images: int = 7
    ) -> np.ndarray:
    """
    Calculate pdf according to method of images. Older implementation,
    this is not vectorized over time so runs slower.

    :param xmesh: x-values for pdf computation
    :param ymesh: y-valuues, should match shape of xmesh
    :param tvec: 1-D array containing times to evaluate pdf
    :param mu: drift rate 2xlen(tvec) array (to incorporate any urgency signal)
    :param bound: bound, length 2 array
    :param num_images: number of images for MOI, default is 7
    :return: pdf at each timepoint (t, x, y)-shape array
    """
    
    sigma, k = _corr_num_images(num_images)

    nx, ny = xmesh.shape
    pdf_result = np.zeros((len(tvec), nx, ny)).squeeze()

    xy_mesh = np.dstack((xmesh, ymesh))

    s0 = -bound

    # skip the first sample (t starts at 1)
    for t in range(1, len(tvec)):

        pdf_result[t, ...] = _pdf_at_timestep(
            tvec[t], mu[t, :], sigma, xy_mesh, k, s0)

    return pdf_result

def _pdf_at_timestep(
    t, 
    mu: np.ndarray, 
    sigma: np.ndarray, 
    xy_mesh: np.ndarray, 
    k: int, 
    s0: np.ndarray
    ):

    pdf = mvn(mean=s0 + mu*t, cov=sigma*t).pdf(xy_mesh)

    # j-values start at 1, go to k*2-1
    for j in range(1, k*2):
        sj = _sj_rot(j, s0, k)
        a_j = _weightj(j, mu.T, sigma, sj, s0)
        pdf += a_j * mvn(mean=sj + mu*t, cov=sigma*t).pdf(xy_mesh)

    return pdf

def _moi_cdf(
    tvec: np.ndarray, 
    mu: np.ndarray,
    bound = np.array([1, 1]),
    margin_width: float = 0.025,
    num_images: int = 7
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the cdf of a 2-D particle accumulator.

    Returns:
        a) the probability of a correct choice
        b) the distribution of response times (bound crossings)
    choices are calculated by evaluating cdf at each boundary separately,
    rt_dist is calculated agnostic to choice.
    :param tvec: 1-D array containing times to evaluate pdf
    :param mu: drift rate 2xlen(tvec) array (to incorporate any urgency signal)
    :param bound: default [1 1]
    :param num_images: number of images for method of images, default 7
    :return: probability of correct choice (p_up), and decision time distribution (rt_dist)

    NOTE bound values are flipped to set the starting point as a negative value from 0, within the lower left quadrant.
    margin_width then determines the area above the bound (now at 0) at which the cdf is calculated (in practice we 
    have to calculate the cdf over some area. The extent to which the value of margin_width affects results
    has not been tested yet.) 

    """
    sigma, k = _corr_num_images(num_images)

    survival_prob = np.ones_like(tvec)
    flux1, flux2 = np.zeros_like(tvec), np.zeros_like(tvec)

    s0 = -bound
    b0, bm = -margin_width, 0
    bound0 = np.array([b0, b0])
    bound1 = np.array([b0, bm])  # top boundary of third quadrant
    bound2 = np.array([bm, b0])  # right boundary of third quadrant

    # for under the hood call to mvnun
    low = np.asarray([-np.inf, -np.inf]) # evaluate cdf from -inf to 0 (bound)
    opts = dict(maxpts=None, abseps=1e-5, releps=1e-5)

    # calling the lower-level Fortran for generating the mv normal distribution is MUCH MUCH faster
    # lots of overhead associated with repeated calls of mvn.cdf...
    # downside is that this is a private function, so have to be more careful as it skips a lot of
    # typical checks e.g. on positive definite-ness of cov matrix. It could also change in
    # future Scipy releases without warning...
    use_mvnun = True

    # skip the first sample (t starts at 1)
    for t in range(1, len(tvec)):

        # why do we need this?, because otherwise cov becomes zero?
        if tvec[t] == 0:
            # tvec[t] = np.finfo(np.float64).eps
            tvec[t] = np.min(tvec[tvec > 0])

        mu_t = mu[t, :].T * tvec[t]

        # define frozen mv object
        if not use_mvnun:
            mvn_0 = _mvn_timestep(mean=s0 + mu_t, cov=sigma * tvec[t])
            mvn_0.maxpts = 10000*2

            # total density within boundaries
            cdf_rest = mvn_0.cdf(bound0)

            # density beyond boundary, in one or other direction
            cdf1 = mvn_0.cdf(bound1) - cdf_rest
            cdf2 = mvn_0.cdf(bound2) - cdf_rest

        else:
            # total density within boundaries
            cdf_rest = _mvn.mvnun(low, bound0, s0 + mu_t, sigma * tvec[t], **opts)[0]

            # density beyond boundary, in one or other direction
            cdf1 = _mvn.mvnun(low, bound1, s0 + mu_t, sigma * tvec[t], **opts)[0] - cdf_rest
            cdf2 = _mvn.mvnun(low, bound2, s0 + mu_t, sigma * tvec[t], **opts)[0] - cdf_rest

        # loop over images
        for j in range(1, k*2):
            sj = _sj_rot(j, s0, k)

            if not use_mvnun:
                mvn_j = _mvn_timestep(mean=sj + mu_t, cov=sigma * tvec[t])
                mvn_j.maxpts = 10000*2

                # total density WITHIN boundaries for jth image
                cdf_add = mvn_j.cdf(bound0)

                # density BEYOND boundary in one or other direction, for jth image
                cdf_add1 = mvn_j.cdf(bound1) - cdf_add
                cdf_add2 = mvn_j.cdf(bound2) - cdf_add

            else:
                # total density WITHIN boundaries for jth image
                cdf_add = _mvn.mvnun(low, bound0, sj + mu_t, sigma * tvec[t], **opts)[0]

                # density BEYOND boundary in one or other direction, for jth image
                cdf_add1 = _mvn.mvnun(low, bound1, sj + mu_t, sigma * tvec[t], **opts)[0] - cdf_add
                cdf_add2 = _mvn.mvnun(low, bound2, sj + mu_t, sigma * tvec[t], **opts)[0] - cdf_add

            a_j = _weightj(j, mu[t, :].T, sigma, sj, s0)
            cdf_rest += (a_j * cdf_add)
            cdf1 += (a_j * cdf_add1)
            cdf2 += (a_j * cdf_add2)

        survival_prob[t] = cdf_rest
        flux1[t] = cdf1
        flux2[t] = cdf2

    p_up = np.sum(flux2) / np.sum(flux1 + flux2)

    # NOTE for correct vs error RT distributions, presumably calculate two survival probs from flux1 and flux2 
    rt_dist = np.diff(np.insert(1-survival_prob, 0, 0))

    # winning and losing pdfs?
    # pdf_up = np.diff(flux2)
    # pdf_lo = np.diff(flux1)

    return p_up, rt_dist, flux1, flux2


def _moi_dv(mu: np.ndarray, s: np.ndarray = np.array([1, 1]), num_images: int = 7) -> np.ndarray:

    sigma, k = _corr_num_images(num_images)

    V = np.diag(s) * sigma * np.diag(s)

    dv = np.zeros_like(mu)

    # FIXME default call to rvs may end up resulting in the same result each time because the numpy random number generator
    # has the same seed on each function call. Not sure why...
    # TODO repeated rvs calls (calling at each timepoint) ends up being quite slow, consider Cholesky alternative below
    for t in range(1, mu.shape[0]):
        dv[t, :] = mvn(mu[t, :].T, cov=V).rvs()

    dv = dv.cumsum(axis=0)

    return dv

# maybe useful for faster dv simulation (replacement for rvs calls)
# def chol_sample(mean, cov):
#     return mean + np.linalg.cholesky(cov) @ np.random.standard_normal(mean.size)


def _multiple_logpdfs_vec_input(xs, means, covs):
    """multiple_logpdfs` assuming `xs` has shape (N samples, P features).

    https://gregorygundersen.com/blog/2020/12/12/group-multivariate-normal-pdf/
    
    Thanks to the above link, this provides a much faster way of computing the pdfs across time
    compared to calling mvn pdf at each timepoint. 
    TODO I did some crude checks using np.allclose that it gives the same results, but a unit test would be much better...
    """
    # NumPy broadcasts `eigh`.
    vals, vecs = np.linalg.eigh(covs)

    # Compute the log determinants across the second axis.
    logdets = np.sum(np.log(vals), axis=1)

    # Invert the eigenvalues.
    valsinvs = 1./vals

    # Add a dimension to `valsinvs` so that NumPy broadcasts appropriately.
    Us = vecs * np.sqrt(valsinvs)[:, None]
    devs = xs[:, None, :] - means[None, :, :]

    # Use `einsum` for matrix-vector multiplications across the first dimension.
    devUs = np.einsum('jnk,nki->jni', devs, Us)

    # Compute the Mahalanobis distance by squaring each term and summing.
    mahas = np.sum(np.square(devUs), axis=2)

    # Compute and broadcast scalar normalizers.
    dim = xs.shape[1]
    log2pi = np.log(2 * np.pi)

    out = -0.5 * (dim * log2pi + mahas + logdets[None, :])
    return out.T


def _urgency_scaling(mu: np.ndarray, tvec: np.ndarray, urg=None) -> np.ndarray:
    """Scale mu according to urgency vector."""
    if len(mu) != len(tvec):
        mu = np.tile(mu, (len(tvec), 1))

    if urg is not None:
        if isinstance(urg, (int, float)):
            urg = np.ones(len(tvec)-1) * urg/(len(tvec)-1)
            urg = np.insert(urg, 0, 0)

        assert len(urg) == len(tvec) == len(mu),\
            "If urgency is a vector, it must match lengths of tvec and mu"

        mu = mu + urg.reshape(-1, 1)

    return mu

def log_odds(pdf1: np.ndarray, pdf2: np.ndarray) -> np.ndarray:
    """
    Calculate log posterior odds of correct choice.

    assumes that drift is the first dimension, which gets marginalized over
    :param pdf1: pdf of losing race for correct trials
    :param pdf2: pdf of losing race for incorrect trials
    :return log_odds_correct: heatmap, log posterior odds of correct choice
    """
    # replaces zeros with tiny value to avoid logarithm issues
    min_val = np.finfo(np.float64).tiny
    pdf1 = np.clip(pdf1, a_min=min_val, a_max=None)
    pdf2 = np.clip(pdf2, a_min=min_val, a_max=None)

    odds = np.sum(pdf1, axis=0) / np.sum(pdf2, axis=0)
    odds = np.clip(odds, a_min=1, a_max=None)
    return np.log(odds)


def log_pmap(pdf: np.ndarray, q: int = 30) -> np.ndarray:
    """Set cut-off on log odds map, for better visualization."""
    pdf = np.clip(pdf, a_min=10**(-q), a_max=None)
    return (np.log10(pdf)+q) / q

