# -*- coding: utf-8 -*-
"""
Low-level helper functions for 2-D self-motion DDM using method of images.
"""

from typing import Optional

import numpy as np
from scipy.stats import multivariate_normal as mvn
from scipy.stats import norm

USE_MVNUN = False

def _sj_rot(j, s0, k):
    """
    Image rotation formalism.

    :param j: jth image
    :param s0: starting_point, length 2 array
    :param k: 2*k-1 is the number of images - k=4 for 7 images
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


def _corr_num_images(num_images: int) -> tuple[np.ndarray, int]:
    """returns 2-D accumulator correlation given a number of images"""

    k = int(np.ceil(num_images / 2))
    rho = -np.cos(np.pi / k)
    sigma = np.array([[1, rho], [rho, 1]])

    return sigma, k

# %% ----------------------------------------------------------------
# CDF/PDF calculations

# 

def moi_pdf(
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

        pdf_result[t, ...] = pdf_at_timestep(
            tvec[t], mu[t, :], sigma, xy_mesh, k, s0)

    return pdf_result


def moi_pdf_vec(
    xmesh: np.ndarray, 
    ymesh: np.ndarray, 
    tvec: np.ndarray,
    mu: np.ndarray,
    bound=np.array([1, 1]),
    num_images: int = 7
    ) -> np.ndarray:
    
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
            # pdf_result[t, ...] += a_j * mvn(mean=sj + mu[t, :], cov=sigma*tvec[t]).pdf(xy_mesh)
            
            aj_all[t] = a_j

        # use vectorized implementation for speed # TODO unit tests to verify correctness
        pdf_result += (aj_all * np.exp(_multiple_logpdfs_vec_input(new_mesh, sj + mu_t, covs)))

    return pdf_result


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


def pdf_at_timestep(
    t, 
    mu: np.ndarray, 
    sigma: np.ndarray, 
    xy_mesh: np.ndarray, 
    k: int, 
    s0: np.ndarray
    ):
    """Calculate pdf at a single timepoint."""
    
    pdf = mvn(mean=s0 + mu*t, cov=sigma*t).pdf(xy_mesh)

    # j-values start at 1, go to k*2-1
    for j in range(1, k*2):
        sj = _sj_rot(j, s0, k)
        a_j = _weightj(j, mu.T, sigma, sj, s0)
        pdf += a_j * mvn(mean=sj + mu*t, cov=sigma*t).pdf(xy_mesh)

    return pdf


def moi_cdf(
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

    NOTE bound values are flipped to set the starting point as a negative value from 0, 
    within the lower left quadrant, with zero as the bound.
    margin_width then determines the area above the bound (now at 0) where the cdf is calculated (because 
    in practice we have to calculate the cdf over some non-zero area).
    The extent to which the value of margin_width affects results has not been tested yet.
    0.025 seems to work ok for now.

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
    if USE_MVNUN:
        low = np.asarray([-np.inf, -np.inf]) # evaluate cdf from -inf to 0 (bound)
        opts = dict(maxpts=None, abseps=1e-5, releps=1e-5)

    # calling the lower-level Fortran for generating the mv normal distribution is MUCH MUCH faster
    # lots of overhead associated with repeated calls of mvn.cdf...
    # downside is that this is a private function, so have to be more careful as it skips a lot of
    # typical checks e.g. on positive definite-ness of cov matrix. It could also change in 
    # future Scipy releases without warning...
    # https://stackoverflow.com/questions/76524943/how-to-compute-faster-scipy-stats-multivariate-normal-cdf-a-large-number-of-ti

    # skip the first sample (t starts at 1)
    for t in range(1, len(tvec)):

        # why do we need this?, because otherwise cov becomes zero?
        if tvec[t] == 0:
            # tvec[t] = np.finfo(np.float64).eps
            tvec[t] = np.min(tvec[tvec > 0])

        mu_t = mu[t, :].T * tvec[t]

        # define frozen mv object
        if not USE_MVNUN:
            mvn_0 = mvn(mean=s0 + mu_t, cov=sigma * tvec[t])
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

            if not USE_MVNUN:
                mvn_j = mvn(mean=sj + mu_t, cov=sigma * tvec[t])
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


def moi_cdf_vec(
    tvec: np.ndarray,
    mu: np.ndarray,
    bound=np.array([1, 1]),
    margin_width: float = 0.025,
    num_images: int = 7,
    bvn_n: int = 64,
):
    """
    Vectorized moi CDF over times (T) and images (J).
    Uses vectorized bivariate normal CDF implementation
    
    Returns same outputs as `moi_cdf`: p_up, rt_dist, flux1, flux2
    """
    sigma, k = _corr_num_images(num_images)
    s0 = -bound
    b0, bm = -margin_width, 0
    bound0 = np.array([b0, b0])
    bound1 = np.array([b0, bm])
    bound2 = np.array([bm, b0])

    # safe copy so we don't mutate input
    tvec_safe = np.array(tvec, dtype=float, copy=True)
    tvec_safe[tvec_safe == 0] = np.min(tvec_safe[tvec_safe > 0])

    T = len(tvec_safe)
    # build all image starting positions (J,2)
    sj_list = [s0] + [_sj_rot(j, s0, k) for j in range(1, k*2)]
    sj_all = np.vstack(sj_list)             # (J,2)
    J = sj_all.shape[0]

    # means: (T, J, 2)
    mu_t = mu * tvec_safe[:, None]         # (T,2)
    means = mu_t[:, None, :] + sj_all[None, :, :]  # (T,J,2)

    # weights: (T, J)
    inv_sigma = np.linalg.inv(sigma)
    delta = (sj_all - s0)                   # (J,2)
    exponent = (mu @ inv_sigma) @ delta.T   # (T,J)
    signs = np.array([1] + [(-1)**j for j in range(1, k*2)])
    weights = signs[None, :] * np.exp(exponent)

    sd = np.sqrt(tvec_safe)[:, None]        # (T,1)
    rho_param = sigma[0, 1]

    # standardized limits shape (T,J)
    h0 = (bound0[0] - means[..., 0]) / sd
    k0 = (bound0[1] - means[..., 1]) / sd
    h1 = (bound1[0] - means[..., 0]) / sd
    k1 = (bound1[1] - means[..., 1]) / sd
    h2 = (bound2[0] - means[..., 0]) / sd
    k2 = (bound2[1] - means[..., 1]) / sd

    # compute cdfs
    cdf0 = _bvn_cdf(h0, k0, rho_param, n=bvn_n)
    cdf1_arr = _bvn_cdf(h1, k1, rho_param, n=bvn_n) - cdf0
    cdf2_arr = _bvn_cdf(h2, k2, rho_param, n=bvn_n) - cdf0
    
    # reduce over images
    cdf_rest = np.sum(weights * cdf0, axis=1)
    cdf1 = np.sum(weights * cdf1_arr, axis=1)
    cdf2 = np.sum(weights * cdf2_arr, axis=1)

    survival_prob = np.ones_like(tvec_safe)
    survival_prob[1:] = cdf_rest[1:]
    flux1 = np.zeros_like(tvec_safe); flux2 = np.zeros_like(tvec_safe)
    flux1[1:] = cdf1[1:]; flux2[1:] = cdf2[1:]

    p_up = np.sum(flux2) / np.sum(flux1 + flux2)
    rt_dist = np.diff(np.insert(1 - survival_prob, 0, 0))

    return p_up, rt_dist, flux1, flux2


def _bvn_cdf(h, k, rho, n=64):
    """Vectorized bivariate normal CDF Phi_2(h,k; rho).

    The naive implementation of the bivariate normal CDF in SciPy
    (scipy.stats.multivariate_normal.cdf) is quite slow because it requires
    a loop over each timepoint and image, to account for the changing mean and covariance
    at each timepoint. This function provides a vectorized alternative using
    Gauss-Legendre quadrature

    Uses a transformation to map the integral from (-inf, h) to (-1,1) and
    applies Gauss-Legendre quadrature. Accepts scalars or 1-D arrays for h,k,rho
    and returns an array of the same shape.
    """
    from numpy.polynomial.legendre import leggauss
    h = np.asarray(h)
    k = np.asarray(k)
    rho = np.asarray(rho)

    # broadcast to common shape
    h, k, rho = np.broadcast_arrays(h, k, rho)

    # get nodes and weights on [-1,1]
    u, w = leggauss(n)
    # reshape weights for broadcasting over the mesh (nodes, *h.shape)
    w = w[:, None, None]

    # transform nodes to x in (-inf, h) via x = h - (1+u)/(1-u)
    # dx/du = -2/(1-u)^2, take absolute value for jacobian (this allows for integration by substitution)
    denom = (1 - u)
    # expand node-dependent terms to shape (n,1,1) so they broadcast correctly over h's shape
    x = h[None, ...] - ((1.0 + u) / denom)[:, None, None]
    jac = (2.0 / (denom ** 2))[:, None, None]

    # compute phi(x) and inner normal CDF Phi((k - rho*x)/sqrt(1-rho^2))
    # phi(x) = exp(-x^2/2)/sqrt(2*pi)
    phi_x = np.exp(-0.5 * x**2) / np.sqrt(2.0 * np.pi)

    denom_r = np.sqrt(1.0 - rho**2)
    inner = (k[None, ...] - rho[None, ...] * x) / denom_r[None, ...]

    # use scipy's univariate normal CDF (vectorized)
    phi_inner = norm.cdf(inner)

    integrand = phi_x * phi_inner * jac

    # integrate using weights (sum over nodes)
    integral = np.sum(w * integrand, axis=0)

    return integral.reshape(h.shape)


def sample_dv(
    mu: np.ndarray,
    s: np.ndarray = np.array([1, 1]),
    num_images: int = 7,
    seed: Optional[int] = None
    ) -> np.ndarray:

    sigma, k = _corr_num_images(num_images)
    V = np.diag(s) * sigma * np.diag(s)
    
    dv = np.zeros_like(mu)
    T = mu.shape[0]

    # Vectorized sampling: draw all increments at once using a standard normal
    # and transform with Cholesky decomposition of V. Cheap trick to simulate 
    # diffusion according to covariance matrix V, without repeated calls to
    # scipy.stats.multivariate_normal.rvs
    
    if T > 1:

        # for t in range(1, T):
        #     dv[t, :] = mvn(mu[t, :].T, cov=V).rvs()
        # dv = dv.cumsum(axis=0)
        
        rng = np.random.default_rng(seed)
        L = np.linalg.cholesky(V)
        Z = rng.standard_normal(size=(T - 1, mu.shape[1]))
        
        increments = mu[1:] + (Z @ L.T)
        dv[1:] = np.cumsum(increments, axis=0)

    return dv




# --- Small helper for quick comparisons ------------------------------------

def compare_moi_cdfs(tvec, mu, **kwargs):
    """
    Quick comparison helper between original `moi_cdf` and `moi_cdf_vec`.
    Returns dict with results and absolute differences for flux arrays.
    """
    out_orig = moi_cdf(tvec.copy(), mu.copy(), **kwargs)
    out_vec = moi_cdf_vec(tvec.copy(), mu.copy(), **kwargs)
    res = {
        'orig': out_orig,
        'vec': out_vec,
        'flux1_err': np.max(np.abs(out_orig[2] - out_vec[2])),
        'flux2_err': np.max(np.abs(out_orig[3] - out_vec[3])),
    }
    return res