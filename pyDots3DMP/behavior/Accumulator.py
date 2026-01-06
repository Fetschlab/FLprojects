import logging
from typing import Optional, Union, Sequence

import numpy as np

# always import from base matplotlib first to avoid backend issues
import matplotlib as mpl
import matplotlib.pyplot as plt

from .moi import moi_cdf, moi_cdf_vec, moi_pdf, moi_pdf_vec, sample_dv

logger = logging.getLogger(__name__)
    
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

    # Use __slots__ to reduce per-instance memory usage and speed attribute access.
    __slots__ = (
        'tvec', 'grid_vec', '_bound', 'drift_rates', 'num_images', 'wager_theta',
        '_is_fitted', 'drift_labels', 'p_corr_', 'rt_dist_', 'pdf3D_', 'up_lose_pdf_',
        'lo_lose_pdf_', 'log_odds_'
    )

    def __init__(
        self,
        tvec: np.ndarray,
        grid_vec: np.ndarray, 
        drift_rates: Optional[list] = None, 
        bound: Optional[Union[float, Sequence, np.ndarray]] = 1.0,
        num_images: Optional[int]=7,
        wager_theta: Optional[float]=1.0
        ):
        self.tvec = tvec
        self.grid_vec = grid_vec
        self.bound = bound
        self.drift_rates = drift_rates if drift_rates is not None else []
        self.num_images = num_images
        self.wager_theta = wager_theta
        self._is_fitted = False

        if self.drift_rates:
            self.drift_labels = self.drift_rates.copy()


    @property
    def is_fitted(self) -> bool:    
        return self._is_fitted

    
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


    def apply_drifts(
        self,
        drifts,
        labels: Optional[list] = None, 
        sensitivity: float = 1.0,
        urgency: Optional[Union[np.ndarray,float]] = None
        ):
        """
        Set accumulator drift rates. Optionally add label for each drift.
        This also adds a mirrored drift rate for the anti-correlated accumulator, and 
        updates drift rates based on sensitivity and urgency parameters.
        """

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
            

    def cdf(self, use_vectorized: bool = True):
        """calculate cdf at boundaries for each drift rate, returns
        probability of correct choice and RT distribution (no NDT)"""
        
        p_corr = np.zeros(len(self.drift_rates))
        rt_dist = np.zeros((len(self.drift_rates), len(self.tvec)))

        cdf_fcn = moi_cdf_vec if use_vectorized else moi_cdf
        
        for d, drift in enumerate(self.drift_rates):
            p_corr[d], rt_dist[d, :], flux1, flux2 = cdf_fcn(
                self.tvec, drift, self.bound, 0.025, self.num_images
                )
            
        self.p_corr_ = p_corr
        self.rt_dist_ = rt_dist


    def pdf(self, use_vectorized=True, full_pdf=False):

        # TODO allow flexible specification of grid_vec, to use mgrid
        if full_pdf:
            xmesh, ymesh = np.meshgrid(self.grid_vec, self.grid_vec)
        else:
            xmesh1, ymesh1 = np.meshgrid(self.grid_vec, self.grid_vec[-1])
            xmesh2, ymesh2 = np.meshgrid(self.grid_vec[-1], self.grid_vec)

        pdfs, marg_up, marg_lo = [], [], []

        for drift in self.drift_rates:

            if full_pdf:
                # return the full pdf in "3D" (x-y space, over time)
                pdf_3d = moi_pdf(xmesh, ymesh, self.tvec, drift, self.bound, self.num_images)
                pdfs.append(pdf_3d)

                pdf_up = pdf_3d[:, :, -1]  # right bound
                pdf_lo = pdf_3d[:, -1, :]  # top bound

            else:
                # sufficient to calculate pdf just at the boundaries, not the full third quadrant
                # and use vectorized version for speed
                pdf_fcn = moi_pdf_vec if use_vectorized else moi_pdf
                pdf_lo = pdf_fcn(xmesh1, ymesh1, self.tvec, drift, self.bound, self.num_images)
                pdf_up = pdf_fcn(xmesh2, ymesh2, self.tvec, drift, self.bound, self.num_images)

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
        return sample_dv(
            mu=drift*self.tvec.reshape(-1, 1),
            s=sigma,
            num_images=self.num_images
            )


    def compute_distrs(self, return_pdf=False):
        """Calculate cdf and pdf for accumulator object. Returns self for chaining commands"""

        self.cdf()
        if return_pdf:
            self.pdf()
        self._is_fitted = True
        
        return self


    def plot(
        self,
        d_ind: int = -1,
        save_path: Optional[str] = None,
        ):
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
        if not hasattr(self, 'p_corr_'):
            raise ValueError('Accumulator distributions have not yet been calculated')

        fig_cdf, axc = plt.subplots(2, 1, figsize=(4, 5))
        axc[0].plot(self.drift_labels, self.p_corr_, marker='o')
        axc[0].set_ylim([0, 1])
        axc[0].set_xlabel('drift')
        axc[0].set_xticks(self.drift_labels)
        axc[0].set_ylabel('prob. correct choice')
        axc[1].set_title('Accumulator CDF/PDF Results')

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
            fig_pdf.tight_layout()

            if has_log_odds:
                vmin, vmax = 0, 3
                contour = axp[2].contourf(self.tvec, self.grid_vec,
                                            self.log_odds_.T, vmin=vmin, vmax=vmax,
                                            levels=100)
                axp[2].set_title("Log Odds of Correct Choice given Losing Accumulator")
                cbar = fig_pdf.colorbar(contour, ax=axp[2])

        if save_path:
            fig_cdf.savefig(f"{save_path}/cdf.png")
            if fig_pdf is not None:
                fig_pdf.savefig(f"{save_path}/pdf.png")

        return fig_cdf, fig_pdf


    def plot_3d(
        self,
        drift_ind=-1,
        save_path: str = 'pdf_animation',
        filetype: str = '.mp4'
        ):
        """Create and save animation of the full PDF over time.

        Notes
        -----
        - This uses `imshow` and blitting to update frames efficiently.
        - Requires `pdf3D_` to be computed with `full_pdf=True` (shape: [n_drifts, n_time, nx, ny]).
        """
        if not hasattr(self, 'pdf3D_'):
            raise ValueError('Full PDF not available: compute pdf(full_pdf=True) first')

        # Precompute log-scaled frames once (shape: n_time, nx, ny)
        frames = log_pmap(self.pdf3D_[drift_ind])
        n_frames = frames.shape[0]

        # determine display range once for consistent color mapping
        vmin, vmax = float(frames.min()), float(frames.max())

        fig, ax = plt.subplots()
        extent = [self.grid_vec[0], self.grid_vec[-1], self.grid_vec[0], self.grid_vec[-1]]

        # draw first frame with imshow (much faster than contourf)
        im = ax.imshow(frames[0].T, origin='lower', extent=extent,
                       aspect='auto', vmin=vmin, vmax=vmax, cmap='viridis', interpolation='nearest')
        drift_str = np.array2string(self.drift_rates[drift_ind][0], formatter={'float_kind': '{:.3e}'.format})
        title = ax.set_title(f"Frame 1 - {self.tvec[0]:.2f}s, drift = {drift_str}")
        fig.colorbar(im, ax=ax)

        def animate(i):
            im.set_data(frames[i].T)
            drift_str = np.array2string(self.drift_rates[drift_ind][i], formatter={'float_kind': '{:.3e}'.format})
            title.set_text(f"Frame {i + 1} - {self.tvec[i]:.2f}, drift = {drift_str}")
            return im, title

        anim = mpl.animation.FuncAnimation(fig, animate, frames=n_frames, blit=True)
        writer = mpl.animation.PillowWriter(fps=10)
        anim.save(f'{save_path}_{self.drift_labels[drift_ind]}.gif', writer=writer)
        plt.close(fig)
        
        return anim


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

