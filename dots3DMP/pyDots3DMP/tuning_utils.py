#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 00:08:10 2023

@author: stevenjerjian
"""

# compute tuning curves and signal correlations
# tuning curve can be a property of neuron class eventually, a an aside

import pdb
import numpy as np
import matplotlib.pyplot as plt

# from scipy.stats import vonmises
from scipy.optimize import curve_fit
from scipy.special import i0

# %%

# fig, ax = plt.subplots(1, 1)
# plt.plot(thesehdgs, temp)
# plt.plot(xhdgs,ypred)

# TODO decorate these to allow variable function inputs


def tuning_vonMises(x, a, k, theta, b):
    return a * np.exp(k * np.cos(np.deg2rad(x - theta))) / \
        (2 * np.pi * i0(k)) + b
    # return a * vonmises.pdf(x, k, loc=theta, scale=1) + b


# could just merge this into the fit_predict function
# def fit_vonMises(y, x):
#     p0 = [np.max(y) - np.min(y), 1.5, x[np.argmax(y)], np.min(y)]
#     popt, pcov = curve_fit(tuning_vonMises, x, y, p0=p0)
#     perr = np.sqrt(np.diag(pcov))  # 1SD error on parameters
#     return popt, perr


def plot_vonMises_fit(func):
    def plot_wrapper(y, x, x_pred):
        y_pred = func(y, x, x_pred)
        plt.plot(x, y, label='y')
        plt.plot(x_pred, y_pred, label='y_pred')
        plt.legend()
        plt.show()
        return y_pred
    return plot_wrapper


# @plot_vonMises_fit
def fit_predict_vonMises(y, x, x_pred, k=1.5):
    # popt, perr = fit_vonMises(y, x)

    if y.any():
        p0 = np.array([np.max(y) - np.min(y), k, x[np.argmax(y)], np.min(y)])
    else:
        p0 = np.array([1, k, 0, 0])

    try:
        popt, pcov = curve_fit(tuning_vonMises, x, y, p0=p0)
    except RuntimeError as re:
        print(re)
        popt = np.full(p0.shape, np.nan)
        perr = np.full(p0.shape, np.nan)
        y_pred = np.full(x_pred.shape, np.nan)
    except ValueError as ve:
        print(ve)
        popt = np.full(p0.shape, np.nan)
        perr = np.full(p0.shape, np.nan)
        y_pred = np.full(x_pred.shape, np.nan)
    except TypeError:
        pdb.set_trace()
    else:
        perr = np.sqrt(np.diag(pcov))  # 1SD error on parameters
        y_pred = tuning_vonMises(x_pred, *popt)

    return y_pred, popt, p0, perr


def tuning_basic(y, x):

    yR = y[x > 0]
    yL = y[x < 0]

    pref_hdg = x[np.argmax(np.abs(y - y.mean()))]
    pref_dir = int(yR.mean() > yL.mean()) + 1  # L, R is 1, 2 by convention
    pref_amp = (yR.mean() - yL.mean()) / np.ravel([yR, yL])

    return pref_hdg, pref_dir, pref_amp


# def tuning_vonMises_err(x, FR, a, k, theta, b):
#    return np.sum((tuning_vonMises(x, a, k, theta, b) - FR) ** 2)


