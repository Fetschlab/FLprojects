Steven Jerjian
2022-09-22

README on dots3DMP behavior codes.
Work in progress, just including what I consider the "main" codes for now.

The idea is to encapsulate as much of the modular analysis and figure plotting in functions as possible, and then use scripts to call these analyses on a given dataset, to cleanly present an overall workflow or analysis for some presentation purpose (paper, talk, poster...). 
Code list below is therefore divided into scripts and functions (processing or plotting)

SCRIPTS found in /offlineTools/dots3DMP/scripts

These give examples of how to load data and run various processing and plotting functions specified below. 
1) dots3DMP_behavioralAnalysis
	- example analyses/basic data viz
	- shows how dots3DMP_parseData_multiConf is used
2) dots3DMP_LucioTrainingAnalysis
3) dots3DMP_performanceByDay
	- work in progress/part-complete. loops over sessions to plot some of the summary behavior, might be useful for assessing performance within/across days. In practice I've found that within day is the lowest level of resolution worth looking at - within block isn't worth it, unless you subselect for blocks with at least, say, 15 trials of each condition, or 500+ trials total.

4) dots3DMP_cleanMonkeyData
	- run on data struct produced by PLDAPS_preprocessing_rigBdots3DMP to clean it up a bit. Removes some fields, all brfix trials, sessions or days the user wants to exclude, trials with no RT or PDW, blocks shorter than N trials (50), blocks with non-specified headings. Standardizes coherence values across sessions.


processing FUNCTIONS in parent folder - /offlineTools/dots3DMP

1a) dots3DMP_parseData             
	- calculate means (e.g. P(right), meanRT, P(HighBet) for each condition
	- fits logistic

1b) dots3DMP_parseData_byConf
	- same as above, but separately for high and low bet

1c) dots3DMP_parseData_multiConf
	- more general version of _byConf - can split trials by high/low/1-targ
	- OR, provided any arbitrary grouping of trials within stimulations conditions 		(e.g. split by offered high bet reward) 
	- this may have some unforeseen bugs but eventual goal is to have one function that allows user to group trials however they want, and so subsume (2).

2a) dots3DMP_fit_cgauss
	- instead of logistic fit, fit cumulative gaussian to choice, inverted gaussian to confidence, and gaussian to RT. inverted Gaussian for confidence in particular seems to have some issues under certain circumstances, particularly if trying to fit on individual session or when little "U" shape is present in PDW

2b) dots3DMP_fit_cgauss_byConf
	- equivalent to parseData_byConf, but for fitting gaussians to high and low bet choice and RT.

2c) dots3DMP_cgauss_bootstrap_func
	- function implementation of dots3DMP_cgauss_bootstrap - fit gaussian to bootstraps of data, for calculating an 'error'. also calculates weights for each iteration


3) dots3DMP_cueWeights
	- uses gfit to calculate weighing of cues at each coherence level

4) dots3DMP_RTquantiles
	- Kiani et al 2014 style, plotting confidence (or accuracy) vs quantiles of RT, separately for each modality/coherence, and each heading


PLOTTING functions - /offlineTools/dots3DMP/plotting

1) dots3DMP_plots, dots3DMP_plots_splitConf, and dots3DMP_plots_multiConf take outputs from 1a-1c functions above

2) dots3DMP_plots_cgauss_byCoh, _byModDelta, _byConf take output from 2a/2b functions




