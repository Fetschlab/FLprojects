README.TXT

-------------------------------------------------------------------
FitDDM_choiceRTconf — Miguel Vivar-Lazo, October 2024
-------------------------------------------------------------------

MATLAB routines (with C source code) for fitting serial and parallel accumulator model variants
to behavioral data consisting of binary choice, reaction time (RT), and binary confidence (i.e. wager).

Algorithms are based on Kiani & Shadlen 2009 and Shan et al. 2019, respectively; see corresponding wrapper scripts for more info,
and the Methods section of the submitted ms for equations and details on the fitting procedures.

Currently shared (under the GNU General Public License) for peer review only, for manuscript submitted to Nature Neuroscience:
M. Vivar-Lazo & C. R. Fetsch, Neural basis of concurrent deliberation toward a choice and degree of confidence (2024)

***NOTE, current version throws an error when running the Parallel model, although the Serial model (default) works fine.
***This will be fixed soon (will unsubmit and resubmit the capsule)


-------------------------------------------------------------------
 Contents
-------------------------------------------------------------------

The package contains the following Matlab scripts:

-- Serial_Wrapper.m		: wrapper script (the thing to run) for serial model
-- Serial_Fitting.m		: fitting function for serial model; uses FMINSEARCHBND and FP4*
-- Parallel_Wrapper.m		: wrapper script (the thing to run) for parallel model
-- Parallel_Fitting.m		: fitting function for parallel model; uses FMINSEARCHBND and FP_Drugowitsch*
*the 'run' file in Code Ocean will compile the required .mex files from source

and three subfolders:

-- FMINSEARCHBND    : bounded version of fminsearch for optimization, from Mathworks File Exchange (John D'Errico, 2024)
-- FP4	            : algorithm for 1D DDM (see Kiani et al. 2009 or Fetsch et al. 2014)
-- FP_Drugowitsch   : algorithm for 2D DDM (see Shan et al., 2019)
-- utils            : some helper functions


-------------------------------------------------------------------
 Instructions
-------------------------------------------------------------------

Note: Code has been tested on 64bit Linux and MacOS, but if running locally you should always compile your own mex files.

1. Choose one of the two models to run and identify the corresponding wrapper script. If in Code Ocean, you need to edit the 'run' file accordingly.

2. Decide if you want to converge precisely on the best-fitting parameters or just want a quick-and-dirty test run. If the former, find the "options = optimset(" code in the wrapper script and comment/uncomment the lines as instructed. This could take anywhere from tens of minutes to several hours! The latter is the default, currently set to maxFunEvals = 100 which takes about 15 min; change it to 10 to go very fast but the fit is bad.

3. Run it (the wrapper, or click 'Reproducible Run').

4. View the figure file in the /results folder. Model parameters and error/BIC values (in fact the entire workspace) will be saved there too. The .fig file is not viewable from within Code Ocean but can be downloaded and opened in MATLAB.