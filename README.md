# OlinPMF 1.0

## Olin Positive Matrix Factorization for PMF analysis on mass spectrometer data



Matlab-based PMF package that is optimized for **large datasets**, is capable of **real-time plotting of the factors** during PMF computation, has **two-way interconnectivity to tofTools**, and includes different PMF-types (**UMR-PMF, HR-PMF, binPMF, PSPMF, and RR-PMF**).



### Compatibility

Works where Matlab works (Windows, Linux, Mac). Needs a Matlab license.

Optional interconnectivity to tofTools requires tofTools (from an external source).



### Theory

It is based on equation X=TS\*MS+R, where X is a data matrix, TS is time series matrix of the resolved factors, MS is the mass spectra of the resolved factors, and R is a residual matrix. Briefly, it minimizes R with its elements normalized with the elements in S, which is an uncertainty matrix.



The minimization can be performed using different minimization algorithm, with the default being ALS5MULT1 algorithm, where alternating least squares algorithm is run with 5 different random seeds, and subsequently, a multiplicative solver, with the initial guess being the best of the 5 ALS-results, is used to output the final result.



Matrix operations employ Matlab's optimized matrix computation methods, such as internal QR-decomposition. Optionally, OlinPMF supports GPU-computation for the matrix operations, which can boost the computation in cases with large datasets and/or a powerful GPU.



### PMF-types

#### UMR-PMF

Unit-mass resolution (UMR)-PMF is the most common and simple PMF-method, where the data matrix consists of the time series of mass spectrometer data binned into unit masses.

#### HR-PMF

High resolution (HR)-PMF is also a commonly used PMF-method. In HR-PMF, mass spectral peaks are first assigned with chemical compositions, which are subsequently used to fit the peaks to provide time series of the assigned ions. The resulting time series matrix is then analyzed with HR-PMF.

#### binPMF

Binned PMF (binPMF) first uses binning of the mass spectra (for example, using a bin width of 0.01 Th, instead of 1 Th as in the case of UMR-PMF), then runs the PMF, and finally the MS-output from binPMF is used to assign the chemical compositions.

#### PSPMF

Pan sharpening PMF (PSPMF) is a method that computationally accelerates the binPMF computation, by first solving PMF as a UMR-PMF and then solving binPMF using a subsampled data matrix, i.e., averaged to longer time resolution. The outputs of these two PMF-solutions are combined and used as an initial guess in the final binPMF run. This method improves the computation time by providing better than random initial guesses for the final binPMF run (which could otherwise be computationally expensive due to the data matrix becoming much larger if binned to very narrow bins).

#### RR-PMF

Special version of binPMF that produces raw data-like PMF factors.



### Installation

Copy all these files into a folder.



### Usage

1. Edit runpmf.m and modify the parameters described therein.

   * Initially choose repetitions=1 and pRange=2 to test your data matrix and to initialize required variables for error matrix construction.
2. Run a quick PMF by running the runpmf.m file from Matlab's editor or by running a command 'runpmf' on Matlab's console.
3. Run 'errorParameterCalculation' and manually fit a function displayed in the console to obtain errorParameters.

   * Modify the default errorParameters in runpmf.m to correspond to the ones you just fit. Note that they are instrument-specific and depend on the time averaging, and that they need to be separately derived for UMR-PMF, binPMF, and HR-PMF.
4. Run 'runpmf' with chosen parameters in runpmf.m file.

   * The outputs of PMF (MS, TS, Q/Qexp, ...) are plotted in real-time.
   * The plots update every time OlinPMF finds a better solution.
   * The plots can be zoomed already during the PMF computation.
   * Each number of factors (p) produces its own figure (with the figure number denoting p).
   * Figure 1, plotting Q/Qexp and unexplained variation vs. p, can be closed any time to abort the PMF run. The best solution for each p found so far will be stored even if aborted.
5. PMF results for each p are stored as variables named like UMRPMFRes5Factors, binPMFRes..., HRPMFRes...
6. The plots can be redrawn by first loading a stored result.

   * For example, 'loadpmfresults(5,'UMR')' loads a UMR-type PMF result that was calculated using 5 factors.
   * 'plotpmfresults' plots the loaded PMF result again.



### Interconnectivity with tofTools

#### Importing the data matrix X from a tofTools preprocessed file

After processing mass calibration in tofTools, you can load the saved paramEval.h5 file in OlinPMF by giving the file name as the 'dataSource' parameter in runpmf.m, instead of giving the variable names as a struct.

In case of HR-PMF, you also need to run the peak fitting in tofTools before importing the data in OlinPMF.

#### Exporting the MS from a binPMF solution to a tofTools file

After finding a suitable binPMF solution in OlinPMF, the resulting MS\_final matrix can be written into a tofTools paramEval.h5 file by running 'binPMF\_to\_tofTools' command (you can decide the name of the new tofTools file inside the binPMF\_to\_tofTools.m file).

Once the new tofTools paramEval.h5 file has been written, you can open it using tofTools 'tof\_gui\_eval\_params' function to post-process the binPMF factors' mass spectra (i.e., peak identification etc.).

#### Exporting the MS from an RR-PMF solution to a tofTools file

After finding a suitable UMR-PMF solution in OlinPMF, the resulting TS\_final matrix contains the time series of the UMR-PMF factors.

Those time series can be used to construct the RR-PMF solution using the Half-ALS algorithm. 

Edit the parameters in RRPMF\_to\_tofTools.m file so that they correspond to the RR-PMF settings you want to use.

Then run the 'RRPMF\_to\_tofTools' command and open it using tofTools 'tof\_gui\_eval\_params' function to post-process the RR-PMF factors' mass spectra (i.e., peak identification etc.).



