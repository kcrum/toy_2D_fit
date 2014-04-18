import physicsPDFs as pdfs
import sys, time, subprocess
import os.path
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.stats as st
import pandas as pd

#
# 1) Evaluate bin fraction vectors (for example see ../gaussPlusConstFit.py)
#     1a?) How low can the population be in highest E, highest T bin (which is 
#          the least populated) before Gaussian assumption breaks down? About 10
#          events? If so, the Poisson CDF at x = 10 for lambda = 20 is ~1%.
#     --> find minimum of numpy array by calling myarr.min()
#     --> create 2-D frac. array by calling np.outer(parabfracs, expfracs)
#
# Loop over: 
# 2) Pull spectra into a histogram
#     2a?) Pull same events into 2-D and 1-D hists for identical experiments
# 3) Minimize least squares for norm. params
#     3a ) At least for 1-D case, have option to plot BF func. on hist.
#     3b?) Use both Neyman and Pearson treatments for errors
# 4) Add best-fit params., errors on params., and chi^2 values to arrays    

# 5) After loop, fit param. dists. with Gaussians and chi^2 curve
#

#########################################################################
# This runs the main fake fit loop. If you'd like to see the output of just one
# fake fit (performing the 1-D and 2-D versions), you can call mainloop from the
# terminal like this: t2d.mainloop(1,1000,100,debug=True)
def mainloop(nexpers, nevents0, nevents1, endpoint0=12.0, endpoint1=8.0, 
             lifetime0=260, lifetime1=170, nEbins=4, nTbins=4, outfilename='',
             minevtsperbin=20, PearsonErrs=True, debug=False):
    nevents = (nevents0, nevents1)
    maxT = max(lifetime0, lifetime1)
    maxE = max(endpoint0, endpoint1)
    # Get energy and time PDFs for our two 'isotopes.'
    pdfsE = [pdfs.ParabolicPDF(endpoint0), pdfs.ParabolicPDF(endpoint1)]
    pdfsT = [pdfs.TruncatedExponentialPDF(lifetime0, maxT),
               pdfs.TruncatedExponentialPDF(lifetime1, maxT)]
    # Get vectors of expected fractional bin content.
    fracsE = [pdfsE[0].binfractionvector(nEbins, (0,maxE)),
              pdfsE[1].binfractionvector(nEbins, (0,maxE))]
    fracsT = [pdfsT[0].binfractionvector(nTbins, (0,maxT)),
                 pdfsT[1].binfractionvector(nTbins, (0,maxT))]
    # Take outer product of vectors to get 2-D array. Rows are energy bins,
    # columns are time bins, meaning all events in the same column occurred in
    # the same deltaT window. 
    fracs2D = [np.outer(fracsE[0],fracsT[0]), np.outer(fracsE[1],fracsT[1])]
    # Make sure no bin has fewer than 'minevtsperbin' events.
    checkminbins(minevtsperbin, nevents, fracsE, fracsT, fracs2D)

    if outfilename and os.path.exists(outfilename):
        print 'Warning: overwriting %s.' % outfilename

    data = pd.DataFrame(outdict(nexpers))
    # Loop over fake experiments
    for i in xrange(nexpers):
        print 'Experiment number %s' % i
        # Create arrays of fake events
        dataE, dataT = throwexperiment(nevents, pdfsE, pdfsT)
        # Bin events
        histE = np.histogram(dataE, bins=nEbins, range=(0,maxE))
        histT = np.histogram(dataT, bins=nTbins, range=(0,maxT))
        # Bin 2-D data and match pattern of fracs2D by having time vary by column
        # and energy by row.
        hist2D = np.histogram2d(dataE, dataT, bins=(nEbins,nTbins),
                                range=((0., maxE), (0., maxT)))        
        # Chi-square min. fit of two 1-D histograms        
        nfit1D, cov1D, chi1D, pval1D = fit1D(histE[0], histT[0], fracsE, fracsT,
                                             nevents, PearsonErrs=PearsonErrs,
                                             debug=debug)
        # Max. likelihood fit of two 1-D histograms
        nfit1Dml, fncmin1D = mlfit1D(histE[0], histT[0], fracsE, fracsT, nevents,
                                     debug=debug)
        # Chi-square min. fit of one 2-d histogram
        nfit2D, cov2D, chi2D, pval2D = fit2D(hist2D[0],fracs2D, nevents,
                                             PearsonErrs=PearsonErrs,
                                             debug=debug)
        # Max. likelihood fit of one 2-D histogram
        nfit2Dml, fncmin2D = mlfit2D(hist2D[0],fracs2D, nevents, debug=debug)
        # Fill outputs into dataframe
        adddata(data, i, nfit1D, cov1D, chi1D, pval1D, nfit2D, cov2D, chi2D,
                pval2D, nfit1Dml, fncmin1D, nfit2Dml, fncmin2D)
        
    if outfilename: data.to_csv(outfilename)
    print 'Main loop finished!'

#########################################################################
# This fills a dataframe built from the output of 'outdict()' with a given fake
# experiment's outputs.
def adddata(df, i, nfit1D, cov1D, chi1D, pval1D, nfit2D, cov2D, chi2D, pval2D,
            nfit1DML, fncmin1D, nfit2DML, fncmin2D): 
    # 1-D chi square outputs
    df['n0_1D'][i], df['n1_1D'][i] = nfit1D[0], nfit1D[1]
    df['var00_1D'][i], df['var01_1D'][i], df['var11_1D'][i] = cov1D[0][0],\
                                                              cov1D[0][1],\
                                                              cov1D[1][1]
    df['chi_1D'][i], df['pval_1D'][i] = chi1D, pval1D
    # 2-D chi square outputs
    df['n0_2D'][i], df['n1_2D'][i] = nfit2D[0], nfit2D[1]
    df['var00_2D'][i], df['var01_2D'][i], df['var11_2D'][i] = cov2D[0][0],\
                                                              cov2D[0][1],\
                                                              cov2D[1][1]
    df['chi_2D'][i], df['pval_2D'][i] = chi2D, pval2D
    # Max. likelihood outputs
    df['n0_1DML'][i], df['n1_1DML'][i] = nfit1DML[0], nfit1DML[1]
    df['n0_2DML'][i], df['n1_2DML'][i] = nfit2DML[0], nfit2DML[1]
    df['fncmin_1DML'][i], df['fncmin_2DML'][i] = fncmin1D, fncmin2D

#########################################################################
# This creates a dict for outputting the data. Each value in the dict is an array
# of length 'nexps.'
def outdict(nexps):
    return {'n0_1D': np.zeros(nexps), 'n1_1D': np.zeros(nexps),
            'var00_1D': np.zeros(nexps), 'var01_1D': np.zeros(nexps),
            'var11_1D': np.zeros(nexps), 'chi_1D': np.zeros(nexps),
            'pval_1D': np.zeros(nexps), 'n0_2D': np.zeros(nexps),
            'n1_2D': np.zeros(nexps), 'var00_2D': np.zeros(nexps),
            'var01_2D': np.zeros(nexps), 'var11_2D': np.zeros(nexps),
            'chi_2D': np.zeros(nexps), 'pval_2D': np.zeros(nexps),
            'n0_1DML': np.zeros(nexps), 'n1_1DML': np.zeros(nexps),
            'fncmin_1DML': np.zeros(nexps), 'n0_2DML': np.zeros(nexps),
            'n1_2DML': np.zeros(nexps), 'fncmin_2DML': np.zeros(nexps)}

#########################################################################
# Find best fit 'isotope' rates for the energy and time variables binned
# separately by maximizing likelihood.
def mlfit1D(binnedE, binnedT, fracsE, fracsT, nevents, debug=False):
    # Concatenate data and prediction vectors
    datavec = np.append(binnedE,binnedT)
    fracvec = [np.append(fracsE[0],fracsT[0]), np.append(fracsE[1],fracsT[1])]
    fnc = lambda p: -np.sum(np.log(st.poisson.pmf(datavec, fracvec[0]*p[0]
                                                  + fracvec[1]*p[1])))
    pfit, fncmin, mingrad, invhess, ncalls, ngradcalls, wflag = \
        sp.optimize.fmin_bfgs(fnc, nevents, full_output=1, disp=True)
    #pfit, fncmin, direc, niter, ncalls, wflag = \
    #    sp.optimize.fmin_powell(fnc, nevents, full_output=True, disp=True)
    if debug:
        print '---------------------- 1-D ML Fit ------------------------------'
        print 'Best fits: %s' % pfit
        print 'Min. fnc. val: %s' % fncmin
        print 'Num. fnc. calls: %s' % ncalls
        print '---------------------------------------------------------------'
    return pfit, fncmin

#########################################################################
# Find best fit 'isotope' rates when time and energy variables of fake data are
# binned together in 2-D histogram
def mlfit2D(binneddata, fracs2D, nevents, debug=False):
    # Concatenate data and prediction vectors
    datavec = binneddata.flatten()
    predfunc = lambda p: p[0]*fracs2D[0].flatten() + p[1]*fracs2D[1].flatten()

    fnc = lambda p: -np.sum(np.log(st.poisson.pmf(datavec, predfunc(p))))
    pfit, fncmin, mingrad, invhess, ncalls, ngradcalls, wflag = \
        sp.optimize.fmin_bfgs(fnc, nevents, full_output=1, disp=True)
    #pfit, fncmin, direc, niter, ncalls, wflag = \
    #    sp.optimize.fmin_powell(fnc, nevents, full_output=True, disp=True)
    if debug:
        print '---------------------- 2-D ML Fit ------------------------------'
        print 'Best fits: %s' % pfit
        print 'Min. fnc. val: %s' % fncmin
        print 'Num. fnc. calls: %s' % ncalls
        print '---------------------------------------------------------------'
    return pfit, fncmin

#########################################################################
# Find best fit 'isotope' rates by minimizing a chi-square (treats Poisson errors
# as Gaussian) when time and energy variables of fake data are binned together in
# 2-D histogram 
def fit2D(binneddata, fracs2D, nevents, PearsonErrs=True, debug=False):
    datavec = binneddata.flatten()
    predfunc = lambda p: p[0]*fracs2D[0].flatten() + p[1]*fracs2D[1].flatten()
    func = lambda : 1
    if PearsonErrs: func = lambda p: (datavec - predfunc(p))/np.sqrt(predfunc(p))
    else: func = lambda p: (datavec - predfunc(p))/np.sqrt(datavec)
    pfit, pcov, infodict, errmsg, success = sp.optimize.leastsq(func, nevents,
                                                                full_output=1)
    chi2 = sum([elem**2 for elem in infodict['fvec']])
    dof = binneddata.size - 2
    pval = st.chisqprob(chi2, dof)
    if debug:
        print '---------------------- 2-D Fit --------------------------------'
        print 'Best fits: %s' % pfit
        print 'Cov. mat.: %s' % pcov
        print 'Chi^2: %s' % chi2
        print 'd.o.f.: %s' % dof
        print 'P-value: %s' % pval
        print '---------------------------------------------------------------'
    return pfit, pcov, chi2, pval

#########################################################################
# Find best fit 'isotope' rates for the energy and time variables binned
# separately by minimizing a chi-square (treats Poisson errors as Gaussian).
def fit1D(binnedE, binnedT, fracsE, fracsT, nevents, PearsonErrs=True,
          debug=False):
    # Concatenate data and prediction vectors
    datavec = np.append(binnedE,binnedT)
    predvec = [np.append(fracsE[0],fracsT[0]), np.append(fracsE[1],fracsT[1])]
    # Define inputs for minimizer
    predfunc = lambda p: p[0]*predvec[0] + p[1]*predvec[1]
    func = lambda : 1
    if PearsonErrs: func = lambda p: (datavec - predfunc(p))/np.sqrt(predfunc(p))
    else: func = lambda p: (datavec - predfunc(p))/np.sqrt(datavec)

    pfit, pcov, infodict, errmsg, success = sp.optimize.leastsq(func, nevents,
                                                                full_output=1)
    chi2 = sum([elem**2 for elem in infodict['fvec']]) 
    #mychi2 = sum([(datavec[i] - predfunc(pfit)[i])**2./predfunc(pfit)[i] for i in xrange(len(datavec))])  ### This just equals 'chi2' calculated above
    dof = datavec.size - 2
    pval = st.chisqprob(chi2, dof)
    if debug:
        print '---------------------- 1-D Fit --------------------------------'
        print 'Best fits: %s' % pfit
        print 'Cov. mat.: %s' % pcov
        print 'Chi^2: %s' % chi2
        print 'd.o.f.: %s' % dof
        print 'P-value: %s' % pval
        print '---------------------------------------------------------------'
    return pfit, pcov, chi2, pval

#########################################################################
# Create arrays of energy and deltaT points drawn from the energy and time PDFs
# of our two 'isotopes'.
def throwexperiment(nevents, pdfsE, pdfsT):
    energyarr = np.array([])
    timearr = np.array([])
    for niso in xrange(len(nevents)):
        energyarr = np.append(energyarr,pdfsE[niso].rvs(size=nevents[niso]))
        timearr = np.append(timearr,pdfsT[niso].rvs(size=nevents[niso]))
    return energyarr, timearr

#########################################################################
# This finds the bins with the fewest expected events in both the 1-D and 2-D 
# cases.
def checkminbins(minevtsperbin, nevents, fracsE, fracsT, twoDfracs):
    min1Dbin = min(np.min(nevents[0]*fracsE[0] + nevents[1]*fracsE[1]),
                   np.min(nevents[0]*fracsT[0] + nevents[1]*fracsT[1]))    
    min2Dbin = np.min(nevents[0]*twoDfracs[0] + nevents[1]*twoDfracs[1]) 

    if min1Dbin <= minevtsperbin:
        print 'Warning: there is a bin in the 1-D fit with %s expected events. Either increase event rate or make binning coarser. Exiting.'
        sys.exit(1)
    if min2Dbin <= minevtsperbin:
        print 'Warning: there is a bin in the 2-D fit with %s expected events. Either increase event rate or make binning coarser. Exiting.'
        sys.exit(1)

#########################################################################
# This allows you to make easy calls to the terminal.
def sh(arg):
    return subprocess.call(arg, shell=True) # This the more modern version of
# os.system(...). The program waits until this finishes before proceding. If you
# don't pass 'shell=True', all whitespace-separated parts of the shell command
# must be passed in like a list. This returns 'returncode'. You could also call:
#     subprocess.Popen(arg, shell=True).wait()
# but this can deadlock (see docs), and it doesn't return 'returncode'.
    
#########################################################################
#########################################################################
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'This requires number of fake experiments, true number of isotope 0 events, and true number of isotope 1 events. Optionally you can pass a filename to save output. Exiting!'
        sys.exit(1)

    nexperiments = int(sys.argv[1])
    nevents0 = int(sys.argv[2])
    nevents1 = int(sys.argv[3])

    outfilename = ''
    if len(sys.argv) == 5: outfilename = sys.argv[4]
    
    starttime = time.clock()
    mainloop(nexperiments, nevents0, nevents1, outfilename=outfilename)
    print 'elapsed time: %s' % (time.clock() - starttime)

    #mainloop(nexpers, nevents0, nevents1, endpoint0=12.0, endpoint1=8.0,
    #         lifetime0=260, lifetime1=170, nEbins=4, nTbins=4, outfilename='',
    #         minevtsperbin=20, PearsonErrs=True, debug=False):
