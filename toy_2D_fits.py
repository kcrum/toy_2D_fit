import physicsPDFs as pdfs
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.stats as st
        
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
def mainloop(nexpers, nevents0, nevents1, endpoint0 = 12.0, endpoint1 = 8.0, 
             lifetime0 = 260, lifetime1 = 170, nEbins = 4, nTbins = 4, 
             minevtsperbin = 20, PearsonErrs = True):
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

    # Loop over fake experiments
    for i in xrange(nexpers):
        # Create arrays of fake events
        dataE, dataT = throwexperiment(nevents, pdfsE, pdfsT)
        # Bin events
        histE = np.histogram(dataE, bins=nEbins, range=(0,maxE))
        histT = np.histogram(dataT, bins=nTbins, range=(0,maxT))
        # Match fracs2D and put have time vary by column and energy by row.
        hist2D = np.histogram2d(dataE, dataT, bins=(nEbins,nTbins),
                                range=((0., maxE), (0., maxT)))
        # Fit two 1-D histograms        
        nfit, cov, chi2 = fit1D(histE[0], histT[0], fracsE, fracsT, nevents,
                                PearsonErrs=PearsonErrs)
        print '---------------------------------------------------------------'
        print 'Best fits: %s' % nfit
        print 'Cov. mat.: %s' % cov
        print 'Chi^2: %s' % chi2
        print '---------------------------------------------------------------'
        
        #print histE[0]
        #print histT[0]
        # Fit one 2-d histogram
        #print hist2D[0]

    print fracs2D[0]*nevents[0] + fracs2D[1]*nevents[1]

#########################################################################
# Find best fit 'isotope' rates for the energy and time variables binned
# separately.
def fit1D(binnedE, binnedT, fracsE, fracsT, nevents, PearsonErrs=True):
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
    return pfit, pcov, chi2

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
#########################################################################
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'This expects an endpoint for the parabolic energy spectrum, a maximum time for the exponential deltaT distribution, and a lifetime for the deltaT distribution (all > 0) . Exiting.'
        sys.exit(1)

    endpoint = float(sys.argv[1])
    maxT = float(sys.argv[2])
    lifetime = float(sys.argv[3])

    if endpoint < 0 or maxT < 0 or lifetime < 0:
        print 'All three parameters must be positive numbers. Exiting.'
        sys.exit(1)
    
    EnergyPDF = pdfs.ParabolicPDF(endpoint)
    DeltaTPDF = pdfs.TruncatedExponentialPDF(maxT, lifetime)

    EnergyPDF.sampleplot2()
    DeltaTPDF.sampleplot()
