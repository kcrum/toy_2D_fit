import physicsPDFs as pdfs
import sys
import numpy as np
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

def mainloop(nevents0, nevents1, endpoint0 = 12.0, endpoint1 = 8.0, 
             lifetime0 = 260, lifetime1 = 170, nEbins = 4, nTbins = 4, 
             minevtsperbin = 20):
    maxT = max(lifetime0, lifetime1)
    maxE = max(endpoint0, endpoint1)
    energypdfs = [pdfs.ParabolicPDF(endpoint0), pdfs.ParabolicPDF(endpoint1)]
    timepdfs = [pdfs.TruncatedExponentialPDF(lifetime0, maxT),
               pdfs.TruncatedExponentialPDF(lifetime1, maxT)]
    # Get vectors of expected fractional bin content.
    energyfracs = [energypdfs[0].binfractionvector(nEbins, (0,maxE)),
                   energypdfs[1].binfractionvector(nEbins, (0,maxE))]
    timefracs = [timepdfs[0].binfractionvector(nTbins, (0,maxT)),
                 timepdfs[1].binfractionvector(nTbins, (0,maxT))]
    # Take outer product of vectors to get 2-D array. Rows are time bins, columns
    # are energy bins. 
    fracs2D = [np.outer(timefracs[0],energyfracs[0]), 
               np.outer(timefracs[1],energyfracs[1])]

    # Make sure no bin has fewer than 'minevtsperbin' events.
    checkminbins(minevtsperbin, [nevents0,nevents1], energyfracs, timefracs, 
                 fracs2D)

    


# This finds the bins with the fewest expected events in both the 1-D and 2-D 
# cases.
def checkminbins(minevtsperbin, nevents, energyfracs, timefracs, twoDfracs):
    min1Dbin = min(np.min(nevents[0]*energyfracs[0] + nevents[1]*energyfracs[1]),
                   np.min(nevents[0]*timefracs[0] + nevents[1]*timefracs[1]))    
    min2Dbin = np.min(nevents[0]*twoDfracs[0] + nevents[1]*twoDfracs[1]) 

    if min1Dbin <= minevtsperbin:
        print 'Warning: there is a bin in the 1-D fit with %s expected events. Either increase event rate or make binning coarser. Exiting.'
        sys.exit(1)
    if min2Dbin <= minevtsperbin:
        print 'Warning: there is a bin in the 2-D fit with %s expected events. Either increase event rate or make binning coarser. Exiting.'
        sys.exit(1)


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
