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

def mainloop(nevents1, nevents2, endpoint1=12.0, endpoint2=8.0, lifetime1=260,
             lifetime2=170, nEbins=4, nTbins=4):
    maxT = max(lifetime1, lifetime2)
    maxE = max(endpoint1, endpoint2)
    energypdfs = [pdfs.ParabolicPDF(endpoint1), pdfs.ParabolicPDF(endpoint2)]
    timepdfs = [pdfs.TruncatedExponentialPDF(lifetime1, maxT),
               pdfs.TruncatedExponentialPDF(lifetime2, maxT)]
    energyfracs = [energypdfs[0].binfractionvector(nEbins, (0,maxE)),
                   energypdfs[1].binfractionvector(nEbins, (0,maxE))]
    timefracs = [timepdfs[0].binfractionvector(nTbins, (0,maxT)),
                 timepdfs[1].binfractionvector(nTbins, (0,maxT))]
    # Take outer product of vectors to get 2-D array. 
    fracs2D = [np.outer(timefracs[0],energyfracs[0]), 
               np.outer(timefracs[1],energyfracs[1])]

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
