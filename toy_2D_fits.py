import physicsPDFs as pdfs
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
        
#def Espec_data(endpoint, nevents):
    
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
