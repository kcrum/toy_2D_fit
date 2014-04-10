import sys, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

# You can make your own PDF by inheriting from rv_continuous in scipy.stats.
# Since the parabolic PDF I'm defining here isn't sensible as a PDF outside of
# [0,endpoint], I'm setting the inherited "a" to 0 and "b" to endpoint.
class ParabolicPDF(st.rv_continuous):
    """
    Create a pdf using the function -(6/endpoint**3)*x*(x - endpoint) on the
    domain 0 <= x <= endpoint. 
    """
    def __init__(self, endpoint):
        if endpoint <= 0:
            print 'Endpoint %s must be greater than 0. Exiting!' % endpoint
            sys.exit()
        self.endpoint = endpoint
        # Make sure to call the base/parent class's __init__, too.
        super(ParabolicPDF, self).__init__(a=0, b=endpoint)
    # Overload default _pdf that comes from rv_continuous. It seems like you must
    # normalize yourself. (Test by removing 6/self.endpoint**3, then trying to
    # call sampleplot().
    def _pdf(self,x):
        return -(6/self.endpoint**3)*x*(x - self.endpoint)
    # You can change the endpoint whenever you like.
    def setendpoint(self,newend):
        self.endpoint = newend
        self.b = newend
    # This is for debugging. It will plot a histogram of 100 draws from the pdf.
    def sampleplot(self,ndraws=100,nbins=10):
        hist, bins = np.histogram(self.rvs(size=ndraws), bins=nbins,
                                  range=(0,self.endpoint))
        barwidth = 0.7*(bins[1]-bins[0])
        bincenters = (bins[:-1] + bins[1:]) / 2
        plt.bar(bincenters, hist, align='center', width = barwidth)
        # Alternatively, you could replace the previous line with:
        #   fig, ax = plt.subplots()
        #   ax.bar(bincenters, hist, align='center', width = barwidth)
        plt.show()

#def Espec_data(endpoint, nevents):
    
