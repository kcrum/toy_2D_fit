import sys, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

# You can make your own PDF by inheriting from rv_continuous in scipy.stats.
# Since the parabolic PDF I'm defining here isn't sensible as a PDF outside of
# [0,endpoint], I'm setting the inherited "a" to 0 and "b" to endpoint.
class ParabolicPDF(st.rv_continuous):
    """
    Create a pdf characterized by the function:
       -(6/endpoint**3)*x*(x - endpoint)
    on the domain 0 <= x <= endpoint. This inherits from the
    scipy.stats.rv_continuous class.
    """
    def __init__(self, endpoint):
        if endpoint <= 0:
            print 'Endpoint %s must be greater than 0. Exiting!' % endpoint
            sys.exit()
        self.endpoint = float(endpoint)
        self.normfactor = 6/(self.endpoint**3)
        # Make sure to call the base/parent class's ctor, too.
        super(ParabolicPDF, self).__init__(a=0, b=endpoint)
    # Overload default _pdf that comes from rv_continuous. It seems like you must
    # normalize the pdf yourself. (Test by removing normfactor, then trying to
    # call sampleplot()).
    def _pdf(self,x):
        return -self.normfactor*x*(x - self.endpoint)
    # You can change the endpoint whenever you like.
    def setendpoint(self,newendpoint):
        self.endpoint = float(newendpoint)
        self.b = newendpoint
        self.normfactor = 6./(newendpoint**3)
    # This is for debugging. It will plot a histogram of 100 draws from the pdf.
    def sampleplot(self,ndraws=100,nbins=10):
        hist, bins = np.histogram(self.rvs(size=ndraws), bins=nbins,
                                  range=(0,self.endpoint))
        histarea = float(ndraws)*(float(self.endpoint)/nbins)
        barwidth = 0.7*(bins[1]-bins[0])
        bincenters = (bins[:-1] + bins[1:]) / 2
        plt.bar(bincenters, hist, align='center', width = barwidth)
        # Alternatively, you could replace the previous line with:
        #   fig, ax = plt.subplots()
        #   ax.bar(bincenters, hist, align='center', width = barwidth)
        xvals = np.linspace(0,self.endpoint,100)
        plt.plot(xvals, histarea*self._pdf(xvals))
        plt.show()

class TruncatedExponentialPDF(st.rv_continuous):
    """
    Create an exponential PDF on a shortened range (x in [0, maxT]). This
    inherits from the scipy.stats.rv_continuous class.
    """
    def __init__(self, lifetime, maxT):
        if maxT <= 0: 
            print 'maxT %s must be greater than 0. Exiting!' % maxT
            sys.exit()
        if lifetime <= 0: 
            print 'lifetime %s must be greater than 0. Exiting!' % lifetime
            sys.exit()
        self.maxT = maxT
        self.lifetime = lifetime
        self.normfactor = 1./(lifetime*(1 - np.exp(-maxT/lifetime)))
        # Make sure to call the base/parent class's ctor, too.
        super(TruncatedExponentialPDF, self).__init__(a=0, b=maxT)
    # Overload default _pdf inherited from rv_continuous. 
    def _pdf(self, T):
        return self.normfactor*np.exp(-T/self.lifetime)
    # You can change maxT
    def setmaxT(self,newmaxT):
        self.maxT = newmaxT
        self.b = newmaxT
        self.normfactor = 1./(self.lifetime*(1 - np.exp(-newmaxT/self.lifetime)))
    # You can change lifetime
    def setlifetime(self,newlifetime):
        self.lifetime = newlifetime
        self.normfactor = 1./(newlifetime*(1 - np.exp(-self.maxT/newlifetime)))
    # This is for debugging. It will plot a histogram of 100 draws from the pdf.
    def sampleplot(self,ndraws=100,nbins=10):        
        hist, bins = np.histogram(self.rvs(size=ndraws), bins=nbins,
                                  range=(0,self.maxT))
        histarea = float(ndraws)*(float(self.maxT)/nbins)
        barwidth = 0.7*(bins[1]-bins[0])
        bincenters = (bins[:-1] + bins[1:]) / 2
        plt.bar(bincenters, hist, align='center', width = barwidth)
        xvals = np.linspace(0,self.maxT,100)
        plt.plot(xvals, histarea*self._pdf(xvals))
        plt.show()
        

#def Espec_data(endpoint, nevents):
    
