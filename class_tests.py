import sys, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

# Always have classes inherit from object if they aren't already extending 
# another class.
class FirstClass(object):
    var1 = 7
    def __init__(self, val):
        self.var2 = val
    def dispvar(self):
        print 'var1: %s' % self.var1
        print 'var2: %s' % self.var2

# Inherit from firstclass
class SecondClass(FirstClass):
    def dispvar(self):
        print '2nd var: %s' % self.var2
    
class MyPDF(st.rv_continuous):
    def __init__(self, endpoint):
        if endpoint <= 0:
            print 'Endpoint %s must be greater than 0. Exiting!' % endpoint
        self.endpoint = endpoint
        print 'Endpoint set to %s' % self.endpoint
        super(MyPDF, self).__init__(a=0, b=endpoint)
        
    def setendpoint(self,newend):
        self.endpoint = newend
        self.b = newend
        print "a: %s" % self.a
        print "b: %s" % self.b
    def _pdf(self,x):
        return -(6/self.endpoint**3)*x*(x - self.endpoint)

# Parabolic PDF class with several histogram/plotting variations
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
    # Return the fractional bin occupancy vector for a given binning
    def binfractionvector(self, nbins, binrange=(0,self.endpoint)):
        binwidth = float(binrange[1] - binrange[0]) / nbins
        binfracvec = []
        prevcdf = 0
        for i in xrange(1,nbins+1): # This counts inclusively over {1,nbins}
            uppercdf = self.cdf(binrange[0] + i*binwidth)
            binfracvec.append(uppercdf - prevcdf)
            prevcdf = uppercdf

    # This is for debugging. It will plot a histogram of 100 draws from the pdf.
    def sampleplot(self,ndraws=100,nbins=10):
        hist, bins = np.histogram(self.rvs(size=ndraws), bins=nbins,
                                  range=(0,self.endpoint))
        histarea = float(ndraws)*(float(self.endpoint)/nbins)
        barwidth = 0.9*(bins[1]-bins[0])
        bincenters = (bins[:-1] + bins[1:]) / 2
        plt.bar(bincenters, hist, width = barwidth)
        xvals = np.linspace(0,self.endpoint,100)
        plt.plot(xvals, histarea*self._pdf(xvals))
        plt.show()
    # This is the same as sampleplot, except now we're using the PyPlot histogram
    # interface instead of the numpy histogram passed to a bar graph.
    def sampleplot2(self,ndraws=100,nbins=10):
        myhist = plt.hist(self.rvs(size=ndraws), bins=nbins,
                        range=(0,self.endpoint), histtype='stepfilled',
                        alpha=0.8, color='green')
        maxval = max(myhist[0]) # myhist is an 2-d array with bin contents at [0]
        # and bin edges at [1].
        plt.ylim([0,maxval + 0.05*maxval])
        histarea = float(ndraws)*(float(self.endpoint)/nbins)
        xvals = np.linspace(0,self.endpoint,100)
        plt.plot(xvals, histarea*self._pdf(xvals))
        plt.show()
    # Yet another plot test.
    def sampleplot3(self,ndraws=100,nbins=10):
        fig, ax = plt.subplots()
        ax.hist(self.rvs(size=ndraws), bins=nbins, range=(0,self.endpoint),
                 histtype='bar', alpha=0.8, color='green')
        histarea = float(ndraws)*(float(self.endpoint)/nbins)
        xvals = np.linspace(0,self.endpoint,100)
        ax.plot(xvals, histarea*self._pdf(xvals))
        plt.show()
