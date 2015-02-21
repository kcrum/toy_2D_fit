import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

np.random.seed(19)

# Both of the PDFs below have nearly identical functionality and could 
# therefore inherit from an intermediate parent class. If you add any more PDFs
# that have the same forms as those here, you should consider adding a base 
# class between the final physics PDFs and their (grand)parent class 
# "rv_continuous."

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
    # Overload default _pdf that comes from rv_continuous. It seems like you
    # must normalize the pdf yourself. (Test by removing normfactor, then 
    # trying to call sampleplot()).
    def _pdf(self,x):
        return -self.normfactor*x*(x - self.endpoint)
    # You can change the endpoint whenever you like.
    def setendpoint(self,newendpoint):
        self.endpoint = float(newendpoint)
        self.b = newendpoint
        self.normfactor = 6./(self.endpoint**3)
    # Return the fractional bin occupancy vector for a given binning
    def binfractionvector(self, nbins, binrange):
        if not isinstance(binrange, tuple) or len(binrange) != 2:
            print 'binfractionvector expects binrange to be a 2-tuple. Exiting!'
            sys.exit(1)
        binwidth = float(binrange[1] - binrange[0]) / nbins
        binfracvec = np.zeros(nbins)
        prevcdf = 0.
        for i in xrange(1,nbins+1): # This counts inclusively over {1,nbins}
            uppercdf = self.cdf(binrange[0] + i*binwidth)
            binfracvec[i-1] = uppercdf - prevcdf
            prevcdf = uppercdf
        # Simple check
        if sum(binfracvec) != 1.0: print 'Warning: sum of binfracvec = %s' % \
                sum(binfracvec)
        return binfracvec
    # This is for debugging. It plots a histogram of 100 draws from the pdf.
    def sampleplot(self,ndraws=100,nbins=10):
        hist, bins = np.histogram(self.rvs(size=ndraws), bins=nbins,
                                  range=(0,self.endpoint))
        histarea = float(ndraws)*(float(self.endpoint)/nbins)
        barwidth = 0.9*(bins[1]-bins[0])
        bincenters = (bins[:-1] + bins[1:]) / 2
        plt.bar(bincenters, hist, width = barwidth)
        xvals = np.linspace(0,self.endpoint,100)
        plt.plot(xvals, histarea*self._pdf(xvals))
        plt.xlabel("Energy [MeV]")
        plt.show()
    # This is the same as sampleplot, except now we're using the PyPlot 
    # histogram interface instead of the numpy histogram passed to a bar graph.
    def sampleplot2(self,ndraws=100,nbins=10):
        myhist = plt.hist(self.rvs(size=ndraws), bins=nbins, 
                        range=(0,self.endpoint), histtype='stepfilled', 
                        alpha=0.8, color='green')
        maxval = max(myhist[0]) # myhist is an 2-d array with bin contents at
        # [0] and bin edges at [1]. We find the maxval here to set the y-axis 
        # limit.
        plt.ylim([0,maxval + 0.05*maxval])
        histarea = float(ndraws)*(float(self.endpoint)/nbins)
        xvals = np.linspace(0,self.endpoint,100)
        plt.plot(xvals, histarea*self._pdf(xvals))
        plt.xlabel("Energy [MeV]")
        plt.show()
    # Yet another plot test.
    def sampleplot3(self,ndraws=100,nbins=10):
        fig, ax = plt.subplots()
        ax.hist(self.rvs(size=ndraws), bins=nbins, range=(0,self.endpoint),
                 histtype='bar', alpha=0.8, color='green')
        histarea = float(ndraws)*(float(self.endpoint)/nbins)
        xvals = np.linspace(0,self.endpoint,100)
        ax.plot(xvals, histarea*self._pdf(xvals))
        plt.xlabel("Energy [MeV]")
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
        self.maxT = float(maxT)
        self.lifetime = float(lifetime)
        self.normfactor = 1./(self.lifetime*(1.-np.exp(-self.maxT/self.lifetime)))
        # Make sure to call the base/parent class's ctor, too.
        super(TruncatedExponentialPDF, self).__init__(a=0, b=maxT)
    # Overload default _pdf inherited from rv_continuous. 
    def _pdf(self, T):
        return self.normfactor*np.exp(-T/self.lifetime)
    # You can change maxT
    def setmaxT(self,newmaxT):
        self.maxT = float(newmaxT)
        self.b = newmaxT
        self.normfactor = 1./(self.lifetime*(1-np.exp(-self.maxT/self.lifetime)))
    # You can change lifetime
    def setlifetime(self,newlifetime):
        self.lifetime = float(newlifetime)
        self.normfactor = 1./(self.lifetime*(1-np.exp(-self.maxT/self.lifetime)))
    # Return the fractional bin occupancy vector for a given binning
    def binfractionvector(self, nbins, binrange):
        if not isinstance(binrange, tuple) or len(binrange) != 2:
            print 'binfractionvector expects binrange to be a 2-tuple. Exiting!'
            sys.exit(1)
        binwidth = float(binrange[1] - binrange[0]) / nbins
        binfracvec = np.zeros(nbins)
        prevcdf = 0.
        for i in xrange(1,nbins+1): # This counts inclusively over {1,nbins}
            uppercdf = self.cdf(binrange[0] + i*binwidth)
            binfracvec[i-1] = uppercdf - prevcdf
            prevcdf = uppercdf
        # Simple check
        if sum(binfracvec) != 1.0: print 'Warning: sum of binfracvec = %s' % \
                sum(binfracvec)
        return binfracvec
    # This is for debugging. It will plot a histogram of 100 draws from the pdf.
    def sampleplot(self,ndraws=100,nbins=10):        
        hist, bins = np.histogram(self.rvs(size=ndraws), bins=nbins,
                                  range=(0,self.maxT))
        histarea = float(ndraws)*(float(self.maxT)/nbins)
        bincenters = (bins[:-1] + bins[1:]) / 2
        barwidth = 1.0*(bins[1]-bins[0])
        plt.bar(bincenters, hist, align='center', alpha=0.8, width=barwidth, 
                color='green')
        xvals = np.linspace(0,self.maxT,100)
        plt.plot(xvals, histarea*self._pdf(xvals))
        plt.xlabel(r"$\Delta T$")
        plt.show()
        

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
    
    EnergyPDF = ParabolicPDF(endpoint)
    DeltaTPDF = TruncatedExponentialPDF(lifetime, maxT)

    EnergyPDF.sampleplot3()
    DeltaTPDF.sampleplot()
