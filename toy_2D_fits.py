import sys, os
import scipy.stats as st

# Inherits  from scipy.stats
class parabolic_pdf(st.rv_continuous):
    """
    Create a pdf using the function -(6/endpoint**3)*x*(x - endpoint) on the
    domain 0 <= x <= endpoint. 
    """
    endpoint = 1.
    def setEndpoint(self, endpoint):
        if endpoint < 0:
            print "Endpoint of parabolic_pdf() must by postive. User passed %s. \
            Exiting!" % endpoint
            sys.exit()
        else: self.endpoint = endpoint
    # Here you're redefining rv_continuous's _pdf method. Ideally you'd force
    # a = 0 and b = 1. Is that possible?
    def _pdf(self, x, endpoint):        
        self.endpoint = endpoint
        return -(6/self.endpoint**3)*x*(x - self.endpoint)
