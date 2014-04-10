import sys, os
import numpy as np
import scipy.stats as st

class FirstClass:
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

