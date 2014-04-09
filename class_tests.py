class firstclass:
    var=7
    def dispvar(self):
        print 'var: %s' % self.var

# Inherit from firstclass
class secondclass(firstclass):
    def dispvar(self):
        print '2nd var: %s' % self.var
    
