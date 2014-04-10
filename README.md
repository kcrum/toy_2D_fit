#Toy 2-D Fit Tests

This will contain code that runs fits to fake data where each event has an energy and deltaT. The fake data will come from two "isotopes," both of which have different energy end points and time constants. The fit will attempt to measure the normalizations of these two spectra (so there will be two free parameters).

I want to test whether simultaneously fitting both of the 1-D distributions is unbiased; I also want to observe the behavior of the chi^2 for many fake fits. I will also fit the (presumably correct) 2-D distribution, checking again for unbiased estimators and proper chi^2 behavior. 
