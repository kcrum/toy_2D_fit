#Toy 2-D Fit Tests

This contains code that runs fits to fake data where each event has an energy and deltaT. The fake data comes from two "isotopes," both of which have different energy end points and exponential decay time constants. The fit will attempt to measure the normalizations of these two spectra (so there will be two free parameters).

I want to test whether simultaneously fitting both of the 1-D distributions is unbiased; I also want to observe the behavior of the chi^2 for many fake fits. I will also fit the (presumably correct) 2-D distribution, checking again for unbiased estimators and proper chi^2 behavior. 

## Running fake fits

The fake data in the fits done here are pulled from energy and time distributions, much like if one were observing Li-9 or He-8 decays. The energy spectrum form is specified in physicsPDFs.py by the "ParabolicPDF" class, and it is a parabola of the form ~x*(x-endpoint), where the user specifies the endpoint. The time spectrum is an exponential from the "TruncatedExponentialPDF" class with user-specified lifetime, as well as a user-specified maximum time. We generate fake data pulling from energy and time distributions for two 'isotopes,' each with its own spectral endpoint and lifetime. 

Fake data is generated and then fit by the "mainloop" function in toy_2D_fits.py. Four different fits are run to extract the normalization of the two isotopes:

- A least-squares fit to the 4 energy bins and 4 time bins, using the Pearson treatment for statistical errors.
- A least-squares fit to the combined 16 (4x4) energy and time bins, using the Pearson treatment for statistical errors.
- A maximum likelihood fit to the 4 energy bins and 4 time bins, using a Poisson pdf for each of the 8 bins.
- A maximum likelihood fit to the combined 16 (4x4) energy and time bins, using a Poisson pdf for each of the 16 bins.

By default the endpoints are 12 and 8 (arb. units) for isotope0 and isotope1, respectively. The four energy bins span the range from 0 to 12. Isotope0 has a default lifetime of 260, while isotope1 has a default lifetime of 170 (arb. units). The time bins span 0 to 260. 

By making the following call from the command line:
```
python toy_2D_fits.py 1000 1000 100 toyfits_1000exp_1000n0_100n1.txt
```
you should be able to reproduce the output in toyfits_1000exp_1000n0_100n1.txt. This contains the fits of 1,000 fake experiments, each having 1,000 events pulled from isotope0 and 100 events pulled from isotope1.