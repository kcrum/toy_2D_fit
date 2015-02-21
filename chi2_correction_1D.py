import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st

# This will take a pval evaluated from a chi^2 assuming 'bad_dof' degrees of 
# freedom and then re-evaluate pval using 'good_dof' degrees of freedom. 
def pval_correction(pval, bad_dof = 6, good_dof = 5):
    # The probability of getting a chi^2 above a certain value x is the same as
    # 1 - cdf(x), which is also called the survival function of x: sf(x). Thus
    # given a pval, we can find the original chi^2 value x by chi2.isf(x), 
    # which gives the inverse of the sf. 
    orig_chi2 = st.chi2.isf(pval, df = bad_dof)
    return st.chi2.sf(orig_chi2, df = good_dof)


#########################################################################
# Make histogram of p-vals for 1-D and 2-D chi^2 fits.
def pval_distributions(filepath = 'toy_fits_1000exp_1000n0_100n1.txt'):

    data = pd.read_csv(filepath)
    # Correct 1-D pvals
    data['corrected_pval_1D'] = data.pval_1D.apply(pval_correction, args=(6,5))

    plt.hist(data.pval_1D, alpha=0.9, hatch='o',
             label=r'n$_{d.o.f.}$ = 6')
    plt.hist(data.corrected_pval_1D, alpha=0.5, hatch='/',
             label=r'n$_{d.o.f.}$ = 5 (corrected)')
    # Uncomment next two lines if you also want a hist of the 2D pvals.
    #plt.hist(data.pval_2D, alpha=0.5, hatch='/',
    #         label='2-dimensional (d.o.f. = 14)')

    plt.legend(loc=2, fontsize='x-large')
    plt.xlabel(r"P-value evaluated from $\chi^2_{min}$", fontsize='x-large')
    plt.title('Histogram of P-values assuming 5 and 6 degrees of freedom')

    plt.show()
