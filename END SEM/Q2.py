import lib

# fitting by modified chebyshev fucntion
# modified fitting code is at the end of the library
# print(lib.fitting('esem4fit.txt', '\t', 6, plot=True))
print(lib.fitting_modified('esem4fit.txt', '\t', 6, plot=True))

'''
{'Coefficients': [[0.07003196671971398], [0.004301685837864321], [-0.010166710608800469], [0.013083743602879215], [0.11411855049286528], [-0.006726972223322477], [-0.012384559712646195]],
'chi2': 0.8808496213958824}
[Finished in 1.658s]

'''
