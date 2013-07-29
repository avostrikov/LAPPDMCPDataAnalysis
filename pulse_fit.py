import numpy as np
from scipy.stats import lognorm
from scipy.optimize import fmin

def fit_pulse(x,y,x0):
    def squares(args):
        c, s, mu, std = args
        #print 'args = ', args
        fy = c*lognorm.pdf(x,s,mu,std)
        fy = np.nan_to_num(fy)

        value = np.sum((fy - y)**2)
        #print 'value = ', value
        return value

    xopt, fopt, niter, funcalls, warnflag = fmin(squares,x0,full_output=True,
                                                 disp=False)

    if warnflag:
        raise Exception('fit did not converge')

    return xopt
