#! /usr/bin/env python
'''
run sigmoid fitting function, record error, plot best fitting sigmoid with raw data
'''
import os
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
import cmap.io.gct as gct
import cmap.util.tool_ops as tool_ops
import cmap.util.progress as progress
import cmap.analytics.fitting as fitting
import cmap.analytics.cluster as cluster
import cmap.util.queue as queue
import cmap.io.plategrp as grp
import cmap.plot.colors as colors

#run with all cell lines/ time points
# cellLst = ['PC3', 'A375', 'MCF7']
# timeLst = ['6H', '24H']
cell = 'PC3'
tim = '6H'
cellLine = cell
timeP = tim
refControl = 'pc' #use pc vs vc controled data
gctfile = glob.glob('/xchip/obelix/pod/brew/%s/PRISM001_%s_%s/by_pert_id_pert_dose/PRISM001_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
gctfile = gctfile[0]
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/ec50_test/%s_%s_%s' % (cell,timeP,refControl)
if not os.path.exists(work_dir):
	os.mkdir(work_dir)
gcto = gct.GCT() #make a gct object
gcto.read(gctfile)

cids = gcto.get_gctx_cid(gctfile)
rids = gcto.get_gctx_rid(gctfile)
pert_descs = gcto.get_column_meta('pert_desc')
doses = [float(x.split(':')[2]) for x in cids]
perts = [x.split(':')[1] for x in cids]
unique_perts = list(set(perts))


#evaluate a unique pert
unique_pert = 'VX-680'
cid_inds = [i for i,x in enumerate(cids) if unique_pert in x]
pert_doses = [doses[x] for x in cid_inds]
tmp_tup = zip(pert_doses, cid_inds)
tmp_tup.sort()
pert_doses, cid_inds = zip(*tmp_tup) #over write variables in order
pert_data = gcto.matrix[:,cid_inds]
EC50 = fitting.EC50(list(pert_doses),pert_data)
EC50.compute_ec50()

errorList = []
for iprobe in range(len(rids)):
	iprobe = 7
	t = EC50.concentrations
	y = EC50.responses[:,iprobe]
	# y = pert_data[iprobe,:]
	v = [EC50.min[iprobe], EC50.max[iprobe], EC50.ec50[iprobe], EC50.hill_coeficient[iprobe]]
	x = np.linspace(0, np.max(pert_doses), num=100)
	ySigFitLong = fitting.EC50._sigmoid(EC50,v,x)
	ySigFitShort = fitting.EC50._sigmoid(EC50,v,t)
	#calculate error for the best fit 
	error1 = fitting.EC50._error_fn(EC50,v,t,y)
	errorList.append(error1)

plt.plot(x,ySigFitLong)
plt.plot(t,ySigFitShort,'co')
plt.plot(t,y,'ro')
plt.xlabel('dose (um)')
plt.ylabel('z-score')
# plot the data
plt.plot(t, y, 'ro')
plt.plot(x,ySigFit)


low_pass = EC50.ec50 >= pert_doses[0]
high_pass = EC50.ec50 <= pert_doses[-1]
passing_inds = low_pass * high_pass
non_passing_inds = np.invert(passing_inds)
rational_rid_inds = np.arange(len(rids))
irrational_rid_inds = rational_rid_inds[non_passing_inds]
rational_ec50 = EC50.ec50[passing_inds]

#use best results to make sigmoid
# http://lister.dulci.duhs.duke.edu/~cliburn/summer-school/python/_build/html/model_fitting.html
# The Hill function generates an S-shaped curve, where  is the variable 
# of interest (time in our example), a is the baseline, a+b is the maximum or 
# saturated value, k is the value of t when the function has 50% of its 
# maximal value, and n, the Hill coefficient, determines how steep the slope 
# of the function is at k.

# t= doses
# a = minVal
# a+b = maxVal 
# k= ec50
# n= Hill Coef
def sigmoid_f(t, a, b, k, n):
  return a + b*(t**n)/(t**n + k**n)



sigOut = sigmoid_f(t, a, b, k, n)

### repeat with sigmoid function from fitting tool

t = pert_doses
y = pert_data[ipert,:]


def _sigmoid(v, x):
	'''
	This is the EC50 sigmoid function

	v is a vector of parameters:
	    v[0] = minimum allowed value
	    v[1] = maximum allowed value
	    v[2] = ec50
	    v[3] = Hill coefficient
	'''
	p_min, p_max, ec50, hill = v
	return p_min + ((p_max - p_min) /
	                (1+(x/ec50)**hill))






import pylab
import numpy
from scipy import optimize

def f(t, a, b, k, n):
    return a + b*(t**n)/(t**n + k**n)

def resid(p, y, t):
    a, b, k, n = p
    return y - f(t, a, b, k, n)

if __name__ == '__main__':
    # load the data
    # 4 such columns, we assign them to t (time), x1 (CFUs in experiment 1), x2 (CFUs in experiment 2) and d (dilurtion factor) respectively.
    # t, x1, x2, d = numpy.loadtxt('cfu.txt', unpack=True)
    #set t to dose
    t = pert_doses
    ipert = 0
    y = pert_data[ipert,:]


    # transform the observed CFU values
    y1 = numpy.log10(x1*10**d)
    y2 = numpy.log10(x2*10**d)

    y = numpy.concatenate([y1, y2])

    a0, b0, k0, n0 = 1, 1, 1, 1

    [a, b, k, n], flag  = optimize.leastsq(resid, [a0, b0, k0, n0], 
                                           args=(y1, t))

    print flag, a, b, k, n

    # plot the data
    pylab.plot(t, y1, 'ro')

    # plot the smooth model fit
    ts = numpy.linspace(t[0], t[-1], 100)
    pylab.plot(ts, f(ts, a, b, k, n))

    pylab.title('Fitted log growth curves of E. coli')
    pylab.text(20, 6.8, 
               r'$f(t) = a +b \frac{t^n}{t^n + k^n}$'  ,
               fontsize=18)
    pylab.text(20, 6.5, 'a=%.2f, b=%.2f' % (a, b))
    pylab.text(20, 6.2, 'k=%.2f, n=%.2f' % (k, n))
    pylab.xlabel('Time (mins)')
    pylab.ylabel('$\log_{10}(\mathrm{CFUs})$')

    # show is necessary to display the plot when
    # not in interactive mode
    pylab.show()

