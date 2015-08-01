#!/usr/bin/env python
import sys,os
sys.path.insert(1,'../../')
import numpy as np
import pylab as py
import lhapdf
from tools import tex, fill_between
import matplotlib.gridspec as gridspec
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
from master import FITPACK, COMPOSER
  
# define x values
X=np.linspace(1e-3,0.9,1000)

# initialize composer
CJ=COMPOSER(name='CJ12min',X=X)
#CJ.map_X() # optional

# get glue 
# note: the output is a dictionary where
# key='f0' -> central PDF
# key='df' -> symmetric error
# key='df+' -> asymmetric error (+)
# key='df0' -> asymmetric error (-)
u=CJ.get_xpdf(flav='u',Q2=10.0)
print u.keys()


# these lines add latex typefont.
# is kind of slow, so uncomment while a plot is been coded
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)  

# create matplotlib instance
ax=py.subplot(111)

# plot error band
T=10.0
UP=u['f0']+T*u['df+']
DO=u['f0']-T*u['df-']
fill_between(X,DO,UP,ax=ax,label=tex('u'),facecolor='r',edgecolor='r',alpha=0.5,hatch=None)

# plot central
ax.plot(X,u['f0'],'r-',label=tex('u'))

# makeup
ax.legend(loc=3,frameon=0,fontsize=20)
ax.semilogx()
ax.set_xlabel(tex('x'),size=20)
ax.set_ylabel(tex('xPDF(x)'),size=20)
py.tight_layout()
py.savefig('plots/advanced.pdf')


