#!/usr/bin/env python
import sys,os
sys.path.insert(1,'../../')
import numpy as np
import pylab as py
import tools
from tools import tex,plot_band,fill_between
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Times-Roman']})
rc('text',usetex=True)

D={}
D['AV18']  ={'h':[-3.0094,1.7526,-2.0895] ,'c':'r','ls':'-'}   
D['CDBonn']={'h':[-2.9851,1.7564,-2.0856] ,'c':'g','ls':'--'}
D['WJC1']  ={'h':[-3.2169,1.8225,-2.0844],'c':'k','ls':'-.'}
D['WJC2']  ={'h':[-3.0403,1.7605,-2.0898] ,'c':'b','ls':':'}

f=lambda x,h: h[0] * x**h[1] * (1+h[2]*x)
X=np.linspace(0,1,100)
ax=py.subplot(111)
for k in ['AV18','CDBonn','WJC1','WJC2']:
  ax.plot(X,f(X,D[k]['h']),\
    color=D[k]['c'],\
    lw=2.0,\
    ls=D[k]['ls'],\
    label=tex(k)
    )
ax.set_ylim(-0.5,2)
ax.set_ylabel(r'$C_{\rm HT}$',size=30)
ax.set_xlabel('$x$',size=30)
ax.axhline(0,color='k',ls='-',alpha=0.2)
ax.legend(frameon=0,loc=2,fontsize=25)
py.tick_params(axis='both',labelsize=25)
py.tight_layout()
py.savefig('gallery/ht.pdf')
