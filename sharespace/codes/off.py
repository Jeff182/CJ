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
D['AV18']  ={'p':[0.098,0.345,0.048] ,'c':'r','ls':'-'}   
D['CDBonn']={'p':[0.138,0.373,0.058] ,'c':'g','ls':'--'}
D['WJC1']  ={'p':[0.132,-0.253,0.152],'c':'k','ls':'-.'}
D['WJC2']  ={'p':[0.076,0.313,0.033] ,'c':'b','ls':':'}

f=lambda x,p:p[0]*(x-p[1])*(x-p[2])*(1+p[1]-x)
X=np.linspace(0,1,100)
ax=py.subplot(111)
for k in ['AV18','CDBonn','WJC1','WJC2']:
  ax.plot(X,f(X,D[k]['p']),\
    color=D[k]['c'],\
    lw=2.0,\
    ls=D[k]['ls'],\
    label=tex(k)
    )
ax.set_ylim(-0.025,0.04)
ax.set_ylabel(r'$\delta f^N$',size=25)
ax.set_xlabel('$x$',size=25)
ax.legend(frameon=0,loc=2,fontsize=20)
py.tick_params(axis='both',labelsize=20)
py.tight_layout()
py.savefig('gallery/off_shell.pdf')