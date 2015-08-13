#!/usr/bin/env python
import sys,os
sys.path.insert(1,'../../')
import numpy as np
import pylab as py
import tools
from tools import tex,plot_band,fill_between
import lhapdf
import matplotlib.gridspec as gridspec
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Times-Roman']})
rc('text',usetex=True)
import pandas as pd

class F2ratio(object):

  def __init__(self):
    self.load_data()
    self.make_plotI()
    self.make_plotII()

  def _load_data(self,name):
    D=self.D
    F=open('../CJ15data/FITS_150701_LO_NLO/calc_NLO_KP_'+name+'.out')
    L=F.readlines()
    F.close()
    L=[l.strip() for l in L]
    L=[l for l in L if l!='']
    data=[l.replace('*','').split() for l in L if 'test_DN' in l if len(l.split())<25]
    for l in L:
      if all([k in l for k in ['X','Q2','W2','THEORY']]):
        H=l.replace('_w','').replace(r'chi^2','').split()
    DF=pd.DataFrame(data,columns=H)
    DF=DF.convert_objects(convert_numeric=True)
    D[name]=DF

  def load_data(self):
    self.D={}
    self._load_data('AV18')
    self._load_data('CDBONN')
    self._load_data('WJC1')
    self._load_data('WJC2')

  def make_plotI(self):
    # retrieve data
    D=self.D

    kmap={}
    kmap['AV18']   = {'c':'r','ls':'-'}
    kmap['CDBONN'] = {'c':'g','ls':'--'}
    kmap['WJC1']   = {'c':'k','ls':'-.'}
    kmap['WJC2']   = {'c':'b','ls':':'}

    ax=py.subplot(111)
    for k in ['AV18','CDBONN','WJC1','WJC2']:
      DF=D[k]
      DF=DF[DF.Q2==10]
      if k=='CDBONN':
        label='CDBonn'
      else:
        label=k
      cls=kmap[k]['c']+kmap[k]['ls']
      ax.plot(DF.X,DF.THEORY,cls,lw=2.0,label=tex(label))

    ax.set_xlabel('$x$',size=25)
    ax.set_ylabel(r'$F_2^d\, /\, F_2^N$',size=25)
    ax.set_ylim(0.97,1.08)
    ax.axhline(1,color='k',ls='-',alpha=0.2)

    ax.legend(frameon=0,loc=2,fontsize=22)
    py.tick_params(axis='both',labelsize=22)
    py.tight_layout()
    py.savefig('gallery/F2d_F2_I.pdf')
    py.close()

  def make_plotII(self):
    # retrieve data
    D=self.D

    kmap={}
    kmap['Q2 = 2']   = {'c':'r','ls':'-'}
    kmap['Q2 = 5']   = {'c':'g','ls':'--'}
    kmap['Q2 = 10']  = {'c':'b','ls':'-.'}
    kmap['Q2 = 100'] = {'c':'k','ls':':'}


    ax=py.subplot(111)
    DF=D['AV18']
    for Q2 in [2,5,10,100]:
      k='Q2 = %d'%Q2
      Q2=float(k.split('=')[1])
      DF=D['AV18'][D['AV18'].Q2==Q2]
      cls=kmap[k]['c']+kmap[k]['ls']
      ax.plot(DF.X,DF.THEORY,cls,lw=2.0,label=r'$Q^2=%0.0f~{\rm GeV}^2$'%Q2)

    ax.set_xlabel('$x$',size=25)
    ax.set_ylabel(r'$F_2^d\, /\, F_2^N$',size=25)
    ax.set_ylim(0.97,1.08)
    ax.axhline(1,color='k',ls='-',alpha=0.2)

    ax.legend(frameon=0,loc=2,fontsize=22)
    py.tick_params(axis='both',labelsize=22)
    py.tight_layout()
    py.savefig('gallery/F2d_F2_II.pdf')

if __name__=='__main__':

  F2ratio()
