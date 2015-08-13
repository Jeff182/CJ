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
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
import pandas as pd

class F2ratio(object):

  def __init__(self):
    self.load_data()
    self.make_plot()

  def load_data(self):
    F=open('../CJ15data/FITS_150701_LO_NLO/calc_NLO_KP_AV18.out')
    L=F.readlines()
    F.close()
    L=[l.strip() for l in L]
    L=[l for l in L if l!='']
    p=[l.split() for l in L if 'slac_p_reb' in l if len(l.split())<25]
    d=[l.split() for l in L if 'slac_d_reb' in l if len(l.split())<25]
    for l in L:
      if all([k in l for k in ['X','Q2','W2','THEORY']]): H=l.replace('_w','').split()
   
    D={}
    D['p']=pd.DataFrame(p,columns=H)
    D['d']=pd.DataFrame(d,columns=H)
    D['p']=D['p'].convert_objects(convert_numeric=True)
    D['d']=D['d'].convert_objects(convert_numeric=True)
    self.D=D

  def make_plot(self):
    # retrieve data
    D=self.D
    Dp=D['p']
    Dd=D['d']
    ax=py.subplot(111)
    ax.plot(Dp.X,Dp.Q2,'r.')
    ax.set_xlabel('X')
    ax.set_ylabel('Q2')
    py.savefig('gallery/tmp_XQ.pdf')

if __name__=='__main__':

  F2ratio()


