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

class LASY(object):

  def __init__(self):
    self.load_data()
    self.define_plotter_geometry()
    self.make_plot()

  def load_data(self):
    F=open('../CJ15data/CJ15_NLO_KP_AV18/CJ15_Wasym.dat')
    L=F.readlines()
    F.close()
    L=[l.strip() for l in L]
    L=[l for l in L if l!='']
    H=L[0].split()
    L=[l.split() for l in L[1:]]
    DF=pd.DataFrame(L,columns=H)
    DF=DF.convert_objects(convert_numeric=True)

    D={}
    D['cdfLasy05'] = DF[DF.ITYPE=='cdfLasy05']
    D['d0Lasy_e15']= DF[DF.ITYPE=='d0Lasy_e15']
    D['d0Lasy13']  = DF[DF.ITYPE=='d0Lasy13']
    self.D=D

  def define_plotter_geometry(self):

    ncols=2
    nrows=2
    py.figure(figsize=(ncols*2,nrows*2))
    gs = gridspec.GridSpec(nrows,ncols)
    gs.update(left=0.13,right=0.98,wspace=0.4,hspace=0.3,\
      top=0.98,bottom=0.12)
    
    AX={}
    AX['cdfLasy05']  = py.subplot(gs[0,0])
    AX['d0Lasy_e15'] = py.subplot(gs[0,1])
    AX['d0Lasy13']   = py.subplot(gs[1,0])
    self.AX=AX

  def plot_dataset(self,dataset,T=10):
    k=dataset
    D=self.D
    ax=self.AX[k]
    data=D[k]['DATA']
    derr=D[k]['DERROR']
    theory=D[k]['THEORY']
    terr=D[k]['ERROR']
    y=D[k]['Y']
    p1,=ax.plot(y,theory,'r-')
    p2=fill_between(y,theory-terr*T,theory+terr*T,
      ax=ax,
      facecolor='yellow',
      edgecolor='yellow')
    p3=ax.errorbar(y,data,yerr=derr,fmt='k.')
    return (p2,p1),p3

  def make_plot(self):
    AX=self.AX
    for k in AX.keys():    
      p21,p3=self.plot_dataset(k)
      AX[k].set_xlabel(r'$y_W$',size=10)
      AX[k].set_ylabel(tex(k.replace('_','')),size=10)

    ax=AX['cdfLasy05']
    ax.legend([p21,p3],[tex('CJ15'),tex('data')]\
      ,frameon=0,loc=3,fontsize=10,numpoints=1)
   
    #ax.text(0.5,0.8,tex('nrep=%d'%nrows),transform=ax.transAxes,size=20)
    py.savefig('gallery/Lasy.pdf')
    py.close()

if __name__=='__main__':

  LASY()
