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
from  scipy.interpolate import  interp1d

class WASY(object):

  def __init__(self):
    self.load_data()
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
    D['CDF_Wasy'] = DF[DF.ITYPE=='CDF_Wasy']
    D['D0_Wasy']  = DF[DF.ITYPE=='D0_Wasy']
    self.D=D

  def make_plot(self):
    D=self.D
    ax=py.subplot(111)

    Y=D['D0_Wasy']['Y']
    T=D['D0_Wasy']['THEORY']
    ET=D['D0_Wasy']['ERROR']

    iT= interp1d(Y,T,kind='cubic')
    iET= interp1d(Y,ET, kind='cubic')
    Y=np.linspace(np.amin(Y),np.amax(Y),100)

    T=10
    p2=fill_between(Y,iT(Y)-iET(Y)*T,iT(Y)+iET(Y)*T,
      ax=ax,
      facecolor='yellow', 
      edgecolor='yellow')
    p1,=ax.plot(Y,iT(Y),'r-')
    H=[(p2,p1)]       
    L=[tex('CJ15')]

    dmap={}
    dmap['CDF_Wasy'] = {'color':'g','marker':'o'}
    dmap['D0_Wasy']  = {'color':'b','marker':'^'}

    for k in D.keys():
      color=dmap[k]['color']
      marker=dmap[k]['marker']
      markersize=4
      p3=ax.errorbar(D[k]['Y'],D[k]['DATA'],\
        yerr=D[k]['DERROR'],fmt=color+marker,mfc=color,mec=color,\
        markersize=markersize,zorder=1,alpha=0.9)
      H.append(p3)
      L=[tex('CJ15'),tex('CDF'),tex('D\O')]

    ax.set_xlabel(r'$y_W$',size=25)    
    ax.set_ylabel(r'$A_W$',size=25)

    ax.legend(H,L,frameon=0,loc=3,fontsize=22,numpoints=1,\
      bbox_to_anchor=(0.02, 0.65))
   
    ##ax.text(0.5,0.8,tex('nrep=%d'%nrows),transform=ax.transAxes,size=20)
    #py.tight_layout()
    py.tick_params(axis='both',labelsize=20)
    py.savefig('gallery/Wasy.pdf')
    py.close()

if __name__=='__main__':

  WASY()
