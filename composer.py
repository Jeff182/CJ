#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from tools import tex, fill_between
import lhapdf

class COMPOSER(object):

  def __init__(self,name,X=[None]):

    self.name=name
    self.SETS=lhapdf.mkPDFs(name)
    if X[0]==None: self.X=np.linspace(1e-3,0.99,100)
    else: self.X=X

  def map_X(self):
    f=lambda x: x**5
    X=self.X
    m=(X[0]-X[-1])/(f(X[0])-f(X[-1]))
    b=X[0]-m*f(X[0])
    self.X=m*f(X)+b

  def _get_xpdf(self,Set,flav,x,Q2):
    if   flav=='g': return Set.xfxQ2(21,x,Q2)
    elif flav=='u': return Set.xfxQ2(2,x,Q2)
    elif flav=='d': return Set.xfxQ2(1,x,Q2)
    elif flav=='s': return Set.xfxQ2(3,x,Q2)
    elif flav=='db+ub': return Set.xfxQ2(-2,x,Q2)+Set.xfxQ2(-1,x,Q2)
    elif flav=='db-ub': return Set.xfxQ2(-1,x,Q2)-Set.xfxQ2(-2,x,Q2)

  def _error(self,message):
    print 'ERR '+message
    sys.exit()

  def _get_symmetric_errors(self,PDFS):
    n=len(PDFS)-1
    feven=np.array([PDFS[2*i] for i in range(1,n/2)])
    fodd=np.array([PDFS[2*i-1] for i in range(1,n/2)])
    df=np.zeros(feven[0].size)
    for i in range(n/2-1):
      df+=(fodd[i]-feven[i])**2
    return df**0.5/2

  def _get_asymmetric_errors(self,PDFS):
    n=len(PDFS)-1
    f0=np.array(PDFS[0])
    feven=np.array([PDFS[2*i] for i in range(1,n/2)])
    fodd=np.array([PDFS[2*i-1] for i in range(1,n/2)])
    dfeven=feven-f0
    dfodd=fodd-f0
    zeros=np.zeros(f0.size)
    dfP=np.zeros(f0.size)
    dfM=np.zeros(f0.size)
    for i in range(n/2-1):
      dfP+=np.amax([dfodd[i],dfeven[i],zeros],0)**2
      dfM+=np.amax([-dfodd[i],-dfeven[i],zeros],0)**2
    return dfP**0.5,dfM**0.5

  def get_xpdf(self,flav=None,X=None,Q2=None):
    if flav==None: self._error('specify flav')
    if X==None: X=self.X
    if Q2==None: self._error('specify Q2')

    D={}
    PDFS=[[self._get_xpdf(Set,flav,x,Q2) for x in X] for Set in self.SETS]
    D['f0']=np.array(PDFS[0])
    D['df']=self._get_symmetric_errors(PDFS)
    D['df+'],D['df-']=self._get_asymmetric_errors(PDFS)
    return D



if __name__=="__main__" :

  CJ=COMPOSER(name='CJ12min')
  CJ.map_X() # optional
  print CJ.get_xpdf('g',Q2=10.0)


