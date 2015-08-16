#!/usr/bin/env python
import sys,os
import numpy as np
import pylab as py
from tools import tex#, fill_between
import lhapdf

class FITPACK(object):

  def get_PDFs(self,fname):
    
    # load file into list
    F=open(fname,'r')
    L=F.readlines()
    F.close()

    # clean data
    L=[l.strip() for l in L]
    L=[l for l in L if l.startswith('#')==False]
    L=[l for l in L if l.startswith('-')==False]
    L=[l for l in L if l!='']

    # get blocks markers
    I=[]
    for i in range(len(L)):
      l=L[i]
      if 'Q2' in l and '=' in l: I.append(i)
      if 'x' in l and '=' in l: I.append(i)
    I.append(len(L))

    # get blocks
    BLOCKS=[]
    for k in range(len(I[:-1])):
      block=[]
      for i in range(I[k],I[k+1]):
        block.append(L[i])
      BLOCKS.append(block)

    # format blocks
    D={'Q2':{},'x':{}}
    for block in BLOCKS:
 
      key,val = block[0].split('=')
      val=float(val)

      H_=block[1].split()
      H=[H_[0]]
      for h in H_[1:]: 
        H.append(h)
        H.append('err-'+h)

      
      data={}
      data_=block[2:]
      data_=np.transpose([[float(x) for x in l.split()] for l in data_])
      for i in range(len(H)):
        data[H[i]]=data_[i]

      D[key][val]=data

    return D

class COMPOSER(object):

  def __init__(self,name,X=[None],central_only=False):

    self.name=name
    self.central_only=central_only
    if central_only==True:
      self.central=lhapdf.mkPDF(name,0)
    else:
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
    elif flav=='ub': return Set.xfxQ2(-2,x,Q2)
    elif flav=='db': return Set.xfxQ2(-1,x,Q2)

  def _get_xpdf_central(self,flav,x,Q2):
    if   flav=='g': return self.central.xfxQ2(21,x,Q2)
    elif flav=='u': return self.central.xfxQ2(2,x,Q2)
    elif flav=='d': return self.central.xfxQ2(1,x,Q2)
    elif flav=='s': return self.central.xfxQ2(3,x,Q2)
    elif flav=='db+ub': return self.central.xfxQ2(-2,x,Q2)+self.central.xfxQ2(-1,x,Q2)
    elif flav=='db-ub': return self.central.xfxQ2(-1,x,Q2)-self.central.xfxQ2(-2,x,Q2)
    elif flav=='ub': return self.central.xfxQ2(-2,x,Q2)
    elif flav=='db': return self.central.xfxQ2(-1,x,Q2)

  def _error(self,message):
    print 'ERR '+message
    sys.exit()

  def _get_symmetric_errors(self,OBS):
    n=len(OBS)-1
    feven=np.array([OBS[2*i] for i in range(1,n/2)])
    fodd=np.array([OBS[2*i-1] for i in range(1,n/2)])
    df=np.zeros(feven[0].size)
    for i in range(n/2-1):
      df+=(fodd[i]-feven[i])**2
    return df**0.5/2

  def _get_asymmetric_errors(self,OBS):
    n=len(OBS)-1
    f0=np.array(OBS[0])
    feven=np.array([OBS[2*i] for i in range(1,n/2)])
    fodd=np.array([OBS[2*i-1] for i in range(1,n/2)])
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
    if self.central_only:
      D['xf0']=np.array([self._get_xpdf_central(flav,x,Q2) for x in X])
      D['dxf']=np.zeros(X.size)
      D['dxf+']=np.zeros(X.size)
      D['dxf-']=np.zeros(X.size)
    else:
      PDFS=[[self._get_xpdf(Set,flav,x,Q2) for x in X] for Set in self.SETS]
      D['xf0']=np.array(PDFS[0])
      D['dxf']=self._get_symmetric_errors(PDFS)
      D['dxf+'],D['dxf-']=self._get_asymmetric_errors(PDFS)
    return D

  def get_dou(self,X=None,Q2=None):
    if X==None: X=self.X
    if Q2==None: self._error('specify Q2')
    D={}
    d=lambda Set,x,Q2: self._get_xpdf(Set,'d',x,Q2) 
    u=lambda Set,x,Q2: self._get_xpdf(Set,'u',x,Q2) 
    OBS=[[d(Set,x,Q2)/u(Set,x,Q2) for x in X] for Set in self.SETS]
    D['central']=np.array(OBS[0])
    D['sym err']=self._get_symmetric_errors(OBS)
    D['asym err +'],D['asym err -']=self._get_asymmetric_errors(OBS)
    return D

class COMPOSER4NNPDF(object):

  def __init__(self,name):

    self.name=name
    self.SETS=lhapdf.mkPDFs(name)

  def _get_xpdf(self,Set,flav,x,Q2):
    if   flav=='g': return Set.xfxQ2(21,x,Q2)
    elif flav=='u': return Set.xfxQ2(2,x,Q2)
    elif flav=='d': return Set.xfxQ2(1,x,Q2)
    elif flav=='s': return Set.xfxQ2(3,x,Q2)
    elif flav=='db+ub': return Set.xfxQ2(-2,x,Q2)+Set.xfxQ2(-1,x,Q2)
    elif flav=='db-ub': return Set.xfxQ2(-1,x,Q2)-Set.xfxQ2(-2,x,Q2)
    elif flav=='ub': return Set.xfxQ2(-2,x,Q2)
    elif flav=='db': return Set.xfxQ2(-1,x,Q2)

  def _error(self,message):
    print 'ERR '+message
    sys.exit()

  def get_xpdf(self,flav=None,X=None,Q2=None):
    if flav==None: self._error('specify flav')
    if X==None: X=self.X
    if Q2==None: self._error('specify Q2')
    D={}
    PDFS=[[self._get_xpdf(Set,flav,x,Q2) for x in X] for Set in self.SETS]
    D['xf0']=np.mean(PDFS,axis=0)
    D['dxf+']=np.var(PDFS,axis=0)**0.5
    D['dxf-']=np.var(PDFS,axis=0)**0.5
    return D

  def get_dou(self,X=None,Q2=None):
    if X==None: X=self.X
    if Q2==None: self._error('specify Q2')
    D={}
    d=lambda Set,x,Q2: self._get_xpdf(Set,'d',x,Q2) 
    u=lambda Set,x,Q2: self._get_xpdf(Set,'u',x,Q2) 
    OBS=[[d(Set,x,Q2)/u(Set,x,Q2) for x in X] for Set in self.SETS]
    D['central']=np.mean(OBS,axis=0)
    D['asym err +']=np.var(OBS,axis=0)**0.5
    D['asym err -']=np.var(OBS,axis=0)**0.5
    return D


if __name__=="__main__" :

  CJ=COMPOSER(name='CJ12min')
  CJ.map_X() # optional
  print CJ.get_xpdf('g',Q2=10.0)








