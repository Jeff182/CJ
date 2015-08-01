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
from master import COMPOSER,FITPACK

# set fig size
nrows=3
ncols=2
py.figure(figsize=(ncols*4,nrows*3))

# prepare grid 
def get_LR(gsL,gsR):
  # construct (L)eft panel
  axL=py.subplot(gsL)
  axL.spines['right'].set_visible(False)
  axL.yaxis.set_ticks_position('left')
  axL.axvline(0.1,color='k',ls='--',alpha=0.5)
  
  # construct (R)ight panel
  axR=py.subplot(gsR)
  axR.spines['left'].set_visible(False)
  axR.axes.yaxis.set_ticklabels([])
  axR.axes.get_yaxis().set_ticks([])
  return [axL,axR]

AX={}

gs = gridspec.GridSpec(3,2)
gs.update(left=0.1,right=0.48,wspace=0,hspace=0.3,top=0.98,bottom=0.05)
AX['u']=get_LR(gs[0,0],gs[0,1])
AX['ub']=get_LR(gs[1,0],gs[1,1])
AX['s']=get_LR(gs[2,0],gs[2,1])

gs = gridspec.GridSpec(3,2)
gs.update(left=0.58,right=0.98,wspace=0,hspace=0.3,top=0.98,bottom=0.05)
AX['d']=get_LR(gs[0,0],gs[0,1])
AX['db']=get_LR(gs[1,0],gs[1,1])
AX['g']=get_LR(gs[2,0],gs[2,1])

def plot(central,other,flav,color,T=1,hatch=None,alpha=0.4,facecolor='none',edgecolor='none'):

  axL,axR=AX[flav]
  XL=10**np.linspace(-4,-1,100)
  if flav=='u' or flav=='d':
    XR=np.linspace(0.1,1,100)
  else:
    XR=np.linspace(0.1,0.9,100)
  LC=central.get_xpdf(flav,X=XL,Q2=Q2)
  RC=central.get_xpdf(flav,X=XR,Q2=Q2)
  LO=other.get_xpdf(flav,X=XL,Q2=Q2)
  RO=other.get_xpdf(flav,X=XR,Q2=Q2)

  axL.plot(XL,LO['xf0']/LC['xf0'],color=color,ls='-')
  p1,=axR.plot(XR,RO['xf0']/RC['xf0'],color=color,ls='-',label=flav)

  fill_between(XL,
    (LO['xf0']-LO['dxf-']*T)/LC['xf0'],
    (LO['xf0']+LO['dxf+']*T)/LC['xf0'],
    ax=axL,
    facecolor=facecolor,
    edgecolor=edgecolor,
    alpha=alpha,
    hatch=hatch)
  p2=fill_between(XR,
    (RO['xf0']-RO['dxf-']*T)/RC['xf0'],
    (RO['xf0']+RO['dxf+']*T)/RC['xf0'],
    ax=axR,
    facecolor=facecolor,
    edgecolor=edgecolor,
    alpha=alpha,
    hatch=hatch)

  # y axis
  if flav=='u' or flav=='d':
    ymin,ymax=0.7,1.3
  else: 
    ymin,ymax=0.7,1.3
    axL.set_yticks(np.arange(0.2,2.0,0.4))
  axL.set_ylim(ymin,ymax)
  axR.set_ylim(ymin,ymax)
  axL.set_xlim(np.amin(XL),0.1)
  axR.set_xlim(0.1,np.amax(XR))

  # x axis
  axL.semilogx()
  axL.set_xticks(10**np.linspace(-4,-2,3))
  if flav=='u' or flav=='d':
    axR.set_xticks(np.arange(0.1,1.1,0.2))
    axR.set_xlim(0.1,1.0)
  else:
    axR.set_xticks(np.arange(0.1,0.4,0.1))
    axR.set_xlim(0.1,0.4)

  # set labels
  if flav=='db': _flav=r'\bar{d}'
  elif flav=='ub': _flav=r'\bar{u}'
  else: _flav=flav
  axL.set_ylabel(r'$%s/%s_{\rm CJ15}$'%(_flav,_flav),size=20)
  axL.set_xlabel('$x$',size=20)
  axL.xaxis.set_label_coords(1.0,-0.08,transform=axL.transAxes)

  return (p2,p1)

def plot2(central,other,flav,color,T=1,hatch=None,alpha=0.4,facecolor='none',edgecolor='none'):

  axL,axR=AX[flav]
  XL=10**np.linspace(-4,-1,100)
  if flav=='u' or flav=='d':
    XR=np.linspace(0.1,1.0,100)
  else:
    XR=np.linspace(0.1,0.9,100)
  LC=central.get_xpdf(flav,X=XL,Q2=Q2)
  RC=central.get_xpdf(flav,X=XR,Q2=Q2)
  LO=np.array([other.get_xpdf_central(flav,_XL,Q2) for _XL in XL])
  RO=np.array([other.get_xpdf_central(flav,_XR,Q2) for _XR in XR])

  axL.plot(XL,LO/LC['xf0'],color=color,ls='-')
  p1,=axR.plot(XR,RO/RC['xf0'],color=color,ls='-',label=flav)

  # y axis
  if flav=='u' or flav=='d':
    ymin,ymax=0.7,1.3
  else: 
    ymin,ymax=0.7,1.3
    axL.set_yticks(np.arange(0.2,2.0,0.4))
  axL.set_ylim(ymin,ymax)
  axR.set_ylim(ymin,ymax)
  axL.set_xlim(np.amin(XL),0.1)
  axR.set_xlim(0.1,np.amax(XR))

  # x axis
  axL.semilogx()
  axL.set_xticks(10**np.linspace(-4,-2,3))
  if flav=='u' or flav=='d':
    axR.set_xticks(np.arange(0.1,1.1,0.2))
    axR.set_xlim(0.1,1.0)
  else:
    axR.set_xticks(np.arange(0.1,0.4,0.1))
    axR.set_xlim(0.1,1.0)

  # set labels
  if flav=='db': _flav=r'\bar{d}'
  elif flav=='ub': _flav=r'\bar{u}'
  else: _flav=flav
  axL.set_ylabel(r'$%s/%s_{\rm CJ15}$'%(_flav,_flav),size=20)
  axL.set_xlabel('$x$',size=20)
  axL.xaxis.set_label_coords(1.0,-0.08,transform=axL.transAxes)

  return p1

CJ12={}
CJ12['min']=COMPOSER('CJ12min')
CJ12['mid']=COMPOSER('CJ12mid')
CJ15={'KP':{},'fmKP':{}}
CJ15['KP']['AV18' ]   =COMPOSER('CJ15_NLO_KP_AV18')
CJ15['KP']['CDBONN']  =COMPOSER('CJ15_NLO_KP_CDBONN',central_only=True)
CJ15['KP']['WJC1' ]   =COMPOSER('CJ15_NLO_KP_WJC1',central_only=True)
CJ15['KP']['WJC2' ]   =COMPOSER('CJ15_NLO_KP_WJC2',central_only=True)
CJ15['fmKP']['AV18']  =COMPOSER('CJ15_NLO_fmKP_AV18',central_only=True)
CJ15['fmKP']['CDBONN']=COMPOSER('CJ15_NLO_fmKP_CDBONN',central_only=True)
CJ15['fmKP']['WJC1']  =COMPOSER('CJ15_NLO_fmKP_WJC1',central_only=True)
CJ15['fmKP']['WJC2']  =COMPOSER('CJ15_NLO_fmKP_WJC2',central_only=True)

#print CJ15['KP']['CDBONN'].get_xpdf_central('u',0.5,100)
#sys.exit()
#'MMHT2014lo68cl'
#'MMHT2014nlo68cl'
#'MSTW2008lo68cl'
#'JR14NLO08VF'
#'MSTW2008nlo68cl'
#'MSTW2008nlo90cl'
#'HERAPDF15LO_EIG'
#'HERAPDF15NLO_EIG'

Q2=100
for flav in ['u','d','ub','db','s','g']:
  print flav
  p1=plot(CJ15['KP']['AV18'],CJ15['KP']['AV18'],flav,'r',
    T=1,hatch=None,alpha=0.4,facecolor='r',edgecolor='none')
  p2=plot2(CJ15['KP']['AV18'],CJ15['fmKP']['AV18'],flav,'k',
    T=1,hatch=None,alpha=0.4,facecolor='k',edgecolor='none')
  p3=plot2(CJ15['KP']['AV18'],CJ15['fmKP']['CDBONN'],flav,'b',
    T=1,hatch=None,alpha=0.4,facecolor='b',edgecolor='none')
  p4=plot2(CJ15['KP']['AV18'],CJ15['fmKP']['WJC1'],flav,'g',
    T=1,hatch=None,alpha=0.4,facecolor='g',edgecolor='none')
  p5=plot2(CJ15['KP']['AV18'],CJ15['fmKP']['WJC2'],flav,'m',
    T=1,hatch=None,alpha=0.4,facecolor='m',edgecolor='none')

  if flav=='u':
    axL,axR=AX[flav]
    H=[p1,p2,p3,p4,p5]
    L=[tex('CJ15'),tex('fmKP\ AV18'),tex('fmKP\ CDBONN'),tex('fmKP\ WJC1'),tex('fmKP\ WJC2')]
    axR.legend(H,L,frameon=0,fontsize=15, bbox_to_anchor=(0.1, 1.0))
    axR.text(0.06,0.1,'$Q^2=%0.0f$'%(Q2)+tex('~GeV^2'),transform=axL.transAxes,size=20)
  #elif flav=='d':
    #axL.text(0.06,0.1,'$T=5$',transform=axL.transAxes,size=20)

py.savefig('gallery/ratio_off.pdf')
