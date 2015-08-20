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
#params={'text.latex.preamble': \
#  [r"\usepackage{upgreek}",
#   r"\usepackage[nice]{units}",
#   r"\usepackage{amstext}"
#  ]}
#py.rcParams.update(params)

from master import COMPOSER,FITPACK,COMPOSER4NNPDF

# plotting functions with LR capabilities

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

def plotI(AX,Q2,central,other,flav,color,ls,T=1,hatch=None,alpha=0.4,facecolor='none',edgecolor='none'):
  """
  This routine plots a ratio  band of 'other' to 'central'  using the error from 'other'
  """

  axL,axR=AX[flav]
  XL=10**np.linspace(-4,-1,100)
  if flav=='u' or flav=='d' or flav=='g':
    XR=np.linspace(0.1,0.95,100)
  else:
    XR=np.linspace(0.1,0.5,100)
  LC=central.get_xpdf(flav,X=XL,Q2=Q2)
  RC=central.get_xpdf(flav,X=XR,Q2=Q2)
  LO=other.get_xpdf(flav,X=XL,Q2=Q2)
  RO=other.get_xpdf(flav,X=XR,Q2=Q2)

  axL.plot(XL,LO['xf0']/LC['xf0'],color=color,ls=ls)
  p1,=axR.plot(XR,RO['xf0']/RC['xf0'],color=color,ls=ls,label=flav)

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

  axL.set_xlim(np.amin(XL),0.1)
  axR.set_xlim(0.1,np.amax(XR))
  return (p2,p1)

def plotII(AX,Q2,central,other,flav,color,ls,T=1,hatch=None,alpha=0.4,facecolor='none',edgecolor='none'):
  """
  This routine plots a ratio (only central) of 'other' to 'central'  using the error from 'other'
  """

  axL,axR=AX[flav]
  XL=10**np.linspace(-4,-1,100)
  if flav=='u' or flav=='d' or flav=='g':
    XR=np.linspace(0.1,0.95,100)
  else:
    XR=np.linspace(0.1,0.5,100)
  LC=central.get_xpdf(flav,X=XL,Q2=Q2)
  RC=central.get_xpdf(flav,X=XR,Q2=Q2)
  LO=other.get_xpdf(flav,X=XL,Q2=Q2)
  RO=other.get_xpdf(flav,X=XR,Q2=Q2)

  axL.plot(XL,LO['xf0']/LC['xf0'],color=color,ls=ls)
  p1,=axR.plot(XR,RO['xf0']/RC['xf0'],color=color,ls=ls,label=flav)

  axL.set_xlim(np.amin(XL),0.1)
  axR.set_xlim(0.1,np.amax(XR))
  return p1

# ploting function without LR 

def plotI0(AX,Q2,central,other,flav,color,ls,T=1,hatch=None,alpha=0.4,facecolor='none',edgecolor='none'):
  """
  This routine plots a ratio  band of 'other' to 'central'  
  using the error from 'other'
  """

  ax=AX[flav]
  if flav=='u' or flav=='d' or flav=='g':
    X=np.linspace(1e-4,0.95,100)
  else:
    X=np.linspace(1e-4,0.5,100)
  C=central.get_xpdf(flav,X=X,Q2=Q2)
  O=other.get_xpdf(flav,X=X,Q2=Q2)

  ax.plot(X,O['xf0']/C['xf0'],color=color,ls=ls)
  p1,=ax.plot(X,O['xf0']/C['xf0'],color=color,ls=ls,label=flav)

  p2=fill_between(X,
    (O['xf0']-O['dxf-']*T)/C['xf0'],
    (O['xf0']+O['dxf+']*T)/C['xf0'],
    ax=ax,
    facecolor=facecolor,
    edgecolor=edgecolor,
    alpha=alpha,
    hatch=hatch)
  return (p2,p1)

def plotII0(AX,Q2,central,other,flav,color,ls,T=1,hatch=None,alpha=0.4,facecolor='none',edgecolor='none'):
  """
  This routine plots a ratio (only central) of 
  'other' to 'central'  using the error from 'other'
  """

  ax=AX[flav]
  if flav=='u' or flav=='d' or flav=='g':
    X=np.linspace(1e-4,0.95,100)
  else:
    X=np.linspace(1e-4,0.5,100)
  C=central.get_xpdf(flav,X=X,Q2=Q2)
  O=other.get_xpdf(flav,X=X,Q2=Q2)

  ax.plot(X,O['xf0']/C['xf0'],color=color,ls=ls)
  p1,=ax.plot(X,O['xf0']/C['xf0'],color=color,ls=ls,label=flav)
  return p1

# main routines

def ratio():

  ###############################
  # plot geometry
  ###############################
  # set gloabl figure dimensions
  nrows=3
  ncols=2
  py.figure(figsize=(ncols*4,nrows*3))

  # construct LR AXs for each flav 
  AX={}

  # left side LR panels
  gs = gridspec.GridSpec(nrows,ncols)
  gs.update(left=0.1,right=0.48,wspace=0,hspace=0.3,\
    top=0.98,bottom=0.05)
  AX['u'] =get_LR(gs[0,0],gs[0,1])
  AX['ub']=get_LR(gs[1,0],gs[1,1])
  AX['s'] =get_LR(gs[2,0],gs[2,1])

  # right side LR panels
  gs = gridspec.GridSpec(nrows,ncols)
  gs.update(left=0.58,right=0.98,wspace=0,hspace=0.3,\
    top=0.98,bottom=0.05)
  AX['d'] =get_LR(gs[0,0],gs[0,1])
  AX['db']=get_LR(gs[1,0],gs[1,1])
  AX['g'] =get_LR(gs[2,0],gs[2,1])

  
  ###############################
  # plot content 
  ###############################
  # initialize composer for each pdf set
  CJ=COMPOSER('CJ15_NLO_KP_AV18')
  HERA15=COMPOSER('HERAPDF15NLO_EIG')
  MMHT14=COMPOSER('MMHT2014nlo68cl')
  NNPDF=COMPOSER4NNPDF('NNPDF30_nlo_as_0118')
  #CT10=COMPOSER('CT10nlo')
  
  Q2=10
  for flav in ['u','d','ub','db','s','g']:
    print 'plotting '+flav

    # the output of plot used to construct specialized 
    # legend marker
    p1=plotI(AX,Q2,CJ,CJ,flav,'yellow','-',T=10,\
      hatch=None,alpha=0.4,facecolor='yellow',edgecolor='yellow')
    p2=plotI(AX,Q2,CJ,CJ,flav,'r','-',T=1,\
      hatch=None,alpha=1.0,facecolor='r',edgecolor='none')
    p3=plotI(AX,Q2,CJ,MMHT14,flav,'b','-',T=1,\
      hatch='...',alpha=0.4,facecolor='none',edgecolor='b')
    p4=plotI(AX,Q2,CJ,HERA15,flav,'g','-',T=1,\
      hatch='...',alpha=0.4,facecolor='none',edgecolor='g')
    p5=plotI(AX,Q2,CJ,NNPDF,flav,'m','-',T=1,\
      hatch='...',alpha=0.4,facecolor='none',edgecolor='m')
    #p5=plotI(AX,Q2,CJ,CT10,flav,'k','-',T=1,\
    #  hatch='\\',alpha=0.4,facecolor='none',edgecolor='k')
  

# print single PDF as a function of x (grafted onto ratio.py code)
    if flav=='db':
      X=np.linspace(0.01,1.0,100)    
      mm=MMHT14.get_xpdf(flav,X=X,Q2=4.0)
      t=1
      for i in range(X.size):
        l='%0.2f %0.2e %0.2e %0.2e'       
        print l%(X[i],X[i]*mm['xf0'][i],X[i]*mm['dxf-'][i]*t,X[i]*mm['dxf+'][i]*t)
        #l='%0.2f %0.2e %0.2e'
        #print l%(X[i],X[i]*mm['xf0'][i],X[i]*mm['dxf'][i]*t)
      ax=py.subplot(111)
      ax.plot(X,X*mm['xf0'],'r-')
      ax.plot(X,X*mm['xf0']-X*mm['dxf-']*t,'r:')
      ax.plot(X,X*mm['xf0']+X*mm['dxf+']*t,'r:')
      #ax.plot(X,X*mm['xf0']-mm['dxf']*t,'r:')
      #ax.plot(X,X*mm['xf0']+mm['dxf']*t,'r:') 
      py.show()
      sys.exit()



    #retrieve LR ax for further proccesing 
    axL,axR=AX[flav]

    # plot specialized legend marker at the specific flav panel
    if flav=='d':
      H=[p1,p2,p3,p4,p5]
      L=[tex('CJ15')+'\ ($T=10$)'\
        ,tex('CJ15')+'\ ($T=1$)'\
        ,tex('MMHT14')\
        ,tex('HERA15')\
        ,tex('NNPDF3.0')]
      #  ,tex('CT10')]
      axR.legend(H,L,frameon=0,fontsize=11,\
        bbox_to_anchor=(0.05, 1.0))
    if flav=='u':
      axR.text(0.1,0.85,'$Q^2=%0.0f$'%(Q2)+tex('~GeV^2'),\
        transform=axL.transAxes,size=15)

    # plot some text at specific flav panel
    #if flav=='d':
    #  axL.text(0.06,0.1,'$T=10$',transform=axL.transAxes,size=20)

    # set ylims
    if flav=='u':
      ymin,ymax=0.795,1.205
    if flav=='d':
      ymin,ymax=0.6,1.7
    if flav=='s':
      ymin,ymax=0.6,1.9
    if flav=='g' or flav=='ub' or flav=='db':
      ymin,ymax=0.7,1.3
    #axL.set_yticks(np.arange(0.2,2.0,0.4))
    axL.set_ylim(ymin,ymax)
    axR.set_ylim(ymin,ymax)

    # set x axis
    axL.semilogx()
    axL.set_xticks(10**np.linspace(-4,-2,3))
    if flav=='u' or flav=='d':
      axR.set_xticks(np.arange(0.1,1.1,0.2))
      axR.set_xlim(0.1,0.95)
    else:
      axR.set_xticks(np.arange(0.1,0.4,0.1))
      axR.set_xlim(0.1,0.4)

    # set labels
    if flav=='db': _flav=r'\bar{d}'
    elif flav=='ub': _flav=r'\bar{u}'
    else: _flav=flav
    axL.set_ylabel(r'$%s/%s$\large$_{\rm CJ15}$'%(_flav,_flav),size=20)
    axL.set_xlabel('$x$',size=20)
    axL.xaxis.set_label_coords(1.0,-0.08,transform=axL.transAxes)

  py.savefig('gallery/ratio.pdf')

def ratio_wfn():

  ###############################
  # plot geometry
  ###############################
  # set gloabl figure dimensions
  nrows=3
  ncols=2
  py.figure(figsize=(ncols*4,nrows*3))

  # construct LR AXs for each flav 
  AX={}
  gs = gridspec.GridSpec(nrows,ncols)
  gs.update(left=0.1,right=0.98,wspace=0.3,hspace=0.3,top=0.98,bottom=0.08)
  AX['u'] =py.subplot(gs[0,0])  
  AX['ub']=py.subplot(gs[1,0])
  AX['s'] =py.subplot(gs[2,0])
  AX['d'] =py.subplot(gs[0,1])
  AX['db']=py.subplot(gs[1,1])
  AX['g'] =py.subplot(gs[2,1])


  ###############################
  # plot content 
  ###############################

  CJ15={'KP':{},'fmKP':{}}
  CJ15['KP']['AV18' ]   =COMPOSER('CJ15_NLO_KP_AV18')
  CJ15['KP']['CDBONN']  =COMPOSER('CJ15_NLO_KP_CDBONN',central_only=True)
  CJ15['KP']['WJC1' ]   =COMPOSER('CJ15_NLO_KP_WJC1',central_only=True)
  CJ15['KP']['WJC2' ]   =COMPOSER('CJ15_NLO_KP_WJC2',central_only=True)
  CJ15['fmKP']['AV18']  =COMPOSER('CJ15_NLO_fmKP_AV18',central_only=True)
  CJ15['fmKP']['CDBONN']=COMPOSER('CJ15_NLO_fmKP_CDBONN',central_only=True)
  CJ15['fmKP']['WJC1']  =COMPOSER('CJ15_NLO_fmKP_WJC1',central_only=True)
  CJ15['fmKP']['WJC2']  =COMPOSER('CJ15_NLO_fmKP_WJC2',central_only=True)
  
  Q2=10
  for flav in ['u','d','ub','db','s','g']:
    print flav
    p1=plotI0(AX,Q2,CJ15['KP']['AV18'],CJ15['KP']['AV18'],flav,'yellow','-',
      T=1,hatch=None,alpha=0.4,facecolor='yellow',edgecolor='none')
    p2=plotII0(AX,Q2,CJ15['KP']['AV18'],CJ15['KP']['CDBONN'],flav,'r','-',
      T=1,hatch=None,alpha=0.4,facecolor='b',edgecolor='none')
    p3=plotII0(AX,Q2,CJ15['KP']['AV18'],CJ15['KP']['WJC1'],flav,'g','--',
      T=1,hatch=None,alpha=0.4,facecolor='g',edgecolor='none')
    p4=plotII0(AX,Q2,CJ15['KP']['AV18'],CJ15['KP']['WJC2'],flav,'b','-',
      T=1,hatch=None,alpha=0.4,facecolor='y',edgecolor='none')
  
    ax=AX[flav]

    if flav=='u':
      H=[p1,p2,p3,p4]
      L=[tex('CJ15\ (AV18)')+r'$\ \ \ T=1$',tex('CDBonn'),tex('WJC1'),tex('WJC2')]
      ax.legend(H,L,frameon=0,fontsize=10, bbox_to_anchor=(0.65, 1.0))
      ax.text(0.06,0.1,'$Q^2=%0.0f$'%(Q2)+tex('~GeV^2'),\
        transform=ax.transAxes,size=15)

    # set ylims
    if flav=='d':
      ymin,ymax=0.85,1.15
    elif flav=='u':
      ymin,ymax=0.95,1.05
    elif flav=='g':
      ymin,ymax=0.9,1.1
    else:
      ymin,ymax=0.95,1.05
    ax.set_ylim(ymin,ymax)

    # set x axis
    if flav=='u' or flav=='d' or flav=='g':
      #ax.set_xticks(np.arange(0.1,1.1,0.2))
      ax.set_xlim(0,0.95)
    elif flav=='g':
      ax.set_xlim(0,0.8)
    else:
      #ax.set_xticks(np.arange(0.1,0.5,0.1))
      ax.set_xlim(0,0.5)

    # set labels
    if flav=='db': _flav=r'\bar{d}'
    elif flav=='ub': _flav=r'\bar{u}'
    else: _flav=flav
    ax.set_ylabel(r'$%s/%s$\large$_{\rm CJ15}$'%(_flav,_flav),size=20)
    ax.set_xlabel('$x$',size=20)
    if flav=='u' or flav=='d' or flav=='g':
      xticklabels=[0,0.2,0.4,0.6,0.8]
  
  py.savefig('gallery/ratio_wfn.pdf')

def ratio_off():

  ###############################
  # plot geometry
  ###############################
  # set gloabl figure dimensions
  nrows=3
  ncols=2
  py.figure(figsize=(ncols*4,nrows*3))

  # construct LR AXs for each flav 
  AX={}
  gs = gridspec.GridSpec(nrows,ncols)
  gs.update(left=0.1,right=0.98,wspace=0.3,hspace=0.3,top=0.98,bottom=0.08)
  AX['u'] =py.subplot(gs[0,0])  
  AX['ub']=py.subplot(gs[1,0])
  AX['s'] =py.subplot(gs[2,0])
  AX['d'] =py.subplot(gs[0,1])
  AX['db']=py.subplot(gs[1,1])
  AX['g'] =py.subplot(gs[2,1])


  ###############################
  # plot content 
  ###############################
  
  CJ15={'KP':{},'fmKP':{}}
  CJ15['KP']['AV18' ]   =COMPOSER('CJ15_NLO_KP_AV18')
  CJ15['KP']['CDBONN']  =COMPOSER('CJ15_NLO_KP_CDBONN',central_only=True)
  CJ15['KP']['WJC1' ]   =COMPOSER('CJ15_NLO_KP_WJC1',central_only=True)
  CJ15['KP']['WJC2' ]   =COMPOSER('CJ15_NLO_KP_WJC2',central_only=True)
  CJ15['fmKP']['AV18']  =COMPOSER('CJ15_NLO_fmKP_AV18',central_only=True)
  CJ15['fmKP']['CDBONN']=COMPOSER('CJ15_NLO_fmKP_CDBONN',central_only=True)
  CJ15['fmKP']['WJC1']  =COMPOSER('CJ15_NLO_fmKP_WJC1',central_only=True)
  CJ15['fmKP']['WJC2']  =COMPOSER('CJ15_NLO_fmKP_WJC2',central_only=True)
  
  Q2=10
  for flav in [ 'u','d','ub','db','s','g']:
    print flav
    p1=plotI0(AX,Q2,CJ15['KP']['AV18'],CJ15['KP']['AV18'],flav,'y','-',
      T=1,hatch=None,alpha=0.4,facecolor='yellow',edgecolor='none')
    p2=plotII0(AX,Q2,CJ15['KP']['AV18'],CJ15['fmKP']['AV18'],flav,'k',':',
      T=1,hatch=None,alpha=0.4,facecolor='k',edgecolor='none')
    p3=plotII0(AX,Q2,CJ15['KP']['AV18'],CJ15['fmKP']['CDBONN'],flav,'r','-',
      T=1,hatch=None,alpha=0.4,facecolor='b',edgecolor='none')
    p4=plotII0(AX,Q2,CJ15['KP']['AV18'],CJ15['fmKP']['WJC1'],flav,'g','--',
      T=1,hatch=None,alpha=0.4,facecolor='g',edgecolor='none')
    p5=plotII0(AX,Q2,CJ15['KP']['AV18'],CJ15['fmKP']['WJC2'],flav,'b','-',
      T=1,hatch=None,alpha=0.4,facecolor='y',edgecolor='none')
  
    ax=AX[flav]

    if flav=='d':
      H=[p1,p2,p3,p4,p5]
      L=[tex('CJ15\ (AV18)'),tex('AV18+OCS'),tex('CDBonn+OCS'),tex('WJC1+OCS'),tex('WJC2+OCS')]
      ax.legend(H,L,frameon=0,fontsize=10, bbox_to_anchor=(0.55, 0.55))
    if flav=='u':
      ax.text(0.06,0.85,'$Q^2=%0.0f$'%(Q2)+tex('~GeV^2'),\
        transform=ax.transAxes,size=15)

    # set ylims
    if flav=='d':
      ymin,ymax=0.8,1.1
    elif flav=='u':
      ymin,ymax=0.95,1.05
    elif flav=='g':
      ymin,ymax=0.9,1.1
    else:
      ymin,ymax=0.95,1.05
    ax.set_ylim(ymin,ymax)

    # set x axis
    if flav=='u' or flav=='d' or flav=='g':
      #ax.set_xticks(np.arange(0.1,1.1,0.2))
      ax.set_xlim(0,0.95)
    elif flav=='g':
      ax.set_xlim(0,0.8)
    else:
      #ax.set_xticks(np.arange(0.1,0.5,0.1))
      ax.set_xlim(0,0.5)

    # set labels
    if flav=='db': _flav=r'\bar{d}'
    elif flav=='ub': _flav=r'\bar{u}'
    else: _flav=flav
    ax.set_ylabel(r'$%s/%s$\large$_{\rm CJ15}$'%(_flav,_flav),size=20)
    ax.set_xlabel('$x$',size=20)
    if flav=='u' or flav=='d' or flav=='g':
      xticklabels=[0,0.2,0.4,0.6,0.8]
  
    # set labels
    if flav=='db': _flav=r'\bar{d}'
    elif flav=='ub': _flav=r'\bar{u}'
    else: _flav=flav
    ax.set_ylabel(r'$%s/%s$\large$_{\rm CJ15}$'%(_flav,_flav),size=20)
    ax.set_xlabel('$x$',size=20)
  
  py.savefig('gallery/ratio_off.pdf')


if __name__=='__main__':
  ratio()
  #ratio_wfn()
  #ratio_off()
