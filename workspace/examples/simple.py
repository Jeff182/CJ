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
from master import FITPACK


# initialize PDF instances
HERA15=lhapdf.mkPDFs('HERAPDF15LO_EIG')
MMHT14=lhapdf.mkPDFs('MMHT2014lo68cl')
MSTW08=lhapdf.mkPDFs('MSTW2008lo68cl')

# dictionary for pdg flav idx
pdgmap= {'u':2,'d':1,'ub':-2,'db': -1,'s':3,'g':21}

# simple function to access lhad pdfs
get_pdfs = lambda flav,x,Q2,iset,grp: grp[iset].xfxQ2(pdgmap[flav],x,Q2)

# load CJ15 from fitpack data 
fname='data/CDBONN_KP_1_LO_6.pdf'
CJ15=FITPACK().get_PDFs(fname)

# dictionary for labels
labmap = {}
labmap['u']  = '$xu(x)$'
labmap['d']  = '$xd(x)$'
labmap['ub'] = r'$x\bar{u}(x)$'
labmap['db'] = r'$x\bar{d}(x)$'
labmap['s']  = '$xs(x)$'
labmap['g']  = '$xg(x)$'

# dictionary for grp == groups 
grpmap = {}
grpmap['HERA15'] = {'grp':HERA15,'color':'g-'}
grpmap['MMHT14'] = {'grp':MMHT14,'color':'b-'}
grpmap['MSTW08'] = {'grp':MSTW08,'color':'b:'}

# dictionary for ylims 
ymap={}
ymap['u'] ={'min':0.0,'max':0.8}
ymap['d'] ={'min':0.0,'max':0.6}
ymap['ub']={'min':0.0,'max':0.6}
ymap['db']={'min':0.0,'max':0.6}
ymap['s'] ={'min':0.0,'max':0.6}
ymap['g'] ={'min':0.0,'max':20.0}

# dictionary for plot location 
gs = gridspec.GridSpec(3,2) # specify plotting grid geometry
gs.update(left=0.1,right=0.98,wspace=0.3,hspace=0.1,top=0.98,bottom=0.1)
grid={}
grid['u'] = gs[0,0]
grid['d'] = gs[0,1]
grid['ub']= gs[1,0]
grid['db']= gs[1,1]
grid['s'] = gs[2,0]
grid['g'] = gs[2,1]


# setup kinematics
Q2=10.0
iset=0


# make plot
for flav in ['u','d','ub','db','s','g']:

  ax=py.subplot(grid[flav])

  for grp in ['CJ15','HERA15','MMHT14','MSTW08']:

    if grp=='CJ15':
      X=CJ15['Q2'][Q2]['X']
      central=CJ15['Q2'][Q2]['x'+flav]
      error=CJ15['Q2'][Q2]['err-x'+flav]*10

      #args={} 
      #args['ax']=ax
      #args['x']=X
      #args['central']=central
      #args['lower']=central-error
      #args['upper']=central+error
      #args['central color']='r'
      #args['central line style']='-'
      #args['band color']='#FFFF00'
      #args['label']=tex(grp)
      #CJ_Legend = plot_band(args)
      p1,=ax.plot(X,central,'r-')
      p2=fill_between(X,central-error,central+error,ax=ax,
        facecolor='#FFFF00',
        edgecolor='#FFFF00',
        alpha=1.0,hatch=None)

    else:
      X=np.linspace(1e-3,0.9,1000)
      grp_=grpmap[grp]['grp']
      col=grpmap[grp]['color']
      ax.plot(X,[get_pdfs(flav,x,Q2,iset,grp_) for x in X],col,label=tex(grp))

  # make legend
  if flav=='g':  
    H_,L_ = ax.get_legend_handles_labels()
    H=[(p2,p1)]
    L=[tex('CJ15')]
    for h in H_: H.append(h)
    for l in L_: L.append(l)
    ax.legend(H,L,loc=1,frameon=0,fontsize=15)

  # setup axis
  ax.semilogx()
  ax.set_xlim(1e-3,1.0)
  ax.set_ylim(ymap[flav]['min'],ymap[flav]['max'])
  if flav!='s' and flav!='g': ax.set_xticks([])
  ax.set_xlabel('$x$',size=20)
  ax.set_ylabel(labmap[flav],size=20)

  ## write info
  if flav=='u': 
    ax.text(0.1,0.1,'$Q^2=$'+tex('~10~GeV^2~LO'),
      transform=ax.transAxes,size=15)


py.savefig('plots/LOfits.pdf')


