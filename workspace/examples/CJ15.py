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
CJ12min=lhapdf.mkPDFs('CJ12min')

# dictionary for PDG flav index
pdgmap= {'u':2,'d':1,'ub':-2,'db': -1,'s':3,'g':21}

# load CJ15 from fitpack data 
CJ15 = {}
CJ15['fmKP_AV18']   = FITPACK().get_PDFs('data/fmKP_AV18.pdf')
CJ15['fmKP_CDBONN'] = FITPACK().get_PDFs('data/fmKP_CDBONN.pdf')
CJ15['fmKP_WJC1']   = FITPACK().get_PDFs('data/fmKP_WJC1.pdf')
CJ15['fmKP_WJC2']   = FITPACK().get_PDFs('data/fmKP_WJC2.pdf')
CJ15['tst_AV18_KP']   = FITPACK().get_PDFs('data/tst_AV18_KP.pdf')
CJ15['tst_CDBONN_KP'] = FITPACK().get_PDFs('data/tst_CDBONN_KP.pdf')
CJ15['tst_WJC1_KP']   = FITPACK().get_PDFs('data/tst_WJC1_KP.pdf')
CJ15['tst_WJC2_KP']   = FITPACK().get_PDFs('data/tst_WJC2_KP.pdf')

# dictionary for labels
labmap = {}
labmap['u']  = '$xu$'
labmap['d']  = '$xd$'
labmap['ub'] = r'$x\bar{u}$'
labmap['db'] = r'$x\bar{d}$'
labmap['s']  = '$xs$'
labmap['g']  = '$xg$'
labmap['du']   = '$d/u$'
labmap['dbub'] = r'$\bar{d}/\bar{u}$'

# dictionary for grp == groups 
grpmap = {}
grpmap['fmKP_AV18']     = {'color':'r-'}
grpmap['fmKP_CDBONN']   = {'color':'g-'}
grpmap['fmKP_WJC1']     = {'color':'b-'}
grpmap['fmKP_WJC2']     = {'color':'k:'}
grpmap['tst_AV18_KP']   = {'color':'r-.'}
grpmap['tst_CDBONN_KP'] = {'color':'k:'}
grpmap['tst_WJC1_KP']   = {'color':'r-.'}
grpmap['tst_WJC2_KP']   = {'color':'k:'}

# dictionary for ylims 
ymap = {}
ymap['u'] ={'min':0.0,'max':0.75}
ymap['d'] ={'min':0.0,'max':0.6}
ymap['du']  ={'min':0.0,'max':1.0}
ymap['dbub']={'min':0.6,'max':1.8}
ymap['ub']={'min':0.0,'max':0.6}
ymap['db']={'min':0.0,'max':0.6}
ymap['s'] ={'min':0.0,'max':0.45}
ymap['g'] ={'min':0.0,'max':10.0}

# dictionary for plot location 
gs = gridspec.GridSpec(4,2) # specify plotting grid geometry
gs.update(left=0.1,right=0.98,wspace=0.2,hspace=0.4,top=0.98,bottom=0.1)
grid = {}
grid['u'] = gs[0,0]
grid['d'] = gs[0,1]
grid['du']  = gs[1,0]
grid['dbub']= gs[1,1]
grid['ub']= gs[2,0]
grid['db']= gs[2,1]
grid['s'] = gs[3,0]
grid['g'] = gs[3,1]


# setup kinematics
Q2=10.0
iset=0


# make plot
for flav in ['u','d','ub','db','s','g','du','dbub']:

  ax=py.subplot(grid[flav])

  for grp in ['fmKP_AV18','fmKP_CDBONN','fmKP_WJC1','fmKP_WJC2']:

    cj15 = CJ15[grp]
    X=cj15['Q2'][Q2]['x']
    if flav!='du' and flav!='dbub':
      central=cj15['Q2'][Q2]['x'+flav]
      error=cj15['Q2'][Q2]['err-x'+flav]*10
    elif flav=='du':
      central=cj15['Q2'][Q2]['d/u']
      error=cj15['Q2'][Q2]['err-d/u']*1
    elif flav=='dbub':
      central=cj15['Q2'][Q2]['db/ub']
      error=cj15['Q2'][Q2]['err-db/ub']*1

    if grp!='fmKP_AV18':
      ax.plot(X,central,grpmap[grp]['color'],label=tex(grp.replace('_','-')))
    else:
      p1,=ax.plot(X,central,grpmap[grp]['color'])
      p2=fill_between(X,central-error,central+error,ax=ax,
        facecolor='#FFFF00',
        edgecolor='#FFFF00',
        alpha=1.0,hatch=None)

  # make legend
  if flav=='d':
    H_,L_ = ax.get_legend_handles_labels()
    H=[(p2,p1)]
    L=[tex('fmKP-AV18')]
    for h in H_: H.append(h)
    for l in L_: L.append(l)
    ax.legend(H,L,loc=1,frameon=0,fontsize=8)
    ax.text(0.1,0.15,tex('NLO'),
      transform=ax.transAxes,size=15)

  # setup axis 
  #if flav!='s' and flav!='g': ax.set_xticks([])
  ax.set_xlim(0.0,1.0)
  if flav=='dbub': ax.set_xlim(0.0,0.4)
  if flav=='ub' or flav=='db' or flav=='s' or flav=='g':
    ax.semilogx()
    ax.set_xlim(1e-3,1.0)

  ax.set_ylim(ymap[flav]['min'],ymap[flav]['max'])
  ax.set_xlabel('$x$',size=20)
  ax.set_ylabel(labmap[flav],size=20)

  # write info
  if flav=='u': 
    ax.text(0.1,0.15,'$Q^2=$'+tex('~10~GeV^2'),
      transform=ax.transAxes,size=15)

#py.show()
py.savefig('plots/CJ15fits.pdf')
