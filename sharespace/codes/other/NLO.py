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
HERA15 =lhapdf.mkPDFs('HERAPDF15NLO_EIG')
MMHT14 =lhapdf.mkPDFs('MMHT2014nlo68cl')
#MSTW08 =lhapdf.mkPDFs('MSTW2008nlo68cl')
JR14   =lhapdf.mkPDFs('JR14NLO08VF')
CJ12min=lhapdf.mkPDFs('CJ12min')
CT10   =lhapdf.mkPDFs('CT10nlo')

# dictionary for PDG flav index
pdgmap= {'u':2,'d':1,'ub':-2,'db': -1,'s':3,'g':21}

# function to access LHAPDF pdfs
def get_pdfs(flav,x,Q2,iset,grp):

  if flav!='du' and flav!='dbub': 
    return grp[iset].xfxQ2(pdgmap[flav],x,Q2)

  elif flav=='du':

    if grp[iset].xfxQ2(pdgmap['u'],x,Q2)!=0.0:
      return grp[iset].xfxQ2(pdgmap['d'],x,Q2)/grp[iset].xfxQ2(pdgmap['u'],x,Q2)
    else:
      return 0.0

  elif flav=='dbub':

    if grp[iset].xfxQ2(pdgmap['ub'],x,Q2)!=0.0:
      return grp[iset].xfxQ2(pdgmap['db'],x,Q2)/grp[iset].xfxQ2(pdgmap['ub'],x,Q2)
    else:
      return 0.0

  else:
    print '*** ERROR: flav not defined ***'
    sys.exit()


# load CJ15 from fitpack data 
CJ15 = {}
CJ15['CJ15_AV18']   = FITPACK().get_PDFs('data/CJ15_NLO_KP_AV18.pdf')

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
grpmap['HERA15'] = {'grp':HERA15,'color':'g--'}
grpmap['JR14']   = {'grp':JR14,'color':'g-'}
grpmap['MMHT14'] = {'grp':MMHT14,'color':'b-'}
#grpmap['MSTW08'] = {'grp':MSTW08,'color':'b:'}
grpmap['CJ12min']= {'grp':CJ12min,'color':'k-.'}
grpmap['CT10']   = {'grp':CT10,'color':'r-.'}

# dictionary for ylims 
ymap = {}
ymap['u'] ={'min':0.0,'max':0.8}
ymap['d'] ={'min':0.0,'max':0.8}
ymap['du']  ={'min':0.0,'max':1.0}
ymap['dbub']={'min':0.6,'max':2.0}
ymap['ub']={'min':0.0,'max':0.6}
ymap['db']={'min':0.0,'max':0.6}
ymap['s'] ={'min':0.0,'max':0.6}
ymap['g'] ={'min':0.0,'max':10.0}

# dictionary for plot location
ncols=2
nrows=4

py.figure(figsize=(ncols*4,nrows*2)) 
gs = gridspec.GridSpec(4,2) # specify plotting grid geometry
gs.update(left=0.1,right=0.98,wspace=0.25,hspace=0.45,top=0.98,bottom=0.1)
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

  for grp in ['CJ15','HERA15','MMHT14','JR14','CT10','CJ12min']:

    if grp=='CJ15':
      cj15 = CJ15['CJ15_AV18']
      X=cj15['Q2'][Q2]['x']
      if flav!='du' and flav!='dbub':
        central=cj15['Q2'][Q2]['x'+flav]
        error=cj15['Q2'][Q2]['err-x'+flav]*10
      elif flav=='du':
        central=cj15['Q2'][Q2]['d/u']
        error=cj15['Q2'][Q2]['err-d/u']*10
      elif flav=='dbub':
        central=cj15['Q2'][Q2]['db/ub']
        error=cj15['Q2'][Q2]['err-db/ub']*10

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
      X=np.linspace(1e-3,0.95,1000)
      grp_=grpmap[grp]['grp']
      col=grpmap[grp]['color']
      ax.plot(X,[get_pdfs(flav,x,Q2,iset,grp_) for x in X],col,label=tex(grp))

  # make legend
  if flav=='d':
    H_,L_ = ax.get_legend_handles_labels()
    H=[(p2,p1)]
    L=[tex('CJ15')]
    for h in H_: H.append(h)
    for l in L_: L.append(l)
    ax.legend(H,L,loc=1,frameon=0,fontsize=9)
    ax.text(0.07,0.15,tex('NLO'),transform=ax.transAxes,size=15)


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
    ax.text(0.55,0.75,'$Q^2=$'+tex('~10~GeV^2'),
      transform=ax.transAxes,size=12)

#py.show()
py.savefig('plots/NLOfits.pdf')