#!/usr/bin/env python
import sys,os
sys.path.insert(1,'../../')
import numpy as np
import pylab as py
# Need to "source setup.sh" before running LHAPDF
import lhapdf
from tools import tex, fill_between
import matplotlib.pyplot as plt
plt.rc('text',usetex=True)


# Initialize PDF instances
#CJ12=lhapdf.mkPDFs('CJ12min')
#CJ15=lhapdf.mkPDFs('CJ15')
HERA15=lhapdf.mkPDFs('HERAPDF15LO_EIG')
MMHT14=lhapdf.mkPDFs('MMHT2014lo68cl')
MSTW08=lhapdf.mkPDFs('MSTW2008lo68cl')
#NNPDF1=lhapdf.mkPDFs('NNPDF30_lo_as_0118')
#NNPDF2=lhapdf.mkPDFs('NNPDF30_lo_as_0130')


# Read CJ15 PDFs from file
# define filepath
fname='CDBONN_KP_1_LO_6.pdf'

# Load file into list L
F=open(fname,'r')
L=F.readlines()
F.close()

# Isolate data at specific Q2 (e.g. 10 GeV^2)
TL=[]
flag=False
for l in L:
  if 'Q2=' in l and '10.00' in l: flag=True
  if 'Q2=' in l and '25.00' in l: flag=False
  if flag==True:
    TL.append(l) 
L=TL

# Remove spaces, newlines etc
L=[l.strip() for l in L]
# Split each line separated by spaces
L=[l.split() for l in L]
# Remove empty lists
L=[l for l in L if l!=[]]

# Get headers
H=L[1]

# Get matrix of data
data=L[2:]
data=[[float(x) for x in l] for l in data ]
data=np.array(data)
data=np.transpose(data)
# data = [[values of 'X'],[ values of 'xu'],...] 

# Construct CJ15 dictionary
CJ15={}
for i in range(len(H)):
  CJ15[H[i]] = data[i]


# Dictionary for PDFs
flavmap = {}
flavmap['u'] = 2
flavmap['d'] = 1
flavmap['ub'] = -2
flavmap['db'] = -1
flavmap['s'] = 3
flavmap['g'] = 21

get_pdfs = lambda flav,x,Q2,iset,grp: grp[iset].xfxQ2(flavmap[flav],x,Q2)


# Dictionary for labels
L = {}
L['u']  = '$xu(x)$'
L['d']  = '$xd(x)$'
L['ub'] = r'$x\bar{u}(x)$'
L['db'] = r'$x\bar{d}(x)$'
L['s']  = '$xs(x)$'
L['g']  = '$xg(x)$'

grpmap = {}
#grpmap['CJ12'] = {'grp':CJ12,'color':'r-'}
grpmap['HERA15'] = {'grp':HERA15,'color':'g-'}
grpmap['MMHT14'] = {'grp':MMHT14,'color':'b-'}
grpmap['MSTW08'] = {'grp':MSTW08,'color':'b:'}

yaxis={}
yaxis['u'] ={'min':0.0,'max':0.8}
yaxis['d'] ={'min':0.0,'max':0.6}
yaxis['ub']={'min':0.0,'max':0.6}
yaxis['db']={'min':0.0,'max':0.6}
yaxis['s'] ={'min':0.0,'max':0.6}
yaxis['g'] ={'min':0.0,'max':20.0}


# Create x-array
X=np.linspace(1e-3,0.9,1000)

Q2=10.0
iset=0


# Make plot
cnt=0
for flav in ['u','d','ub','db','s','g']:
  cnt+=1
  ax=py.subplot(3,2,cnt)    # (rows,columns,index)

  for grp in ['CJ15','HERA15','MMHT14','MSTW08']:

    if grp=='CJ15':
      central=CJ15['x'+flav]
      error=CJ15['Dx'+flav]*10
      fill_between(CJ15['X'],central-error,central+error,ax=ax,\
        facecolor='#FFFF00',edgecolor='#FFFF00',alpha=1.0,hatch=None)
      ax.plot(CJ15['X'],central,'r-',label=tex(grp))

    else:

      grp_=grpmap[grp]['grp']
      col=grpmap[grp]['color']
      ax.plot(X,[get_pdfs(flav,x,Q2,iset,grp_) for x in X],col,label=tex(grp))

    if flav!='u' and flav!='d':
      ax.semilogx()
      ax.set_xlim(1e-3,1.0)

    ax.set_ylim(yaxis[flav]['min'],yaxis[flav]['max'])

    if cnt==1: ax.legend(loc=1,frameon=0,fontsize=20)
    if cnt==1: ax.text(0.1,0.1,tex('LO'),transform=ax.transAxes,size=20)
    if cnt==2: ax.text(0.5,0.7,r'$Q^2=10\;GeV^2$',transform=ax.transAxes,size=20)
  ax.set_xlabel('$x$',size=20)
  ax.set_ylabel(L[flav],size=20)

py.tight_layout()
py.savefig('plots/LOfits.pdf')
#py.show()
