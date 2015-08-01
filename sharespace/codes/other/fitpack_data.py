#!/usr/bin/env python
import sys,os
sys.path.insert(1,'../../')
import numpy as np
from master import FITPACK


fitpack=FITPACK()

data=fitpack.get_PDFs('CDBONN_KP_1_LO_6.pdf')



#import numpy as np
#import pylab as py
#
## define filepath
#path='/u/group/cteqX/wmelnitc/fits/CJ15'
#fname='CDBONN_KP_1_LO_6.pdf'
#
## load file into list L
#F=open(fname,'r')
#L=F.readlines()
#F.close()
#
## isolate Q2=10.0 data
#TL=[]
#flag=False
#for l in L:
#  if 'Q2=' in l and '10.00' in l: flag=True
#  if 'Q2=' in l and '25.00' in l: flag=False
#  if flag==True:
#    TL.append(l) 
#L=TL
#
## remove spaces, newlines etc
#L=[l.strip() for l in L]
## split each line separated by spaces
#L=[l.split() for l in L]
## rm empty lists
#L=[l for l in L if l!=[]]
#
## get headers
#H=L[1]
#
## get matrix of data
#data=L[2:]
#data=[[float(x) for x in l] for l in data ]
#data=np.array(data)
#data=np.transpose(data)
## data = [[values of 'X'],[ values of 'xu'],...] 
#
## construct dictionary
#D={}
#for i in range(len(H)):
#  D[H[i]] = data[i]
#
##print D.keys()
##print D['xg']
#
#
#ax=py.subplot(111)
#ax.errorbar(D['X'],D['xu'],yerr=10*D['Dxu'],fmt='k.')
#py.show()

















