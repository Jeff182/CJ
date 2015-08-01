#!/usr/bin/env python
import sys,os
sys.path.insert(1,'../../')
import numpy as np
import pylab as py
import pandas as pd
from tools import tex
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import lhapdf

class StructFunc(object):

  def __init__(self):
    self.path2exp='../../expdata/'
    self.D={}
    self.get_expdata_fnames()
    self.load_expdata()
    self.split_data()
    #self.get_errors()

  def get_expdata_fnames(self):
    D=self.D

    D['HERA p']={'fname':'HERA.dat', 'label':'HERA'}
    D['BCDMS p']={'fname':'BcdF2pCor', 'label':'BCDMS'}
    #D['BCDMS d']={'fname':'BcdF2dCor'}
    D['NMC p']  ={'fname':'NmcF2pCor', 'label':'NMC'}
    D['SLAC p'] ={'fname':'slac_p_reb', 'label':'SLAC'}
    #D['SLAC d'] ={'fname':'slac_d_reb'}
    D['JLab p'] ={'fname':'jl00106F2p', 'label':'JLab'}
    #D['JLab d'] ={'fname':'jl00106F2d'}
    ##D['BNS d']  ={'fname':'BNS_F2nd'}

    D['HERA p']['color']='g'
    D['BCDMS p']['color']='b'
    #D['BCDMS d']['color']='g'
    D['NMC p']['color']  ='c'
    D['SLAC p']['color'] ='y'
    #D['SLAC d']['color'] ='y'
    D['JLab p']['color'] ='r'
    #D['JLab d']['color'] ='k' 
    ##D['BNS d']['color'] ='r' 

    D['HERA p']['symbol']='*'
    D['BCDMS p']['symbol']='.'
    #D['BCDMS d']['symbol']='s'
    D['NMC p']['symbol']  ='*'
    D['SLAC p']['symbol'] ='>'
    #D['SLAC d']['symbol'] ='<'
    D['JLab p']['symbol'] ='^'
    #D['JLab d']['symbol'] ='v' 
    ##D['BNS d']['symbol']  ='*' 

    # erros
    D['HERA p']['ERR-key']=['STAT+SYST']
    D['BCDMS p']['ERR-key']=['STAERR','SYSERT']
    #D['BCDMS d']['ERR-key']=['STAERR','SYSERT']
    D['NMC p']['ERR-key']=['STAERR','SYSERT']
    D['SLAC p']['ERR-key']=['STAT.','SYS.']
    #D['SLAC d']['ERR-key']=['STAT','SYS']
    D['JLab p']['ERR-key']=['STAT','SYST']
    #D['JLab d']['ERR-key']=['STAT','SYST']

    keys=[]
    keys.append('HERA p')
    keys.append('BCDMS p')
    #keys.append('BCDMS d')
    keys.append('NMC p')
    keys.append('SLAC p')
    #keys.append('SLAC d')
    keys.append('JLab p')
    #keys.append('JLab d')
    self.ordered_keys=keys

  def load_expdata(self):  
    D=self.D

    for k in D.keys():
      print 'loading ',k

      # open file
      F=open(self.path2exp+D[k]['fname'])
      L=F.readlines()
      F.close()
      L=[l.strip() for l in L]
      L=[l.split() for l in L if l!='']

      # construct table
      table=[]
      flag=False
      for i in range(len(L)): 
        try:
          l=[float(x) for x in L[i]]
          if flag==False:
            ih=i-1
            flag=True
        except:
          continue
        table.append(l)
      table=np.transpose(table)

      # construct headers
      H=[x.upper() for x in L[ih]]
      for i in range(len(H)): H[i]=H[i].replace('Q**2','Q2')
      for i in range(len(H)): H[i]=H[i].replace('Q^2','Q2')
      for i in range(len(H)): H[i]=H[i].replace('F2P','F2')
      for i in range(len(H)): H[i]=H[i].replace('F2D','F2')
      for i in range(len(H)): H[i]=H[i].replace('F2D','F2')

      # construct pandas data frame
      d={}
      for i in range(len(H)):
        d[H[i]]=table[i]
      d['W2']=0.9389185**2 + d['Q2']/d['X'] - d['Q2']

      ERR=np.zeros(d['W2'].size)
      for kk in D[k]['ERR-key']:
        ERR+=d[kk]**2
      ERR=ERR**0.5
      d['ERR']=ERR

      DF=pd.DataFrame(d)
      #DF=DF[DF.W2>4.0]
      DF=DF[DF.Q2>1.0]

      # store DF in global dic
      D[k]['DF']=DF

  def get_xbins(self):
    xbins=[]


    xbins.append([3.8e-6,4.5e-6])
    xbins.append([4.9e-6,5.7e-6])
    xbins.append([6.3e-6,7.2e-6])
    xbins.append([7.9e-6,9.0e-6])
    xbins.append([9.5e-6,11.1e-6])
    xbins.append([12e-6,13.9e-6])
    xbins.append([15.2e-6,17.0e-6])
    xbins.append([19e-6,21e-6])
    xbins.append([31e-6,33e-6])
    xbins.append([38e-6,42e-6])
    xbins.append([48e-6,52e-6])
    xbins.append([58e-6,67e-6])
    xbins.append([77e-6,82e-6])
    xbins.append([96e-6,11e-5])
    xbins.append([12.5e-5,13.5e-5])
    xbins.append([19.5e-5,20.5e-5])
    xbins.append([24.5e-5,25.5e-5])
    xbins.append([30.5e-5,32.5e-5])
    xbins.append([47.5e-5,50.5e-5])
    xbins.append([78.5e-5,82.5e-5])
    xbins.append([12.5e-4,13.5e-4])
    xbins.append([19.5e-4,20.5e-4])

    xbins.append([3.1e-3,3.8e-3])
    xbins.append([4.8e-3,5.7e-3])
    xbins.append([7.2e-3,9.3e-3])
    xbins.append([1.15e-2,1.4e-2])
    xbins.append([1.65e-2,1.9e-2])
    xbins.append([1.95e-2,2.05e-2])
    xbins.append([2.3e-2,2.9e-2])
    xbins.append([3.15e-2,3.25e-2])
    xbins.append([3.4e-2,3.8e-2])
    xbins.append([4.65e-2,5.4e-2])
    xbins.append([6.5e-2,7.3e-2])
    xbins.append([7.8e-2,8.2e-2])
    xbins.append([8.5e-2,9.2e-2])
    xbins.append([9.8e-2,10.3e-2])
    xbins.append([10.8e-2,11.3e-2])
    xbins.append([12.8e-2,13.2e-2])
    xbins.append([13.6e-2,14.6e-2])
    xbins.append([17.1e-2,18.7e-2])
    xbins.append([19.7e-2,20.7e-2])
    xbins.append([21.7e-2,23.7e-2])
    xbins.append([24.8e-2,25.2e-2])
    xbins.append([26.0e-2,29.0e-2])
    xbins.append([33.0e-2,36.0e-2])
    xbins.append([39.5e-2,40.5e-2])
    xbins.append([42.0e-2,48.0e-2])

    #xbins.append([48.1e-2,49.4e-2])
    #xbins.append([49.4e-2,50.8e-2])
    #xbins.append([50.8e-2,51.8e-2])
    #xbins.append([51.8e-2,52.6e-2])
    #xbins.append([52.6e-2,53.4e-2])
    #xbins.append([53.4e-2,53.9e-2])
    #xbins.append([53.9e-2,54.5e-2])
    #xbins.append([54.5e-2,55.6e-2])
    #xbins.append([55.6e-2,56.6e-2])
    #xbins.append([56.6e-2,57.3e-2])
    #xbins.append([53.0e-2,56.0e-2])

    xbins.append([48.0e-2,63.0e-2])

    xbins.append([63.0e-2,66.0e-2])
    xbins.append([73.0e-2,76.0e-2])
    xbins.append([84.0e-2,86.0e-2])
    return xbins[::-1]

  def split_data(self):
    xbins=self.get_xbins()
    D=self.D
    for k in D.keys():
      d=D[k]['DF']
      dbinned={}
      for i in range(len(xbins)):
        xmin,xmax=xbins[i]
        dbinned[i]=d[d.X>xmin]
        dbinned[i]=dbinned[i][dbinned[i].X<xmax]
      D[k]['dbinned']=dbinned

  def make_XQ2_plot(self):
    D=self.D
    ax=py.subplot(111)
    for k in self.D.keys():
      d=D[k]['DF']
      if 'HERA' in k: color='k'
      else: color='r'
      ax.plot(d['X'],d['Q2'],color+'.',markersize=3)
    ax.semilogx()
    ax.semilogy()
    xbins=np.array(self.get_xbins()).flatten()
    ax.set_xticks(xbins)
    ax.set_xlim(1e-1,1e-0)
    ax.grid()
    py.savefig('XQ2.pdf')

  def add_subplot_axes(self,ax,rect,axisbg='w'):
    fig = py.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax
 
  def make_F2_plot(self):
    D=self.D
    xbins=self.get_xbins()
    py.figure(figsize=(15,15))

    gs = gridspec.GridSpec(1,1)
    gs.update(left=0.1,right=0.98,top=0.98,bottom=0.1)
    ax=py.subplot(gs[0,0])

    for k in self.ordered_keys:
      dbinned=D[k]['dbinned']
      color=D[k]['color']
      sym=D[k]['symbol']
      flag=False
      for i in range(len(xbins)):
        dbin=dbinned[i]
        if dbin['X'].size!=0:
          if flag==False:
            ax.errorbar(dbin['Q2'],dbin['F2']*2**(i+1)\
              ,yerr=dbin['ERR'].values\
              ,fmt=color+sym\
              ,mec=color
              ,label=tex(D[k]['label']))
            flag=True
          else:
            ax.errorbar(dbin['Q2'],dbin['F2']*2**(i+1)\
              ,yerr=dbin['ERR'].values\
              ,fmt=color+sym\
              ,mec=color)

    for i in range(len(xbins)):
      Q2=1
      for k in D.keys():
        dbin=D[k]['dbinned'][i]
        if dbin['X'].size==0:continue
        imax=np.argmax(dbin['Q2'].values)
        if dbin['Q2'].values[imax]>Q2: 
          Q2=dbin['Q2'].values[imax]
          F2=dbin['F2'].values[imax]

      if any([i==k for k in [10,34,37,39,41,44,49,47,43,35,32,30,28,26,45,46]])!=True:
	if i>30:
          text='$x=%0.1e$'%(np.mean([xbins[i][0],xbins[i][1]]))
	  exp=int(text.split('-')[1].replace('$',''))
	  text=text.split('e')[0]+'e'
          text=text.replace('e',r'\times 10^{-%d}'%exp)+'\ (i=%d)$'%i
	elif i<16:
	  text='$x=%0.2f'%(np.mean([xbins[i][0],xbins[i][1]]))+'\ (i=%d)$'%i
	else:
	  text='$x=%0.3f'%(np.mean([xbins[i][0],xbins[i][1]]))+'\ (i=%d)$'%i

        ax.text(Q2*1.2,F2*2**(i+1),text)
        #ax.text(Q2*1.2,F2*2**(i+1),text+str(i))

    ax.legend(frameon=0,fontsize=20,numpoints=1)

    ax.set_xlim(8e-5,3e5)
    ax.set_ylim(1e-3,2e13)
    ax.semilogy()
    ax.semilogx()
    ax.set_xticks([1e-1,1,1e1,1e2,1e3,1e4,1e5])
    #ax.set_ylabel(tex('F_2')+'$(x,Q^2)$',size=30)
    ax.set_ylabel('$F_2^p(x,Q^2)\ * 2^{\ i}$',size=30)
    ax.set_xlabel('$Q^2$'+tex('\ (GeV^2)'),size=30)
    py.tick_params(axis='both',labelsize=20)

    xsq=0.45
    ysq=0.07
    ax.add_patch(patches.Rectangle((xsq,ysq),0.1,0.15\
      ,fill=False,transform=ax.transAxes))
    ax.plot([0.36,xsq],[0.07,ysq],'k:'
      ,transform=ax.transAxes)
    ax.plot([0.36,xsq],[0.47,ysq+0.15],'k:'
      ,transform=ax.transAxes)


    rect1 = [0.06,0.07,0.3,0.3]
    ax1 = self.add_subplot_axes(ax,rect1)
    rect2 = [0.06,0.47,0.3,0.3]
    ax2 = self.add_subplot_axes(ax,rect2)

    for k in D.keys():
      if 'JLab' not in k: continue
      d=D[k]['DF']
      color=D[k]['color']
      sym=D[k]['symbol']

      ax1.errorbar(d['Q2'],d['F2']\
        ,yerr=d['ERR'].values\
        ,fmt=color+sym\
        ,mec=color)

      ax2.errorbar(d['X'],d['F2']\
        ,yerr=d['ERR'].values\
        ,fmt=color+sym\
        ,mec=color)

    ax1.locator_params(nbins=5) 
    ax2.locator_params(nbins=5) 

    ax1.tick_params(axis='both',labelsize=12)
    ax2.tick_params(axis='both',labelsize=12)

    ax1.set_xlabel(r'$Q^2$'+tex('\ (GeV^2)'),size=20)
    ax2.set_xlabel(r'$x$',size=20)
    ax1.set_ylim(1e-3,0.25)
    ax2.set_ylim(1e-3,0.25)
    ax1.semilogy()
    ax2.semilogy()
    #py.tight_layout()
    py.savefig('F2p.pdf')

if __name__=='__main__':

  SF=StructFunc()
  #SF.make_XQ2_plot()
  SF.make_F2_plot()

  #CJ150=lhapdf.mkPDF('CJ15_NLO',0)
  #print CJ150.xfxQ2(2,0.5,100)

