#!/usr/bin/env python
import sys,os
import numpy as np
import time
import fnmatch
import cPickle 
from operator import mul
import pylab as py
import matplotlib.patches as mpatches
#from line_profiler import LineProfiler

def tex(x):
  return r'$\mathrm{'+x+'}$'

def add_subplot_axes(ax,rect,axisbg='w'):
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

class SPLIT_AX(object):

  def __init__(self,ax):
    self.ax=ax
    self.get_LR(ax)
    self.is_ylim_set=False

  def get_LR(self,ax,xlabel='x-label',ylabel='y-label'):

    ax.axis('off')
    axisbg='w'
  
    fig = py.gcf()
    box = ax.get_position()
    transFigure = fig.transFigure.inverted()
    width = box.width/2
    height = box.height
    
    # create axL
    inax_position  = ax.transAxes.transform([0,0])
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    axL = fig.add_axes([x,y,width,height],axisbg=axisbg)
    axL.spines['right'].set_visible(False)
    axL.get_yaxis().tick_left()
  
    # create axR
    inax_position  = ax.transAxes.transform([0.5,0])
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    axR = fig.add_axes([x,y,width,height],axisbg=axisbg)
    axR.get_yaxis().tick_left()
    axR.spines['left'].set_visible(False)
    axR.axes.yaxis.set_ticklabels([])
    axR.axes.get_yaxis().set_ticks([])
  
    self.axL=axL
    self.axR=axR

  def plot(self,X,Y,*args,**kwargs):

    # break the arrays for L&R
    I=-1
    for i in range(len(X)):
      if X[i]>=0.1: 
        I=i
        break
    XL,YL=X[:I+1],Y[:I+1]
    XR,YR=X[I:],Y[I:]

    # plot arrays
    self.axR.plot(XR,YR,*args,**kwargs)
    self.axL.plot(XL,YL,*args,**kwargs) 

    # set y-limits
    y1=np.amin(Y)
    y2=np.amax(Y)

    if self.is_ylim_set==False:
      self.y1_=y1
      self.y2_=y2
      self.is_ylim_set=True
    else:
      self.y1_=np.amin([y1,self.y1_])
      self.y2_=np.amax([y2,self.y2_])

    self.axL.set_ylim(self.y1_,self.y2_)
    self.axR.set_ylim(self.y1_,self.y2_)

    # set x-limits
    self.axL.set_xlim(XL[0],0.1)
    self.axR.set_xlim(0.1,XR[-1])

    self.axR.set_xticks([0.3,0.5,0.7,0.9])
    self.axL.semilogx()

  def set_ylabel(self,text,displace=-0.15,**kwargs):
    self.axL.set_ylabel(text)
    self.axL.yaxis.set_label_coords(displace,0.5)

  def set_xlabel(self,text,displace=-0.1,**kwargs):
    self.axL.set_xlabel(text)
    self.axL.xaxis.set_label_coords(1.0,displace)

  def tick_params(self,*args,**kwargs):
    self.axL.tick_params(*args,**kwargs)
    self.axR.tick_params(*args,**kwargs)

  def set_title(self,*args,**kwargs):
    self.axL.set_title(*args,**kwargs)

  def legend(self,*args,**kwargs):
    if any([k=='loc' for k in kwargs.keys()]):
      if kwargs['loc']==1 or kwargs['loc']==4: self.axR.legend(**kwargs)
      if kwargs['loc']==2 or kwargs['loc']==3: self.axL.legend(**kwargs)
    else:
      self.axR.legend(**kwargs)

  def set_ylim(self,*args):
    self.axL.set_ylim(*args)
    self.axR.set_ylim(*args)

  def axhline(self,**kwargs):
    self.axL.axhline(**kwargs)
    self.axR.axhline(**kwargs)


if  __name__=="__main__":

  ax=py.subplot(111)
  SA=SPLIT_AX(ax)
  X=10**np.linspace(-5,-1,100)
  X=np.append(X,np.linspace(0.1,1,100))
  Y=X*(1-X)
  SA.plot(X,Y)
  py.savefig('plot.pdf')













