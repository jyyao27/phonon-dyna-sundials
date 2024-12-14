#!/usr/bin/env python3

"""
Set of utils for the PERTURBO postprocessing
"""

import os, re, sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from collections    import defaultdict
from scipy import stats
import scipy.integrate as integrate
from scipy.signal import argrelextrema
from yaml        import load,dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper



####################################################################################################
def kpts_fold(kpts,bg):

   print('Folding k points into the first BZ...')

   b1 = bg[0,:]
   b2 = bg[1,:]
   b3 = bg[2,:]

   #b1 = np.array([-1.0,-1.0,1.0])
   #b2 = np.array([1.0,1.0,1.0])
   #b3 = np.array([-1.0,1.0,-1.0])
   
   loop_region = np.array([-1,0,1])
   
   kpts_out = kpts.copy()
   
   for i in loop_region:
      for j in loop_region:
         for k in loop_region:
            
            # linear combination of vectors
            tmp_array = kpts - i*b1 - j*b2 - k*b3

            # create a logical array, where the norm is smaller
            condition_array = np.linalg.norm(kpts_out,axis=1) < np.linalg.norm(tmp_array,axis=1)
            
            # duplicate the columns to match the kpts array shape
            condition_array = np.array([condition_array]*3).T

            kpts_out = np.where( condition_array, kpts_out, tmp_array )


   return kpts_out

   for ik in range(len(kpts)):
       kpt = kpts[ik]
       dismin = 1e6
       for i in loop_region:
           for j in loop_region:
               for k in loop_region:
                   tmp = kpt - i*b1 - j*b2 - k*b3
                   dis = np.linalg.norm(tmp)
                   if (dis < dismin):
                      dismin = dis
                      i_sel = i
                      j_sel = j
                      k_sel = k
       kpts_out[ik] = kpts_out[ik] - i_sel*b1 - j_sel*b2 - k_sel*b3

