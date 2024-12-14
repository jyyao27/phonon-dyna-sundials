#!/usr/bin/env python3

import os, sys, datetime
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter

from mpl_toolkits.axes_grid1 import make_axes_locatable

import scipy.interpolate as interpolate
import scipy.optimize    as optimize
from scipy.interpolate import griddata

import perturbo_postproc_utils as ppu

"""
Associate the k-points to a given valley using cc3d

cc3d needs the 3D array with some values 1 and other 0. It will connect the 1s. 
The scipt converts the array of k-point vectors to this 3D array and then uses cc3d

Run in:
/Users/imaliyov/run_external_field/Si/vel_field/cdyna-elec/NERSC/EPWAN_TUTORIAL/k_160_q_160_dt_10.0_dE_200_dir_100_c_1.0e18_emax_7.1
"""

#####################################################################################
def gaussian(x,mu,sig): # Gaussian, mu = energy center
    
    return np.exp(-0.5*((x - mu)/sig)**2 ) / ( sig * np.sqrt(2*np.pi) )

#####################################################################################
def get_popu(dist,energies,e_scale,eta=5e-2): 
   # eta is smearing in eV (typical: 5e.-3)
   
   popu = np.zeros_like(e_scale) 
   
   for iband in range(energies.shape[1]):
      for ien,en in enumerate(e_scale):
          popu[ien] = 1/N * np.sum(dist[:] * gaussian( (energies[:,iband] - en ), 0.0, eta ))
            
   return popu

#####################################################################################################
# This is called everytime you release the mouse button
def on_click(event):
    azim, elev = ax.azim, ax.elev
    print(azim, elev)

#####################################################################################################
def border_minus_one(k,n):
   k_out = k
   for i,kk in enumerate(k):
      if kk == n:
         k_out[i] = kk-1
   return k_out

#####################################################################################################
def draw_axes(ax,**kwargs):

   x1 = -1.5
   x2 = 1.5

   ax.plot([x1,x2],[0,0],[0,0],**kwargs)
   ax.plot([0,0],[x1,x2],[0,0],**kwargs)
   ax.plot([0,0],[0,0],[x1,x2],**kwargs)

#####################################################################################################
def draw_BZ_border_2D(ax,BZ_array, center, omit_dir, order_dir = None, **kwargs):
   BZ_array = BZ_array + center
   #BZ_array = BZ_array[np.where(BZ_array[:,omit_dir] == 0)[0],:]
   #BZ_array  = np.array([
   #   [0.5, 1.0, 0.0],
   #   [1.0, 0.5, 0.0],
   #   [1.0, -0.5, 0.0],
   #   [0.5, -1.0, 0.0],
   #   [-0.5, -1.0, 0.0],
   #   [-1.0, -0.5, 0.0],
   #   [-1.0, 0.5, 0.0],
   #   [-0.5, 1.0, 0.0],
   #   [0.5, 1.0, 0.0]]) 
   x1 = np.mod(omit_dir + 1, 3)
   x2 = np.mod(omit_dir + 2, 3)
   #x1 = 0
   #x2 = 1
   if order_dir is not None:
      
      for i in range(BZ_array.shape[0]-1):
         
         # foreground
         zf = 10
         lsf = '-'

         # background
         zb = 0
         lsb = '--'
         
         if order_dir.shape[0] == 3:
            dirvec = order_dir/np.linalg.norm(order_dir)
         elif order_dir.shape[0] == 2:
            azim = order_dir[0]/180.0*np.pi
            elev = order_dir[1]/180.0*np.pi
            r = 1.0

            x = r * np.cos(elev) * np.cos(azim)
            y = r * np.cos(elev) * np.sin(azim)
            z = r * np.sin(elev)

            dirvec = np.array([x,y,z])

            #ax.plot([0,dirvec[0]], [0,dirvec[1]], [0,dirvec[2]], c='red')

         else:
            sys.exit('draw_BZ_border: Wrong order_dir')

         lim = -0.4

         r1 = BZ_array[i,:]
         r2 = BZ_array[i+1,:]

         p1 = np.dot(dirvec,r1)
         p2 = np.dot(dirvec,r2)

         if p1 >= lim and p2 >= lim:
            zorder = zf
            ls = lsf
            c = 'red'
         else:
            zorder = zb
            ls = lsb
            c = 'blue'
         ax.plot( [r1[x1],r2[x1]], [r1[x2],r2[x2]], zorder=zorder, ls = ls,**kwargs )
         #ax.plot( [r1[0],r2[0]], [r1[1],r2[1]], [r1[2],r2[2]], zorder=zorder, ls = ls,**kwargs )
         
         # debug
         #ax.plot( [r1[0],r2[0]], [r1[1],r2[1]], [r1[2],r2[2]], zorder=zorder, c=c )

   else:
      
      ax.plot(BZ_array[:,x1],BZ_array[:,x2],**kwargs)


#####################################################################################################

def plot_border_2D(fig, ax, folder, omit_dir, centers):
   
   highsym_dict_Si = {
      'Xc1': {
         'color' : 'blue', 'ls':'-',
         #'color' : 'black',
         'radius': 0.24,
         'coord' : np.array([
                  [-0.5, 0.0,-0.5 ],
                  ])
      },
      
      'Xc2': {
         #'color' : 'black',
         'color' : 'blue', 'ls':'--',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.5, 0.0, 0.5 ],
                  ])
      },
      
      'Xh1': {
         #'color' : 'black', 
         'ls':'-',
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.5, 0.5, 0.0 ],
                  ])
      },
      
      'Xh2': {
         #'color' : 'black',
         'ls':'--',
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [-0.5,-0.5, 0.0 ],
                  ])
      },
      
      'Xh3': {
         #'color' : 'black',
         'ls':'-.',
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.0,-0.5,-0.5 ],
                  ])
      },
      
      'Xh4': {
         #'color' : 'black',
         'ls':'dotted',
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.0, 0.5, 0.5 ],
                  ])
      },
   }
   
   BZ_array = np.array([[ 0.25, -0.25, -0.5 ],
          [-0.25, -0.5 , -0.75],
          [-0.5 , -0.25, -0.75],
          [-0.25,  0.25, -0.5 ],
          [ 0.25,  0.5 , -0.25],
          [ 0.5 ,  0.25, -0.25],
          [ 0.25, -0.25, -0.5 ],
          [ 0.25, -0.5 , -0.25],
          [-0.25, -0.75, -0.5 ],
          [-0.25, -0.5 , -0.75],
          [-0.25, -0.75, -0.5 ],
          [-0.5 , -0.75, -0.25],
          [-0.25, -0.5 ,  0.25],
          [ 0.25, -0.25,  0.5 ],
          [ 0.5 , -0.25,  0.25],
          [ 0.25, -0.5 , -0.25],
          [ 0.5 , -0.25,  0.25],
          [ 0.75,  0.25,  0.5 ],
          [ 0.75,  0.5 ,  0.25],
          [ 0.5 ,  0.25, -0.25],
          [ 0.75,  0.5 ,  0.25],
          [ 0.5 ,  0.75,  0.25],
          [ 0.25,  0.5 , -0.25],
          [ 0.5 ,  0.75,  0.25],
          [ 0.25,  0.75,  0.5 ],
          [-0.25,  0.5 ,  0.25],
          [-0.5 ,  0.25, -0.25],
          [-0.25,  0.25, -0.5 ],
          [-0.5 ,  0.25, -0.25],
          [-0.75, -0.25, -0.5 ],
          [-0.5 , -0.25, -0.75],
          [-0.75, -0.25, -0.5 ],
          [-0.75, -0.5 , -0.25],
          [-0.5 , -0.75, -0.25],
          [-0.25, -0.5 ,  0.25],
          [-0.5 , -0.25,  0.25],
          [-0.75, -0.5 , -0.25],
          [-0.5 , -0.25,  0.25],
          [-0.25,  0.25,  0.5 ],
          [ 0.25,  0.5 ,  0.75],
          [ 0.5 ,  0.25,  0.75],
          [ 0.25, -0.25,  0.5 ],
          [ 0.5 , -0.25,  0.25],
          [ 0.75,  0.25,  0.5 ],
          [ 0.5 ,  0.25,  0.75],
          [ 0.25,  0.5 ,  0.75],
          [ 0.25,  0.75,  0.5 ],
          [-0.25,  0.5 ,  0.25],
          [-0.25,  0.25,  0.5 ]])
   
   
   #
   # Read k-points of the full BZ with the energy window
   #
   #prefix='gaas'
   prefix='si'
   
   highsym_dict = highsym_dict_Si
   epwan_filename = 'si_epr.h5' 
   epwan_file = h5py.File(epwan_filename, 'r')
   # bg - reciprocal lattice vectors in unit of $2\pi/a$ where $a$ is the lattice constant
   bg = epwan_file['basic_data']['bg'][()]
   alat = epwan_file['basic_data']['alat'][()]
   #bg = bg * 2 * np.pi / alat
   epwan_file.close()
   
   #
   # !!! transform kpts and BZ to cartesian coordinates
   #
   
   #tet_filename = f'{folder}/si_tet.h5'
   tet_filename = f'{folder}/si_qg_tet.h5'

   #tet_filename = 'r40-mri-50-10ps/si_qg_tet.h5'
   #tet_filename = 'r80-mri-50-30ps/si_qg_tet.h5'
   tet_file = h5py.File(tet_filename, 'r')
   kpts = tet_file['kpts_all_crys_coord'][()]
   tet_file.close()

   kpts = np.matmul(kpts, bg)
   BZ_array = np.matmul(BZ_array, bg)

   # Fold kpts into the first BZ
   kpts = ppu.kpts_fold(kpts,bg)
   #
   # 3D array of zeros
   #
   ngrid = 160
   dim=np.zeros(3,dtype=np.int32)
   #labels_in = np.zeros((ngrid,ngrid,ngrid), dtype=np.int32)
   #
   # !!! Because of the folding, the indicies can be < 0 or > ngrid !!!
   #
   labels_in = np.zeros(tuple(dim), dtype=np.int32)
   
   #
   # Convert to integers with the origin in 0,0,0
   #
   shift_int_array = np.array([ngrid/2,ngrid/2,ngrid/2],dtype=np.int32)
   
   connectivity = 6 # only 4,8 (2D) and 26, 18, and 6 (3D) are allowed
   
   
   plot_type = 'border'
   #plot_type = 'threshold'
   #plot_type = 'double_threshold'
   #ax.xaxis.pane.fill = False
   #azim, elev = -82.63, 18.83
   azim, elev = 45, 19
   #azim, elev = -82.78104903555192, 8.105518475814904
    
   for i_center in range(len(centers[:,0])):
      draw_BZ_border_2D(ax,BZ_array, centers[i_center,:], omit_dir, order_dir=None, c='grey', alpha=0.8, lw=1)

   
   return fig, ax#, kpts

#####################################################################################################
def draw_BZ_border(ax,BZ_array,order_dir = None, **kwargs):

   if order_dir is not None:
      
      for i in range(BZ_array.shape[0]-1):
         
         # foreground
         zf = 10
         lsf = '-'

         # background
         zb = 0
         lsb = '--'
         
         if order_dir.shape[0] == 3:
            dirvec = order_dir/np.linalg.norm(order_dir)
         elif order_dir.shape[0] == 2:
            azim = order_dir[0]/180.0*np.pi
            elev = order_dir[1]/180.0*np.pi
            r = 1.0

            x = r * np.cos(elev) * np.cos(azim)
            y = r * np.cos(elev) * np.sin(azim)
            z = r * np.sin(elev)

            dirvec = np.array([x,y,z])

            #ax.plot([0,dirvec[0]], [0,dirvec[1]], [0,dirvec[2]], c='red')

         else:
            sys.exit('draw_BZ_border: Wrong order_dir')

         lim = -0.4

         r1 = BZ_array[i,:]
         r2 = BZ_array[i+1,:]

         p1 = np.dot(dirvec,r1)
         p2 = np.dot(dirvec,r2)

         if p1 >= lim and p2 >= lim:
            zorder = zf
            ls = lsf
            c = 'red'
         else:
            zorder = zb
            ls = lsb
            c = 'blue'

         ax.plot( [r1[0],r2[0]], [r1[1],r2[1]], [r1[2],r2[2]], zorder=zorder, ls = ls,**kwargs )
         
         # debug
         #ax.plot( [r1[0],r2[0]], [r1[1],r2[1]], [r1[2],r2[2]], zorder=zorder, c=c )

   else:
      ax.plot(BZ_array[:,0],BZ_array[:,1],BZ_array[:,2],**kwargs)


#####################################################################################################

def normk_num(r_array, y_array):
   rmax = np.max(r_array)

   ngrid = 100
   rgrid = np.linspace(0,rmax,ngrid)
   sum_grid = np.zeros(ngrid)

   for i,r in enumerate(rgrid):
      idx_tmp = np.where(r_array < r)

      sum_tmp = np.sum(y_array[idx_tmp])
      sum_grid[i] = sum_tmp

   return rgrid, sum_grid

#####################################################################################################

def get_unique_list(array,atol=1e-6):
   array_sort = np.sort(array)
   unique_list = [array_sort[0]]
   for x in array_sort:
      x_last = unique_list[-1]
      if abs(x_last - x) < atol:
         continue
      unique_list.append(x)

   return unique_list
#####################################################################################################

def plot_border(folder):
   highsym_dict_GaAs = {
   
      'Gamma': {
         'color' : 'black',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.0, 0.0, 0.0]
                  ])
      },
   
      'L': {
         'color' : 'blue',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.0, 0.0,-0.5 ],
                  [ 0.5, 0.0, 0.0 ],
                  [ 0.0, 0.0, 0.5 ],
                  [-0.5, 0.0, 0.0 ],
                  [ 0.0, 0.5, 0.0 ],
                  [ 0.5, 0.5, 0.5 ],
                  [ 0.0,-0.5, 0.0 ],
                  [-0.5,-0.5,-0.5 ],
                  ])
      },
      
      'X': {
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [-0.5, 0.0,-0.5 ],
                  [ 0.5, 0.5, 0.0 ],
                  [ 0.5, 0.0, 0.5 ],
                  [-0.5,-0.5, 0.0 ],
                  [ 0.0,-0.5,-0.5 ],
                  [ 0.0, 0.5, 0.5 ],
                  ])
      }
   }
   
   highsym_dict_Si = {
      'Xc1': {
         'color' : 'blue', 'ls':'-',
         #'color' : 'black',
         'radius': 0.24,
         'coord' : np.array([
                  [-0.5, 0.0,-0.5 ],
                  ])
      },
      
      'Xc2': {
         #'color' : 'black',
         'color' : 'blue', 'ls':'--',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.5, 0.0, 0.5 ],
                  ])
      },
      
      'Xh1': {
         #'color' : 'black', 
         'ls':'-',
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.5, 0.5, 0.0 ],
                  ])
      },
      
      'Xh2': {
         #'color' : 'black',
         'ls':'--',
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [-0.5,-0.5, 0.0 ],
                  ])
      },
      
      'Xh3': {
         #'color' : 'black',
         'ls':'-.',
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.0,-0.5,-0.5 ],
                  ])
      },
      
      'Xh4': {
         #'color' : 'black',
         'ls':'dotted',
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.0, 0.5, 0.5 ],
                  ])
      },
   }
   
   BZ_array = np.array([[ 0.25, -0.25, -0.5 ],
          [-0.25, -0.5 , -0.75],
          [-0.5 , -0.25, -0.75],
          [-0.25,  0.25, -0.5 ],
          [ 0.25,  0.5 , -0.25],
          [ 0.5 ,  0.25, -0.25],
          [ 0.25, -0.25, -0.5 ],
          [ 0.25, -0.5 , -0.25],
          [-0.25, -0.75, -0.5 ],
          [-0.25, -0.5 , -0.75],
          [-0.25, -0.75, -0.5 ],
          [-0.5 , -0.75, -0.25],
          [-0.25, -0.5 ,  0.25],
          [ 0.25, -0.25,  0.5 ],
          [ 0.5 , -0.25,  0.25],
          [ 0.25, -0.5 , -0.25],
          [ 0.5 , -0.25,  0.25],
          [ 0.75,  0.25,  0.5 ],
          [ 0.75,  0.5 ,  0.25],
          [ 0.5 ,  0.25, -0.25],
          [ 0.75,  0.5 ,  0.25],
          [ 0.5 ,  0.75,  0.25],
          [ 0.25,  0.5 , -0.25],
          [ 0.5 ,  0.75,  0.25],
          [ 0.25,  0.75,  0.5 ],
          [-0.25,  0.5 ,  0.25],
          [-0.5 ,  0.25, -0.25],
          [-0.25,  0.25, -0.5 ],
          [-0.5 ,  0.25, -0.25],
          [-0.75, -0.25, -0.5 ],
          [-0.5 , -0.25, -0.75],
          [-0.75, -0.25, -0.5 ],
          [-0.75, -0.5 , -0.25],
          [-0.5 , -0.75, -0.25],
          [-0.25, -0.5 ,  0.25],
          [-0.5 , -0.25,  0.25],
          [-0.75, -0.5 , -0.25],
          [-0.5 , -0.25,  0.25],
          [-0.25,  0.25,  0.5 ],
          [ 0.25,  0.5 ,  0.75],
          [ 0.5 ,  0.25,  0.75],
          [ 0.25, -0.25,  0.5 ],
          [ 0.5 , -0.25,  0.25],
          [ 0.75,  0.25,  0.5 ],
          [ 0.5 ,  0.25,  0.75],
          [ 0.25,  0.5 ,  0.75],
          [ 0.25,  0.75,  0.5 ],
          [-0.25,  0.5 ,  0.25],
          [-0.25,  0.25,  0.5 ]])
   
   
   #
   # Read k-points of the full BZ with the energy window
   #
   #prefix='gaas'
   prefix='si'
   
   highsym_dict = highsym_dict_Si
   epwan_filename = 'si_epr.h5' 
   epwan_file = h5py.File(epwan_filename, 'r')
   # bg - reciprocal lattice vectors in unit of $2\pi/a$ where $a$ is the lattice constant
   bg = epwan_file['basic_data']['bg'][()]
   print(bg)
   alat = epwan_file['basic_data']['alat'][()]
   #bg = bg * 2 * np.pi / alat
   epwan_file.close()
   
   #
   # !!! transform kpts and BZ to cartesian coordinates
   #
   
   
   tet_filename = f'{folder}/si_tet.h5'
   #tet_filename = f'{folder}/si_qg_tet.h5'

   #tet_filename = 'r40-mri-50-10ps/si_qg_tet.h5'
   #tet_filename = 'r80-mri-50-30ps/si_qg_tet.h5'
   tet_file = h5py.File(tet_filename, 'r')
   kpts = tet_file['kpts_all_crys_coord'][()]
   tet_file.close()

   kpts = np.matmul(kpts, bg)
   BZ_array = np.matmul(BZ_array, bg)

   # Fold kpts into the first BZ
   kpts = ppu.kpts_fold(kpts,bg)
   #
   # 3D array of zeros
   #
   ngrid = 160
   dim=np.zeros(3,dtype=np.int32)
   #labels_in = np.zeros((ngrid,ngrid,ngrid), dtype=np.int32)
   #
   # !!! Because of the folding, the indicies can be < 0 or > ngrid !!!
   #
   labels_in = np.zeros(tuple(dim), dtype=np.int32)
   
   #
   # Convert to integers with the origin in 0,0,0
   #
   shift_int_array = np.array([ngrid/2,ngrid/2,ngrid/2],dtype=np.int32)
   
   connectivity = 6 # only 4,8 (2D) and 26, 18, and 6 (3D) are allowed
   
   #labels_out = labels_out.astype(np.int32)
   
   plot_3D = True
   
   if plot_3D:
      print('Plotting...')
   
      t0 = datetime.datetime.now()
   
      plot_type = 'border'
      #plot_type = 'threshold'
      #plot_type = 'double_threshold'
      
      fig = plt.figure(figsize=(8,8),constrained_layout=True)
      ax = fig.add_subplot(111, projection='3d', proj_type='persp') # 'persp' 'ortho'
      #ax.xaxis.pane.fill = False
      ax.set_axis_off()
      
      #azim, elev = -82.63, 18.83
      azim, elev = 45, 19
      #azim, elev = -82.78104903555192, 8.105518475814904
      
      ax.view_init(azim=azim,elev=elev)
   
      draw_BZ_border(ax,BZ_array,order_dir=np.array([azim,elev]), c='grey', alpha=0.8, lw=1)
      
      incr=1
      #incr=30
      #ax.scatter(kpts[::incr,0],kpts[::incr,1],kpts[::incr,2],color='black',s=0.1)
   
      #draw_axes(ax,c='black',lw=0.4)
   
   #plt.show()
   return fig, ax, kpts, bg

##################################################
def plot_border2(fig, ax, folder):
   
   highsym_dict_Si = {
      'Xc1': {
         'color' : 'blue', 'ls':'-',
         #'color' : 'black',
         'radius': 0.24,
         'coord' : np.array([
                  [-0.5, 0.0,-0.5 ],
                  ])
      },
      
      'Xc2': {
         #'color' : 'black',
         'color' : 'blue', 'ls':'--',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.5, 0.0, 0.5 ],
                  ])
      },
      
      'Xh1': {
         #'color' : 'black', 
         'ls':'-',
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.5, 0.5, 0.0 ],
                  ])
      },
      
      'Xh2': {
         #'color' : 'black',
         'ls':'--',
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [-0.5,-0.5, 0.0 ],
                  ])
      },
      
      'Xh3': {
         #'color' : 'black',
         'ls':'-.',
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.0,-0.5,-0.5 ],
                  ])
      },
      
      'Xh4': {
         #'color' : 'black',
         'ls':'dotted',
         'color' : 'red',
         'radius': 0.24,
         'coord' : np.array([
                  [ 0.0, 0.5, 0.5 ],
                  ])
      },
   }
   
   BZ_array = np.array([[ 0.25, -0.25, -0.5 ],
          [-0.25, -0.5 , -0.75],
          [-0.5 , -0.25, -0.75],
          [-0.25,  0.25, -0.5 ],
          [ 0.25,  0.5 , -0.25],
          [ 0.5 ,  0.25, -0.25],
          [ 0.25, -0.25, -0.5 ],
          [ 0.25, -0.5 , -0.25],
          [-0.25, -0.75, -0.5 ],
          [-0.25, -0.5 , -0.75],
          [-0.25, -0.75, -0.5 ],
          [-0.5 , -0.75, -0.25],
          [-0.25, -0.5 ,  0.25],
          [ 0.25, -0.25,  0.5 ],
          [ 0.5 , -0.25,  0.25],
          [ 0.25, -0.5 , -0.25],
          [ 0.5 , -0.25,  0.25],
          [ 0.75,  0.25,  0.5 ],
          [ 0.75,  0.5 ,  0.25],
          [ 0.5 ,  0.25, -0.25],
          [ 0.75,  0.5 ,  0.25],
          [ 0.5 ,  0.75,  0.25],
          [ 0.25,  0.5 , -0.25],
          [ 0.5 ,  0.75,  0.25],
          [ 0.25,  0.75,  0.5 ],
          [-0.25,  0.5 ,  0.25],
          [-0.5 ,  0.25, -0.25],
          [-0.25,  0.25, -0.5 ],
          [-0.5 ,  0.25, -0.25],
          [-0.75, -0.25, -0.5 ],
          [-0.5 , -0.25, -0.75],
          [-0.75, -0.25, -0.5 ],
          [-0.75, -0.5 , -0.25],
          [-0.5 , -0.75, -0.25],
          [-0.25, -0.5 ,  0.25],
          [-0.5 , -0.25,  0.25],
          [-0.75, -0.5 , -0.25],
          [-0.5 , -0.25,  0.25],
          [-0.25,  0.25,  0.5 ],
          [ 0.25,  0.5 ,  0.75],
          [ 0.5 ,  0.25,  0.75],
          [ 0.25, -0.25,  0.5 ],
          [ 0.5 , -0.25,  0.25],
          [ 0.75,  0.25,  0.5 ],
          [ 0.5 ,  0.25,  0.75],
          [ 0.25,  0.5 ,  0.75],
          [ 0.25,  0.75,  0.5 ],
          [-0.25,  0.5 ,  0.25],
          [-0.25,  0.25,  0.5 ]])
   
   
   #
   # Read k-points of the full BZ with the energy window
   #
   #prefix='gaas'
   prefix='si'
   
   highsym_dict = highsym_dict_Si
   epwan_filename = 'si_epr.h5' 
   epwan_file = h5py.File(epwan_filename, 'r')
   # bg - reciprocal lattice vectors in unit of $2\pi/a$ where $a$ is the lattice constant
   bg = epwan_file['basic_data']['bg'][()]
   alat = epwan_file['basic_data']['alat'][()]
   #bg = bg * 2 * np.pi / alat
   epwan_file.close()
   
   #
   # !!! transform kpts and BZ to cartesian coordinates
   #
   
   
   #tet_filename = f'{folder}/si_tet.h5'
   tet_filename = f'{folder}/si_qg_tet.h5'

   #tet_filename = 'r40-mri-50-10ps/si_qg_tet.h5'
   #tet_filename = 'r80-mri-50-30ps/si_qg_tet.h5'
   tet_file = h5py.File(tet_filename, 'r')
   kpts = tet_file['kpts_all_crys_coord'][()]
   tet_file.close()

   kpts = np.matmul(kpts, bg)
   BZ_array = np.matmul(BZ_array, bg)

   # Fold kpts into the first BZ
   kpts = ppu.kpts_fold(kpts,bg)
   #
   # 3D array of zeros
   #
   ngrid = 160
   dim=np.zeros(3,dtype=np.int32)
   #labels_in = np.zeros((ngrid,ngrid,ngrid), dtype=np.int32)
   #
   # !!! Because of the folding, the indicies can be < 0 or > ngrid !!!
   #
   labels_in = np.zeros(tuple(dim), dtype=np.int32)
   
   #
   # Convert to integers with the origin in 0,0,0
   #
   shift_int_array = np.array([ngrid/2,ngrid/2,ngrid/2],dtype=np.int32)
   
   connectivity = 6 # only 4,8 (2D) and 26, 18, and 6 (3D) are allowed
   
   
   plot_type = 'border'
   #plot_type = 'threshold'
   #plot_type = 'double_threshold'
   
   #ax.xaxis.pane.fill = False
   ax.set_axis_off()
   
   #azim, elev = -82.63, 18.83
   azim, elev = 45, 19
   #azim, elev = -82.78104903555192, 8.105518475814904
   
   ax.view_init(azim=azim,elev=elev)
   
   draw_BZ_border(ax,BZ_array,order_dir=np.array([azim,elev]), c='grey', alpha=0.8, lw=1)
   
   #draw_axes(ax,c='black',lw=0.4)
   
   return fig, ax, kpts, bg

