#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 10:41:54 2019

@author: nikonalesheo
"""
import glob
import sys
import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from openfits import fits_open
from openfits import geomspace
from openfits import ax_set
from astropy.stats import mad_std
from mpl_toolkits.axes_grid1 import make_axes_locatable

from difmap_cc import cc_fluxes 

#from matplotlib import rc
#rc('text', usetex=True)


no_release = True

input_pars = sys.argv
input_pars_len = len(input_pars)

if no_release:
    input_pars = ['test', '0923+392']
    input_pars_len = len(input_pars)
    input_pars_EVPA = {'0851+202':153, '1308+326': 70, '0923+392': 152}

#current_dir = os.getcwd()

obj_n = glob.glob("*_i.fits")
obj_n_len = len(obj_n)
objects = []

for i in range(obj_n_len):
    obj_name = obj_n[i][:-7]
    if input_pars[1] != obj_name:
        objects.append(obj_name)

objects.append(input_pars[1])
objects.sort()

objects_len = len(objects)

print(objects)

    
if input_pars[1]+'_i.fits' in obj_n:
    print('Source to callibrate: '+input_pars[1])
else:
    sys.exit("ERROR: enter correct source for EVPA callibration. Check //$python3 pol_maps.py// help for details")
    
if input_pars == 1:
    sys.exit("ERROR: Enter source to callibrate. Check //$python3 pol_maps.py// help for details")

if input_pars[1] == 'help':
    sys.exit("Three parameters needed to run script: filenames of stokes i, u, q" +
             " maps in partucular order. Same mapsize, beam etc." )
    
    
#print("Stokes I: "+glob.glob("*_i.fits"), "Stokes U: " + glob.glob("*_u.fits"), "Stokes Q: " + 
#      glob.glob("*_q.fits"))
#print(filenames[1], type(filenames[0])) Аргументы, которые при запуске использ
# 0 элемент списка - первый аргумент - название скрипта, 1 - первый параметр


def delta_EVPA_finder(obj_name,mode, **kwargs):
    """
    This function needed to calculate EVPA difference between 
    
    **kwargs
    contours_min = 1; default = 1 sigma (in sigma calculated via mad_std from astropy)
    """    
        
    i_map, i_head = fits_open(obj_name+'_i.fits') #opening fits with my lib
    u_map, u_head = fits_open(obj_name+'_u.fits')
    q_map, q_head = fits_open(obj_name+'_q.fits')
    
    obj = i_head['OBJECT'] #just name from header
    map_size = np.shape(i_map)
    im_scale = abs(i_head['CDELT1']*60*60*1000) #for ticks in plots
    
    if mode == 'cal_target':
        if no_release == True: #this is for auto enternig the params in script
            file_save = obj
        else:
            file_save = input("Enter polarization map name (.pdf): ")
    
    i_sigma = mad_std(i_map) #the right way to caculate std
    u_sigma = mad_std(u_map)
    q_sigma = mad_std(q_map)
    
    print(i_sigma)
    
    nongauss = 1.36 #quantille: the stokes u and q maps have non-gauss type of noise. 
    
    i_map = ma.masked_less(i_map, 5*i_sigma) # cutt-off 5sigma level of full intensity
    
    I_pol = np.sqrt(u_map**2+q_map**2) # polaraized intensity
    
    I_pol_mask = ma.getmask(ma.masked_less(I_pol, nongauss*5*np.sqrt(u_sigma**2
                                                                     +q_sigma**2)))
    I_pol = ma.array(I_pol, mask = I_pol_mask)
    u_map = ma.array(u_map, mask = I_pol_mask)
    q_map = ma.array(q_map, mask = I_pol_mask)
    
    pol_angle_map = np.arctan(u_map/q_map)  # + n*pi for doing EVPA 
    
       
    mask = ma.getmask(pol_angle_map) + i_map.mask
    
    q_sum = np.sum(cc_fluxes(obj[:4], 'q')) # Loading q and u stokes clean componets
    u_sum = np.sum(cc_fluxes(obj[:4], 'u')) # fluxes and sum
        

    if no_release == True:
        EVPA = input_pars_EVPA[obj]
        print('NO RELEASE then default EVPA for '+obj+
              str(round(i_head['CRVAL3']/1000000000,1))+ ' GHz '
              +' =',EVPA)
    else:
        EVPA = float(input("Enter EVPA for "+obj+ ' ' +
                       str(round(i_head['CRVAL3']/1000000000,1))+ ' GHz '+
                       " callibrator: "))
      
    print('Total Flux Q = ',q_sum,' Jy, Total Flux U = ', u_sum, ' Jy')
    
    EVPA_before = np.arctan(u_sum/q_sum)/2 # It`s important to caculate EVPA
    # from clean components fluxes to avoid non-real objects on the maps
    
    delta_EVPA = EVPA_before - np.deg2rad(EVPA)
    print('DELTA EVPA = ', np.rad2deg(delta_EVPA))

    if mode == 'cal_target':

        delta_EVPA = kwargs['delta_EVPA'] #loading the result of caculation of
        #median of delta EVPAs for all calibrators 
        
        EVPA_after = EVPA_before - delta_EVPA
        
        vec_scale_input = 5000 #parameters to draw polarization lines
        step_x = 4
        step_y = 4
        vec_scale = vec_scale_input*im_scale
        
        fig, ax = plt.subplots(figsize=(10,10)) 
        
        left, width = .1, .5 #for text location (left, bottom, right, top)
        bottom, height = .1, .5
        right = left + width
        top = bottom + height
        
        ax.text(left, bottom, '$EVPA_{before}$ = '+
                '{:.0f}'.format(np.rad2deg(EVPA_before))+ '$^{\circ}$'+
                '\n$EVPA_{true (NRAO+UMRAO data)}$ = '+
                str(input_pars_EVPA[obj])+ '$^{\circ}$'
                +
                ' \n$\\Delta EVPA$ = '+
                np.array2string(np.rad2deg(delta_EVPA_tot).astype(int))+ '$^{\circ}$'+
                ' \n$EVPA_{after}$ = ' + 
                str(int(np.rad2deg(EVPA_after)))+ '$^{\circ}$'
                +' \n$P_{polarizaion}$'+'={:.2f}'.format(np.sqrt(q_sum**2+u_sum**2))
                + ' Jy', 
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes)
    
        ax.set_aspect('equal')
        title = 'Clean I contours. Colour P map. Linear polarization map.\nOBJ:' + obj +' '+str(round(
            i_head['CRVAL3']/1000000000,1))+ ' GHz '
    
        ax_set(ax, map_size[0], im_scale, title = title, 
               xlabel = 'realtive RA, mas', ylabel = 'relative DEC, mas')
        
        ### countours drawing
        if 'contours_min' in kwargs:
            contours_min = kwargs['contours_min'] * i_sigma
        else:
            contours_min =  i_sigma
        levels = geomspace(contours_min ,1, 'true')
        ax.contour(i_map, levels, colors = 'black', linewidths = 0.5)
        ###
        
        pos = ax.imshow(I_pol)
        
        ### colorbar
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", "5%", pad="3%")
        cbar = fig.colorbar(pos, cax = cax)
        cbar.ax.set_ylabel('Jy/beam')
        ###
              
        EVPA_map = pol_angle_map - delta_EVPA # correcting EVPA in each pixel
        
        ### drawing lines
        count_y = 0
        for y in range(map_size[1]):
            
            if count_y >= step_y:
                count_x = 0
                for x in range(map_size[0]):
        
                    if mask[y,x] == False:
                        if count_x >= step_x:
                            
                            vec_len = I_pol[y,x]*vec_scale
                            vec_len = 5
                            
                            x_len = vec_len*np.sin(EVPA_map[y,x])/2
                            y_len = vec_len*np.cos(EVPA_map[y,x])/2
                            
                            x_l = [x-x_len, x+x_len] 
                            y_l = [y-y_len, y+y_len]
                            
                            pol_line = Line2D(x_l,y_l, color = 'black', linewidth=0.8)
                            
                            ax.add_line(pol_line)

                            count_x = 0
                        
                        count_x = count_x + 1 
                    else:
                        count_x = count_x + 1 
                count_y = 0
            count_y = count_y + 1

        x = int(ax.get_xlim()[1])#just for zoom
        
        ax.set_xlim(0.25*x, 0.75*x)
        ax.set_ylim(0.25*x, 0.75*x)
        
        fig.savefig(file_save+".pdf", dpi = 300)
    
    return delta_EVPA, EVPA_before




d_EVPAs = np.empty(objects_len)
EVPA_means = np.empty(objects_len)

for i in range(objects_len):
    d_EVPAs[i], EVPA_means[i] = delta_EVPA_finder(objects[i], 'calibrator')

delta_EVPA_tot = np.median(d_EVPAs)

print('Total EVPA correction: ', delta_EVPA_tot)

for i in range(len(objects)):
    delta_EVPA_finder(objects[i], 'cal_target', delta_EVPA = delta_EVPA_tot)