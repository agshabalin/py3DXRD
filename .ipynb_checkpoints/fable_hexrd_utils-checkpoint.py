# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 16:52:49 2021

@author: shabalin

Utils to work with fable and hexrd functions.
"""

import sys, os
import numpy as np
import yaml, subprocess
#import cbftiffmxrdfix

def run_peaksearch(par_file=None):
    """ Wrapper for the ImageD11 peaksearch.py script"""             
    with open(par_file) as f:
        pars = yaml.safe_load(f)
    if pars['stem_out'] == None:
        pars['stem_out'] = ''
    first_im = int(pars['first_image'])
    last_im = int(pars['first_image']) + int(pars['nbr_images']) - 1
    ndigits = pars['ndigits']
    path_inp = os.path.join(pars['image_path'],pars['image_stem'])
    path_out = os.path.join(pars['output_dir'], pars['stem_out']+pars['det_code']+'_peaks')
    # construct the command for peaksearch.py
    command = ('peaksearch.py -n {} -F {} -f {:d} -l {:d} -o {} -d {} -p Y --ndigits {:d} -S {:.3f} -T {:.3f} '.format(
        path_inp,pars['filetype'],first_im,last_im,path_out,
        pars['dark_image'],ndigits,pars['omegastep'], pars['startomega']
        ))
    # Adds threshold values to command
    for t in pars['thresholds']:
        command += '-t {:d} '.format(t)
    # Adds keyword args
    if 'kwargs' in pars:
        command += '{} '.format(pars['kwargs'])
    # modify command for lunarc
    if 'lunarc' in pars:
        command = lunarc_path + command
    print('Running peaksearch with the following command:')
    print(command)
    try:
        subprocess.call(command, shell=True)
    except AttributeError as a:
        print('peaksearch.py ended with error. It seems to work nonetheless.', a)
    
    del pars, first_im, last_im, ndigits, path_inp, path_out, command
    return
      
def merge_peaks(par_file, config_file):
    # Wrapper for ImageD11 merge_flt.py
    if (par_file is None):
        raise ValueError('Must supply par_file to run_peaksearcher')
    with open(par_file) as f:
        pars = yaml.safe_load(f)
    if pars['stem_out'] == None:
        pars['stem_out'] = ''
    if 'merged_name' in pars:
        file_out = os.path.join(pars['output_dir'],pars['stem_out']+pars['merged_name'])
    else:
        file_out = os.path.join(pars['output_dir'],pars['stem_out']+pars['det_code']+'_peaks_merged.flt')
    inp = os.path.join(pars['output_dir'], pars['stem_out']+pars['det_code']+'_peaks')
    print('Merging flt files matching {}'.format(inp))
    if not config_file:
        config_file = 'junk'
    command = 'merge_flt.py {} {} {} {:d} '.format(config_file,inp,file_out,pars['pixel_tol']) + ('{:d} '*len(pars['thresholds'])).format(*pars['thresholds'])
    # modify command for lunarc
    if 'lunarc' in pars:
        command = lunarc_path + command
    print(command)
    subprocess.call(command, shell=True)
    
    del pars, file_out, inp, command
    return

def hexrd_to_fable(path_to_hexrd_yml, path_to_fable_par, det=1, mat='Nb'):
    detname = 'detector_{:d}'.format(det)
    if mat=='ruby':
        cell_params = { "a": 4.7608, "b": 4.7608, "c": 12.99568, "alpha": 90.0, "beta": 90.0, "gamma": 120.0, "lattice": 'R'}
    elif mat=='Nb':
        cell_params = { "a": 3.3042, "b": 3.3042, "c": 3.3042, "alpha": 90.0, "beta": 90.0, "gamma": 90.0, "lattice": 'I'}
    elif mat=='CeO2':
        cell_params = { "a": 5.41153, "b": 5.41153, "c": 5.41153, "alpha": 90.0, "beta": 90.0, "gamma": 90.0, "lattice": 'F'}
    elif mat=='Ti':
        cell_params = { "a": 2.9505, "b": 2.9505, "c": 4.6826, "alpha": 90.0, "beta": 90.0, "gamma": 120.0, "lattice": 'P'}
    else:
        print('ERROR! Incorrect material!')
    with open(path_to_hexrd_yml) as f:
        pars = yaml.safe_load(f)
    wavelength = 12.39842/pars['beam']['energy']
    translation = pars['detectors'][detname]['transform']['translation']
    tilt = pars['detectors'][detname]['transform']['tilt']
    frame_size = [pars['detectors'][detname]['pixels']['columns'], pars['detectors'][detname]['pixels']['rows']]
    pix_size = pars['detectors'][detname]['pixels']['size']
    if os.path.exists(path_to_fable_par):
        if input('File %s already exist! Overwrite it? (y/n):' % path_to_fable_par) != 'y':
            print('Aborted!')
            return
        else:
            pass
    else:
        pass
    f = open(path_to_fable_par,'w') 
    f.write( 'cell__a {}'.format(cell_params['a']) )
    f.write( '\ncell__b {}'.format(cell_params['b']) )
    f.write( '\ncell__c {}'.format(cell_params['c']) )
    f.write( '\ncell_alpha {}'.format(cell_params['alpha']) )
    f.write( '\ncell_beta {}'.format(cell_params['beta']) )
    f.write( '\ncell_gamma {}'.format(cell_params['gamma']) )
    f.write( '\ncell_lattice_[P,A,B,C,I,F,R] {}'.format(cell_params['lattice']) )
    f.write( '\nchi {}'.format(0.0) )
    f.write( '\ndistance {}'.format((-translation[2]*1000)) )
    f.write( '\nfit_tolerance {}'.format(0.5) )
    f.write( '\nmin_bin_prob {}'.format(1e-05) )
    f.write( '\nno_bins {}'.format(10000) )
    f.write( '\no11 {}'.format(0) )
    f.write( '\no12 {}'.format(-1) )
    f.write( '\no21 {}'.format(1) )
    f.write( '\no22 {}'.format(0) )
    f.write( '\nomegasign {}'.format(1.0) )
    f.write( '\nt_x {}'.format(0) )
    f.write( '\nt_y {}'.format(0) )
    f.write( '\nt_z {}'.format(0) )
    f.write( '\ntilt_x {}'.format(tilt[2]) )
    f.write( '\ntilt_y {}'.format(tilt[1]) ) # -?
    f.write( '\ntilt_z {}'.format(tilt[0]) )
    f.write('\nwavelength {:0.6f}'.format(wavelength) )
    f.write( '\nwedge {}'.format(0.0) )
    f.write( '\nweight_hist_intensities {}'.format(0) )
    f.write( '\ny_center {}'.format((translation[1]/pix_size[1] + frame_size[1]/2)) )
    f.write( '\ny_size {}'.format((pix_size[1]*1000)) )
    f.write( '\nz_center {}'.format((translation[0]/pix_size[0] + frame_size[0]/2)) )
    f.write( '\nz_size {}'.format((pix_size[0]*1000)) )
    f.close()
    
    del detname, cell_params, pars, wavelength, translation, tilt, frame_size, pix_size
    return

def fable_to_hexrd(path_to_fable_par, path_to_hexrd_yml):
    y_frm_size = 2880
    z_frm_size = 2880
    with open(path_to_fable_par) as f:
        for line in f:
            if ('distance' in line):
                dist = float(line.split()[1])/1000
            elif ('tilt_x' in line):
                tilt_1 = float(line.split()[1])
            elif ('tilt_y' in line):
                tilt_2 = float(line.split()[1])
            elif ('tilt_z' in line):
                tilt_3 = float(line.split()[1])
            elif ('wavelength' in line):
                wavelength = float(line.split()[1])
            elif ('y_center' in line):
                y_cen = float(line.split()[1])
            elif ('y_size' in line):
                y_pix_size = float(line.split()[1])/1000
            elif ('z_center' in line):
                z_cen = float(line.split()[1])
            elif ('z_size' in line):
                z_pix_size = float(line.split()[1])/1000
    f.close()
    pars = {'beam':
            {'energy': 12.39842/wavelength, 'vector': {'azimuth': 90.0, 'polar_angle': 90.0}},
            'detectors':
            {'detector_1':
             {'buffer': None,
              'pixels': {'columns': y_frm_size, 'rows': z_frm_size, 'size': [z_pix_size, y_pix_size]},
              'saturation_level': 14000.0,
              'transform': {'tilt': [tilt_1, tilt_2, tilt_3], 'translation': [(z_cen-z_frm_size/2)*z_pix_size, (y_cen-y_frm_size/2)*y_pix_size, -dist]}}},
             'id': 'instrument',
             'oscillation_stage': {'chi': 0.0, 'translation': [0.0, 0.0, 0.0]}}
    if os.path.exists(path_to_hexrd_yml):
        if input('File %s already exist! Overwrite it? (y/n):' % path_to_hexrd_yml) != 'y':
            print('Aborted!')
            return
        else:
            pass
    else:
        pass
    with open(path_to_hexrd_yml, 'w') as f:
        yaml.dump(pars, f)
        
    del y_frm_size, z_frm_size, pars, dist, tilt_1, tilt_2, tilt_3, wavelength, y_cen, y_pix_size, z_cen, z_pix_size
    return