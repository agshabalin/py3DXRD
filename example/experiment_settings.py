# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:15:49 2022
@author: shabalin
"""

import os, sys, subprocess, pickle
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
sys.path.insert(0, '/home/shabalin/')
import py3DXRD
from py3DXRD.p212_tools import load_p212_log, load_p212_fio, parse_p212_command

single_separator = "--------------------------------------------------------------"
double_separator = "=============================================================="

### MATERIAL TABLE. SETTING WORKING DIRECTORY.
path_gen = '/asap3/petra3/gpfs/p21.2/2022/data/11013744/'
materials_table = {'Au'  : {'name': 'Au'  , 'spacegroup': 225, 'symmetry':'F', 'unitcell': [4.0720, 4.0720,  4.0720, 90., 90.,  90.]},
                   'CeO2': {'name': 'CeO2', 'spacegroup': 225, 'symmetry':'F', 'unitcell': [5.4115, 5.4115,  5.4115, 90., 90.,  90.]},
                   'Cu'  : {'name': 'Cu'  , 'spacegroup': 225, 'symmetry':'F', 'unitcell': [3.5942, 3.5942,  3.5942, 90., 90.,  90.]},
                   'LaB6': {'name': 'LaB6', 'spacegroup': 221, 'symmetry':'P', 'unitcell': [4.1568, 4.1568,  4.1568, 90., 90.,  90.]},
                   'MgCa': {'name': 'MgCa', 'spacegroup': 194, 'symmetry':'P', 'unitcell': [3.1980, 3.1980,  5.1900, 90., 90., 120.]},
                   'Nb'  : {'name': 'Nb'  , 'spacegroup': 229, 'symmetry':'I', 'unitcell': [3.3042, 3.3042,  3.3042, 90., 90.,  90.]},
                   'Ni'  : {'name': 'Ni'  , 'spacegroup': 225, 'symmetry':'F', 'unitcell': [3.5100, 3.5100,  3.5100, 90., 90.,  90.]},
                   'Pd'  : {'name': 'Pd'  , 'spacegroup': 225, 'symmetry':'F', 'unitcell': [3.8907, 3.8907,  3.8907, 90., 90.,  90.]},
                   'Ruby': {'name': 'Ruby', 'spacegroup': 167, 'symmetry':'R', 'unitcell': [4.7608, 4.7608, 12.9957, 90., 90., 120.]},
                   'Ti'  : {'name': 'Ti'  , 'spacegroup': 194, 'symmetry':'P', 'unitcell': [2.9505, 2.9505,  4.6826, 90., 90., 120.]}}

### SETTING THE SWEEP PARAMETERS:
def set_p212_sweep(i_slow, i_fast, det_num, default_xyz = [0,0,0], meta_key = None):
    # meta_key is either full path to .fio file or .log file or a string with sweep command.
    x,y,z = None,None,None
    if '.log' in meta_key.split()[-1]:
        log_meta = load_p212_log(meta_key)
        command = parse_p212_command(log_meta['command'])
        if command['slow']:
            if 'x' in command['slow']['motor']:
                x  =  command['slow']['start'] + i_slow*command['slow']['step']
            if 'y' in command['slow']['motor']:
                y  =  command['slow']['start'] + i_slow*command['slow']['step']
            if 'z' in command['slow']['motor']:
                z  =  command['slow']['start'] + i_slow*command['slow']['step']
        if command['fast']:
            if 'x' in command['fast']['motor']:
                x  =  command['fast']['start'] + i_fast*command['fast']['step']
            if 'y' in command['fast']['motor']:
                y  =  command['fast']['start'] + i_fast*command['fast']['step']
            if 'z' in command['fast']['motor']:
                z  =  command['fast']['start'] + i_fast*command['fast']['step']

        x_fn, y_fn, z_fn = default_xyz # As it is in filenames.
        if x: x_fn = x
        if y: y_fn = y
        if z: z_fn = z
        fio_file_path = meta_key.replace('.log', f'_y_{y_fn:.3f}_z_{z_fn:.3f}.fio')
        fio_meta = load_p212_fio(fio_file_path)
        command = parse_p212_command(fio_meta['command'])
    elif '.fio' in meta_key.split()[-1]:
        log_meta = None
        fio_file_path = meta_key
        fio_meta = load_p212_fio(fio_file_path)
        command = parse_p212_command(fio_meta['command'])
#         print('COMMAND:', command)
#         print('fio meta command:', fio_meta['command'])
    else:
        log_meta = None
        fio_meta = None
        command = parse_p212_command(meta_key)
        
    if command['slow']:
        if 'x' in command['slow']['motor']:
            x  =  command['slow']['start'] + i_slow*command['slow']['step']
        if 'y' in command['slow']['motor']:
            y  =  command['slow']['start'] + i_slow*command['slow']['step']
        if 'z' in command['slow']['motor']:
            z  =  command['slow']['start'] + i_slow*command['slow']['step']
    if command['fast']:
        if 'x' in command['fast']['motor']:
            x  =  command['fast']['start'] + i_fast*command['fast']['step']
        if 'y' in command['fast']['motor']:
            y  =  command['fast']['start'] + i_fast*command['fast']['step']
        if 'z' in command['fast']['motor']:
            z  =  command['fast']['start'] + i_fast*command['fast']['step']
    
    if not x:
        try:    x = fio_meta['idtx2'] # Define here which motor should be used as x.
        except: x = default_xyz[0]
    if not y:
        try:    y = fio_meta['idty1'] # Define here which motor should be used as y.
        except: y = default_xyz[1]
    if not z:
        try:    z = fio_meta['idtz2'] # Define here which motor should be used as z.
        except: z = default_xyz[2]    
    
    path_in_raw   = command['directory'].split('raw')[1].replace('%D',f'{det_num}')
    if path_in_raw[-1] != '/': path_in_raw+='/'
    raw_data_path = path_gen + 'raw'       + path_in_raw
    sweep_path    = path_gen + 'processed' + path_in_raw
    det_folders =  [f'{det_num}', f'varex{det_num}', f'perk{det_num}']
    if sweep_path.split('/')[-2].replace('_', '').lower() in det_folders:
        sweep_path = sweep_path.replace('/'+sweep_path.split('/')[-2]+'/','/')
    
    ## Sweep initialization:
    print(single_separator)
    SP = py3DXRD.SweepProcessor(directory = sweep_path, name = f's{i_slow:03d}_f{i_fast:03d}_d{det_num}')
    SP.log_meta             = log_meta
    SP.fio_meta             = fio_meta
    SP.sweep['omega_start'] = command['sweep']['start'] + command['sweep']['step'] - 180 # -180 because Grainspotter uses omega range [-180;+180]
    SP.sweep['omega_step' ] = command['sweep']['step']
    SP.sweep['directory'  ] = raw_data_path
    SP.sweep['stem'       ] = command['file_stem']
    SP.sweep['ndigits'    ] = 'auto'
    SP.sweep['ext'        ] = command['file_ext']
    SP.chunk['frames'     ] = list(range(command['sweep']['points']))[1:-2] # Drop out the first frame and the last two frames
    SP.position             = [x,y,z]
    SP.processing['options']= None # [('flip', 'r270')]
    SP.spline_file          = None # Perfect detector or spline needed
    
#     print(f'Sweep directory: {SP.directory}')
#     print(f'Sweep name: {SP.name}')
#     if SP.log_meta:
#         print('Log file:', SP.log_meta['directory'] + SP.log_meta['file'])
#         print('Command :', SP.log_meta['command'])
#     if SP.fio_meta:
#         print('Fio file:', SP.fio_meta['directory'] + SP.fio_meta['file'])
#         print('Command :', SP.fio_meta['command'])
#     else:
#         print('Command string:', meta_key)
#     print(f'Position inidices: slow = {i_slow}, fast = {i_fast}')
#     print('Position (x,y,z):'  , SP.position)
#     print('Data path template:', SP.sweep['directory'], SP.sweep['stem'], SP.sweep['ext'])
    
    return SP

### SETTING THE GRAINSPOTTER PARAMETERS;
def set_grainspotter(material=None, domega = None):
    GS = py3DXRD.GrainSpotter(directory = path_gen + 'processed/')
    if material: GS.set_attr('spacegroup'   , material['spacegroup'])
    if domega: GS.set_attr('domega'       , domega)
    GS.set_attr('tth_ranges'   , [   [ 8.0, 17.0] ] ) # 12.7]])
    GS.set_attr('ds_ranges'    , [   [0.45, 1.2] ] ) # [0.5, 1.0]]) # GV.ds_ranges)
    GS.set_attr('eta_ranges'   , [   [ 5, 85], [275, 355]] ) # GV.eta_ranges)
    GS.set_attr('omega_ranges' , [   [-179.5,  179.5]] ) # GV.omega_ranges)
    GS.set_attr('cuts'         , [ 12, 0.6,  0.6] )
    GS.set_attr('uncertainties', [0.2, 1.5,  1.0] ) # [sigma_tth sigma_eta sigma_omega] in degrees
    GS.set_attr('nsigmas'      , 1)
    GS.set_attr('eulerstep'    , 6)
    GS.set_attr('Nhkls_in_indexing', None)
    GS.set_attr('random', 10000)
    GS.set_attr('positionfit', True)
    return GS

### SETTING THE POLYXSIM PARAMETERS:
def set_polyxsim(grainspotter, material = None):
    if not material: material = {'name': ''}
    PS = py3DXRD.PolySim(directory = grainspotter.directory+'sim_'+material['name']+'/')
    PS.set_attr('inp_file', grainspotter.log_file.strip('.log'))
    PS.set_attr('beamflux', 1e12)
    PS.set_attr('beampol_factor', 1)
    PS.set_attr('beampol_direct', 0)
    PS.set_attr('direc' , grainspotter.directory+'sim_'+material['name']+'/')
    PS.set_attr('stem'  , grainspotter.log_file.replace('.log', '_sim'))
    PS.set_attr('grains', [])
    PS.set_attr('omega_start', grainspotter.omega_ranges[0][0])
    PS.set_attr('omega_step' , abs(grainspotter.domega))
    PS.set_attr('omega_end'  , grainspotter.omega_ranges[-1][1])
    PS.set_attr('theta_min'  , grainspotter.tth_ranges[0][0]/2)
    PS.set_attr('theta_max'  , grainspotter.tth_ranges[-1][1]/2)
    PS.set_attr('no_grains'  , 1)
    PS.set_attr('gen_U'   , 0)
    PS.set_attr('gen_pos' , [0, 0])
    PS.set_attr('gen_eps' , [1, 0, 0 ,0, 0])
    PS.set_attr('gen_size', [0.0, 0.0, 0.0 ,0.0])
    PS.set_attr('make_image', 0)
    PS.set_attr('output', ['.tif', '.par', '.gve'])
    PS.set_attr('bg' , 0)
    PS.set_attr('psf', 0.7)
    PS.set_attr('peakshape', [1, 4, 0.5])
    return PS

print(single_separator+'\nSETTINGS LOADED!')
# print(f'Path:', path_gen)
# print('Material:', material)