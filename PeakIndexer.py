# -*- coding: utf-8 -*-
"""
Created on Wed Feb 01 11:00:49 2022

@author: shabalin

Class to work with a list of peaks and *.flt files.
"""

import os
import numpy as np
from datetime import datetime
from ImageD11 import indexing
from ImageD11 import transformer
from Geometry import Geometry
from GvectorEvaluator import GvectorEvaluator
# from SweepProcessor import SweepProcessor

single_separator = "--------------------------------------------------------------\n"
double_separator = "==============================================================\n"

class PeakIndexer:
    
    def __init__(self, directory = None):
        self.directory = None
        self.name = None
        self.gve_file = None
        self.header  = []
        self.peaks = []
        self.log = []
        self.absorbed = []
        self.geometry = None
        self.spot3d_id_reg = 100000
        self.add_to_log('Initialized PeakIndexer object.', True)
        if directory: self.set_attr('directory', directory)
        return
 

    def add_to_log(self, str_to_add, also_print = False):
        self.log.append( str(datetime.now()) + '> ' + str_to_add )
        if also_print: print(str_to_add)
        return

    
    def set_attr(self, attr, value):
        old = getattr(self, attr)
        setattr(self, attr, value)
        new = getattr(self, attr)
        if attr in ['absorbed', 'header', 'peaks']:
            old, new = f'list of {len(old)}', f'list of {len(new)}'
        if attr == 'geometry' and old is not None: old = old.__dict__
        if attr == 'geometry' and new is not None: new = new.__dict__
        self.add_to_log(attr+': '+str(old)+' -> '+str(new))
        return
        
        
    def add_to_attr(self, attr, value):
        old_list = getattr(self, attr)
        if type(old_list) == list: 
            setattr(self, attr, old_list+[value])
            new_list = getattr(self, attr)
            self.add_to_log(attr+': += '+str(new_list[-1]))
        else:
            raise AttributeError('This attribute is not a list!')
        return
    
    
    def print(self, also_log = False):
        print(double_separator+'PeakIndexer object:')
        print('absorbed:' , len(self.absorbed))
        print('directory:', self.directory)
        print('name:'     , self.name )
        print('gve_file:' , self.gve_file)
        [print(f'header {i:2}:', r) for i,r in enumerate(self.header)]
        print('peaks:', len(self.peaks))
        print('spot3d_id_reg:', self.spot3d_id_reg)
        if self.geometry:
            print('geometry:')
            self.geometry.print()
        if also_log:
            print(single_separator + 'Log:')
            for record in self.log: print(record)
        return
    
    
    def load_flt(self, directory = None, flt_file = None):
        if directory: self.set_attr('directory', directory)
        if flt_file: self.set_attr('name', flt_file.replace('.flt', '') )
        self.add_to_log('Reading file:'+self.directory+self.name+'.flt', True)
        if not os.path.isfile(self.directory+self.name+'.flt'): raise FileNotFoundError  
        header, peaks, inds = [], [], []
        with open(self.directory+self.name+'.flt', "r") as f:
            for line in f:
                if line[0] == '#':
                    if 'sc  fc  omega' in line:
                        peak_keys = line[2:-1].split()
                        inds.append(peak_keys.index("Number_of_pixels"))
                        inds.append(peak_keys.index("IMax_s"))
                        inds.append(peak_keys.index("IMax_f"))
                        inds.append(peak_keys.index("Min_s"))
                        inds.append(peak_keys.index("Max_s"))
                        inds.append(peak_keys.index("Min_f"))
                        inds.append(peak_keys.index("Max_f"))
                        inds.append(peak_keys.index("onfirst"))
                        inds.append(peak_keys.index("onlast"))
                        inds.append(peak_keys.index("spot3d_id"))
                    else: header.append(line[:-1])
                else:
                    words = line.split()
                    if len(words) == len(peak_keys):
                        p = [float(v) for v in words]
                        for i in inds: p[i] = int(words[i])
                        peaks.append( dict(zip(peak_keys,p)) )
        f.close()
        self.set_attr('header', header)
        self.set_attr('peaks', peaks)
        print(f'{len(peaks)} peaks loaded.')
        self.add_to_log('File closed!', False)
        return
 

    def save_flt(self, directory = None, flt_file = None, overwrite = False):
        if directory: self.set_attr('directory', directory)
        if flt_file: self.set_attr('name', flt_file.split('.')[0])
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            self.add_to_log('Created directory: '+self.directory, True)
        
        self.add_to_log('Writing file: '+self.directory+self.name+'.flt', True)
        # Check if file exists, if yes then ask for the permission to overwrite
        while os.path.isfile(self.directory+self.name+'.flt'):
            self.add_to_log('File already exist!', True)
            if overwrite:
                self.add_to_log('Overwriting...', True)
                break
            else:
                x = input('Type new name or ! to overwrite, a - abort:')
                if x in ['!']:
                    self.add_to_log('Overwriting...', True)
                    break
                elif x in ['a', 'A']:
                    self.add_to_log('Aborted!', True)
                    return
                else:
                    self.set_attr('name', x.split('.')[0])
        
        f = open(self.directory+self.name+'.flt' ,"w") 
        for line in self.header: f.write(line + '\n')
        if len(self.peaks)>0: f.write('#  ' + '  '.join(self.peaks[0].keys()) + '\n' )
        for p in self.peaks:
            s = [f'{v:d}' if type(v)==type(1) else f'{v:.4f}' for v in list(p.values())]
            f.write('  ' + '  '.join(s) + '\n' )
        f.close()
        self.add_to_log('File closed!', False)
        return
    
    
    def remove_not_in_range(self, int_range):
        peaks_in_range = []
        for p in self.peaks:
            if p['IMax_int'] < int_range[0] or p['IMax_int'] > int_range[1]:
                pass
            elif p['Number_of_pixels'] < 3:
                pass
            else:
                peaks_in_range.append(p)
        self.set_attr('peaks', peaks_in_range)
    
    
    def absorb(self, x, spot3d_id_reg = None):
        if type(spot3d_id_reg) == type(1): self.set_attr('spot3d_id_reg', spot3d_id_reg)
        if self.spot3d_id_reg < 1: raise ValueError(f'spot3d_id_reg must be > 0 (recommended 100000). Provided: {self.spot3d_id_reg}')
        if min([p['spot3d_id'] for p in x.peaks])+self.spot3d_id_reg*len(self.absorbed) < max([p['spot3d_id'] for p in self.peaks])+1:
            raise ValueError(f'spot3d_id_reg = {self.spot3d_id_reg} is too low for assigning unique spot3d_ids!')
            
        self.add_to_log(f'Absorbing PeakIndexer object... shifting its spot3d_ids by {self.spot3d_id_reg*len(self.absorbed)}.')
        self.header.insert(0, '#  Primary PeakIndexer object: ' + self.directory + self.name)
        self.header.append('# Absorbed PeakIndexer object: ' + x.directory + x.name)
        self.header.append(f'# Spot3d_id shift {add_to_id} computed as multiple of spot3d_id_reg = {self.spot3d_id_reg}')
        add_to_header = [line for line in x.header if line not in self.header]
        for line in add_to_header: self.add_to_attr('header', line)
        common = [k for k in self.peaks[0].keys() if k in x.peaks[0].keys()]
        peaks = []
        for p in self.peaks:
            c = {key:value for key,value in p.items() if key in common}
            peaks.append(c)
                        
        for p in x.peaks:
            c = {key:value for key,value in p.items() if key in common}
            c['spot3d_id'] += self.spot3d_id_reg*len(self.absorbed)
            peaks.append(c)
            
        ind = np.argsort([p['spot3d_id'] for p in peaks]) # Sort in asceding spot3d_id 
        self.set_attr('peaks', [peaks[i] for i in ind])
        self.add_to_attr('absorbed', x)
        return


    def run_indexer(self, directory = None, par_file = None, gve_file = None):
        if directory: self.set_attr('directory', directory)
        if par_file : self.geometry.set_attr('par_file' , par_file)
        if gve_file : self.set_attr('gve_file' , gve_file)
        else: self.set_attr('gve_file' , self.name+'.gve')
        self.geometry.save_par(directory = self.directory, par_file = self.name+'.par', overwrite = True)
        self.geometry.save_yml(directory = self.directory, yml_file = self.name+'.yml', overwrite = True)
        self.add_to_log('Running indexer on: '+self.directory+self.name+'.flt', True)
        self.add_to_log('Using par_file: '+self.geometry.directory+self.geometry.par_file, True)
        obj = transformer.transformer()
        obj.loadfiltered( self.directory+self.name+'.flt' )
        obj.loadfileparameters( self.geometry.directory+self.geometry.par_file )
        obj.compute_tth_eta( )
        obj.addcellpeaks( )
        obj.computegv( )
        obj.savegv( self.directory + self.gve_file )
        self.add_to_log('Resulted gvectors saved in: '+self.directory+self.gve_file, True)
        return
    
    
    def generate_GvectorEvaluator(self, directory = None, name = None):
        if not directory: directory = self.directory
        if not name: name = self.name
        if not os.path.exists(directory):
            os.makedirs(directory)
            self.add_to_log('Created directory: '+directory, True)
        GE = GvectorEvaluator(directory = directory)
        GE.set_attr('geometries', [self.geometry])
        GE.load_gve(gve_file = name+'.gve') # merged
        self.add_to_log('Generated GvectorEvaluator in '+directory+' named '+name, True)
        return GE    
#     def import_metadata(self, sweep_processor):
#         SP = sweep_processor
#         self.add_to_log(f'Importing metadata from SweepProcesor: {SP.directory}{SP.name}', True)
#         self.set_attr('directory', SP.directory)
#         self.set_attr('name'     , SP.name)
#         self.set_attr('geometry' , SP.geometry)
#         return