# -*- coding: utf-8 -*-
"""
Created on Wed Feb 01 11:00:49 2022

@author: shabalin

Class to work with a list of peaks and *.flt files.
"""

import os
import numpy as np
from datetime import datetime

single_separator = "--------------------------------------------------------------\n"
double_separator = "==============================================================\n"

class FltFile:
    
    def __init__(self, directory = None):
        self.log = []
        self.absorbed = []
        self.directory = None
        self.flt_file  = None
        self.header  = []
        self.peaks = []
        self.add_to_log('Created FltFile object.', True)
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
        self.add_to_log(attr+': '+str(old)+' -> '+str(new))
        return
        
        
    def add_to_attr(self, attr, value):
        old_list = getattr(self, attr)
        setattr(self, attr, old_list+[value])
        new_list = getattr(self, attr)
        self.add_to_log(attr+': += '+str(new_list[-1]))
        return
    
    
    def print(self, also_log = False):
        print(double_separator+'FltFile object:')
        print('absorbed:', len(self.absorbed))
        print('directory:', self.directory)
        print('flt_file:' , self.flt_file )
        [print(f'header {i:2}:', r) for i,r in enumerate(self.header)]
        print('peaks:', len(self.peaks))
        if also_log:
            print(single_separator + 'Log:')
            for record in self.log: print(record)
        return
    
    
    def load_flt(self, directory = None, flt_file = None):
        if directory: self.set_attr('directory', directory)
        if flt_file: self.set_attr('flt_file', flt_file)
        self.add_to_log(f'Reading file: {self.directory+self.flt_file}', True)
        if not os.path.isfile(self.directory+self.flt_file): raise FileNotFoundError  
        header, peaks, inds = [], [], []
        with open(self.directory+self.flt_file, "r") as f:
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
        if flt_file: self.set_attr('flt_file', flt_file)
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            self.add_to_log('Created directory: '+self.directory, True)
        
        self.add_to_log(f'Writing file: {self.directory+self.flt_file}', True)
        # Check if file exists, if yes then ask for the permission to overwrite
        while os.path.isfile(self.directory+self.flt_file):
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
                    self.set_attr('flt_file', x)
        
        f = open(self.directory+self.flt_file ,"w") 
        for line in self.header: f.write(line + '\n')
        f.write('#  ' + '  '.join(self.peaks[0].keys()) + '\n' )
        for p in self.peaks:
            s = [f'{v:d}' if type(v)==type(1) else f'{v:.4f}' for v in list(p.values())]
            f.write('  ' + '  '.join(s) + '\n' )
        f.close()
        self.add_to_log('File closed!', False)
        return
    
    
    def absorb(self, x, add_to_id = 100000):
        if min([p['spot3d_id'] for p in x.peaks])+add_to_id < max([p['spot3d_id'] for p in self.peaks])+1:
            raise ValueError(f'ind_to_add={add_to_id} is too low for assigning unique spot3d_ids!')
            
        self.add_to_log('Absorbing FltFile object...')
        self.header.insert(0, '#  Primary flt_file: ' + self.directory + self.flt_file)
        self.header.append('# Absorbed flt_file: ' + x.directory + x.flt_file)
        add_to_header = [line for line in x.header if line not in self.header]
        for line in add_to_header: self.add_to_attr('header', line)
        common = [k for k in self.peaks[0].keys() if k in x.peaks[0].keys()]
        peaks = []
        for p in self.peaks:
            c = {key:value for key,value in p.items() if key in common}
            peaks.append(c)
                        
        for p in x.peaks:
            c = {key:value for key,value in p.items() if key in common}
            c['spot3d_id'] += add_to_id #100000*(1+int(max_id/100000))
            peaks.append(c)
            
        ind = np.argsort([p['spot3d_id'] for p in peaks]) # Sort in asceding spot3d_id 
        self.set_attr('peaks', [peaks[i] for i in ind])
        self.add_to_attr('absorbed', x)
        return
    
#### TESTS    ###############
# G = FltFile(directory='/asap3/petra3/gpfs/p21.2/2021/data/11008399/processed/CS_1/load1/z_00/')
# G.load_flt(flt_file = 'V4_peaks_t32000.flt')
# G2 = FltFile(directory='/asap3/petra3/gpfs/p21.2/2021/data/11008399/processed/CS_1/load1/z_01/')
# G2.load_flt(flt_file = 'V4_peaks_merged.flt')
# G.absorb(G2)
# G.save_flt(flt_file = 'V4_peaks_t32000_test.flt')
# G.print(True)