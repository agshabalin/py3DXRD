# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 11:00:00 2022

@author: shabalin

Class to work with a grain object. Loads .log  files.
"""

import sys, os, subprocess, pdb, re
import numpy as np
from numpy import float32
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

single_separator = "--------------------------------------------------------------\n"
double_separator = "==============================================================\n"

class Grain:
    
    def __init__(self, directory=None, grain_id=None):
        self.log = []
#         self.absorbed =[]
        self.directory = None
        self.grain_id = None
        self.log_file = None
        self.gff_file = None
        self.spacegroup = None
        self.unitcell = []
        self.u = None
        self.ubi = None
        self.eps = None
        self.size = None
        self.mean_IA = None
        self.position = None
        self.pos_chisq = None
        self.r = None
        self.phi = None
        self.quaternion = None
        self.summary = None
        self.gvectors_report = []
        self.measured_gvectors = []
        self.expected_gvectors = []
        self.add_to_log('Created Grain object.', True)
        if directory: self.set_attr('directory', directory)
        if grain_id : self.set_attr('grain_id' , grain_id)
        return
 

    def add_to_log(self, str_to_add, also_print = False):
        self.log.append( str(datetime.now()) + '> ' + str_to_add )
        if also_print: print(str_to_add)
        return
    
    
    def set_attr(self, attr, value):
        old = getattr(self, attr)
        setattr(self, attr, value)
        new = getattr(self, attr)
        if attr in ['gvectors_report', 'measured_gvectors', 'expected_gvectors']:
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
        print(double_separator+'Grain object:')
        print('directory:' , self.directory )
        print('grain_id:'  , self.grain_id  )
        print('log_file:'  , self.log_file  )
        print('gff_file:'  , self.gff_file  )
        print('spacegroup:', self.spacegroup)
        print('unitcell:'  , self.unitcell  )
        print('u:'         , self.u         )
        print('ubi:'       , self.ubi       )
        print('eps:'       , self.ubi       )
        print('size:'      , self.size      )
        print('mean_IA:'   , self.mean_IA   )
        print('position:'  , self.position  )
        print('pos_chisq:' , self.pos_chisq )
        print('r:'         , self.r         )
        print('phi:'       , self.phi       )
        print('quaternion:', self.quaternion)
        print('summary:'   , self.summary   )
        print('gvectors_report:'  , len(self.gvectors_report  ))
        print('measured_gvectors:', len(self.measured_gvectors))
        print('expected_gvectors:', len(self.expected_gvectors))
#         print('absorbed:', len(self.absorbed))
        if also_log:
            print(single_separator + 'Log:')
            for record in self.log: print(record)
        return

    
    def load_log(directory = None, log_file = None):
        print(double_separator + 'Reading file: ' + directory+log_file)
        header = []
        f = open(directory+log_file,"r")
        line = f.readline()
        while line and '#  gvector_id' not in line:
            if 'Found' in line and 'grains' in line: # the very 1st line in file
                n_grains = int(line.split()[1])
            elif 'Syntax' in line or line == '.':
                pass
            else:
                header.append(line)
            line = f.readline()
        if '#  gvector_id' in line:
            titles = line
        else:
            raise ValueError('Title line not found!') 
        list_of_grains = []
        line = f.readline()
        while line:
            if 'Grain ' in line:
                num_grain = line.split()[1].replace(',','')
                g = Grain(directory, int(num_grain))
                g.set_attr('log_file', log_file)

                line = f.readline()
                keys = header[0].replace('#', '').split()
                values = [int(x) for x in line.split()]
                g.set_attr('summary', dict(zip( keys, values )) )

                line = f.readline()
                values = [float(x) for x in line.split()]
                g.set_attr('mean_IA'  , values[0]  )
                g.set_attr('position' , values[1:4])
                g.set_attr('pos_chisq', values[4]  )

                line = f.readline()
                raw_0 = [float(x) for x in line.split()]
                line = f.readline()
                raw_1 = [float(x) for x in line.split()]
                line = f.readline()
                raw_2 = [float(x) for x in line.split()]
                g.set_attr('u', np.array([raw_0,raw_1,raw_2]) )

                line = f.readline()

                line = f.readline()
                raw_0 = [float(x) for x in line.split()]
                line = f.readline()
                raw_1 = [float(x) for x in line.split()]
                line = f.readline()
                raw_2 = [float(x) for x in line.split()]
                g.set_attr('ubi', np.array([raw_0,raw_1,raw_2]) )

                line = f.readline()
                line = f.readline()
                g.set_attr('r', [float(x) for x in line.split()])

                line = f.readline()
                line = f.readline()
                g.set_attr('phi', [float(x) for x in line.split()])

                line = f.readline()
                line = f.readline()
                g.set_attr('quaternion', [float(x) for x in line.split()])               

                line = f.readline()
                keys = titles.replace('#', '').split()
                gvectors_report = []
                while True:
                    line = f.readline()
                    values = [float(x) if '.' in x else int(x) for x in line.split()]
                    if len(keys) == len(values)-1:
                        gvectors_report.append( dict(zip(keys,values[1:])) )
                    else:
                        break
                g.set_attr('gvectors_report', gvectors_report)
                list_of_grains.append(g)
            else:
                line = f.readline()
        f.close()
        print(f'{len(list_of_grains)} grains loaded.\n')
        return list_of_grains
    

    def load_gff(directory = None, gff_file = None):
        print(double_separator + 'Reading file: ' + directory+gff_file)
        list_of_grains = []
        titles = "grain_id mean_IA chisq x y z U11 U12 U13 U21 U22 U23 U31 U32 U33 UBI11 UBI12 UBI13 UBI21 UBI22 UBI23 UBI31 UBI32 UBI33".split()
        with open(directory+gff_file, "r") as f:
            for line in f:
                if line[0] == '#' and 'grain_id' in line:
                    titles = line[1:-1].split()
                elif len(line.split()) == len(titles):
                    x = dict(zip(titles,line.split()))
                    g = Grain(directory, int(x['grain_id']))
                    g.set_attr('gff_file' , gff_file)
                    g.set_attr('mean_IA'  , float(x['mean_IA']) )
                    g.set_attr('pos_chisq', float(x['mean_IA']) )
                    g.set_attr('position' , [float(x['x']), float(x['y']), float(x['z'])] )
                    raw_0 = [float(v) for v in [x['U11'], x['U12'], x['U13']] ]
                    raw_1 = [float(v) for v in [x['U21'], x['U22'], x['U23']] ]
                    raw_2 = [float(v) for v in [x['U31'], x['U32'], x['U33']] ]
                    g.set_attr('u', np.array([raw_0,raw_1,raw_2]) )
                    raw_0 = [float(v) for v in [x['UBI11'], x['UBI12'], x['UBI13']] ]
                    raw_1 = [float(v) for v in [x['UBI21'], x['UBI22'], x['UBI23']] ]
                    raw_2 = [float(v) for v in [x['UBI31'], x['UBI32'], x['UBI33']] ]
                    g.set_attr('ubi', np.array([raw_0,raw_1,raw_2]) )
                    list_of_grains.append(g)
        f.close()                             
        print(f'{len(list_of_grains)} grains loaded.\n')
        return list_of_grains
    
    
    def save_gff(directory = None, gff_file = None, list_of_grains = None, overwrite = False):
        print(double_separator + 'Writing file: ' + directory+gff_file)
        if not os.path.exists(directory): os.makedirs(directory)

        while os.path.isfile(directory+gff_file):
            print('File already exist!')
            if overwrite:
                break
            else:
                x = input('Type new name or ! to overwrite, a - abort:')
                if x in ['!']:
                    print('Overwriting...')
                    break
                elif x in ['a', 'A']:
                    print('Aborted!')
                    return
                else:
                    gff_file = x
                    
        titles = "grain_id mean_IA chisq x y z U11 U12 U13 U21 U22 U23 U31 U32 U33 UBI11 UBI12 UBI13 UBI21 UBI22 UBI23 UBI31 UBI32 UBI33".split()
        f = open(directory+gff_file ,"w") 
        f.write('# '+' '.join(titles) + '\n')
        for g in list_of_grains:
            values = [g.grain_id, g.mean_IA, g.pos_chisq] + g.position + g.u.flatten().tolist() + g.ubi.flatten().tolist()
            s = [f'{v:d}' if type(v)==type(1) else f'{v:f}' for v in values]
            f.write('  ' + '  '.join(s) + '\n' )
        f.close()
        print('File closed!')
        return
    
    
    def identify_measured(self, gvectors):
        measured_gvectors = []
        for gv in self.gvectors_report:
            gv_list = list(filter(lambda g: g['spot3d_id'] == gv['peak_id'], gvectors))
            if len(gv_list) < 0:
                print(gv)
                raise ValueError('this g-vector not found in the provided list of gvectors!')
            elif len(gv_list) > 1:
                print(gv)
                raise ValueError('more than 1 g-vector was found in the provided list of gvectors!')
            else:
                measured_gvectors.append(gv_list[0])
        self.set_attr('measured_gvectors', measured_gvectors)
        return 
        
    def simulate_gvectors(self, geometry, omega_range, tth_range, beamflux, bckg, psf, peakshape):
        P = py3DXRD.PolySim(directory = self.directory)
        P.set_attr('geometry', geometry)
        P.set_attr('beamflux', beamflux)
        P.set_attr('direc', './sim/')
        P.set_attr('stem', 'sim')
        P.set_attr('grains', [self])
         
        P.set_attr('omega_start', omega_range.start)
        P.set_attr('omega_step', omega_range.step)
        P.set_attr('omega_end', omega_range.end)
        P.set_attr('theta_min', tth_range.start)
        P.set_attr('theta_max', tth_range.end)
        P.set_attr('no_grains', 1)
        P.set_attr('gen_U', 0)
        P.set_attr('gen_pos', [0, 0])
        P.set_attr('gen_eps', [0, 0, 0 ,0, 0])
        P.set_attr('gen_size', [0.0, 0.0, 0.0 ,0.0])
        P.set_attr('sample_cyl', [0.17, 0.2])
        if not P.grains[0].size:
            P.grains[0].set_attr('size',0.05)
        P.set_attr('make_image', 0)
        P.set_attr('output', ['.tif', '.par', '.gve'])
        P.set_attr('bg', bckg)
        P.set_attr('psf', psf)
        P.set_attr('peakshape', peakshape) #[1, 4, 0.5])
        P.save_inp(inp_file = 'sim.inp')
        return
        
        
    def plot_measured_vs_expected(self):
        tth_measured = [gv['tth']   for gv in self.measured_gvectors]
        eta_measured = [gv['eta']   for gv in self.measured_gvectors]
        omg_measured = [gv['omega'] for gv in self.measured_gvectors]
        tth_expected = [gv['tth']   for gv in self.expected_gvectors]
        eta_expected = [gv['eta']   for gv in self.expected_gvectors]
        omg_expected = [gv['omega'] for gv in self.expected_gvectors]
        
        fig = plt.figure(figsize=(16, 10))
        sub1 = fig.add_subplot(131)
        sub1.scatter(omg_expected, eta_expected, s=10, c='k', marker="o", label='Expected')
        sub1.scatter(omg_measured, eta_measured, s=5 , c='r', marker="v", label='Measured')
        # sub1.legend(bbox_to_anchor=(0.4, 1.1))
        plt.title(f'Expected {len(self.expected_gvectors)} (black), measured {len(self.measured_gvectors)} (red)')
        plt.xlabel('omega (deg)')
        plt.ylabel('eta (deg)')

        sub2 = fig.add_subplot(132)
        sub2.scatter(tth_expected, omg_expected, s=10, c='k', marker="o", label='Expected')
        sub2.scatter(tth_measured, omg_measured, s=5 , c='r', marker="v", label='Measured')
        # sub1.legend(bbox_to_anchor=(0.4, 1.1))
        plt.title(f'Expected {len(self.expected_gvectors)} (black), measured {len(self.measured_gvectors)} (red)')
        plt.xlabel('tth (deg)')
        plt.ylabel('omega (deg)')

        sub3 = fig.add_subplot(133)
        sub3.scatter(tth_expected, eta_expected, s=10, c='k', marker="o", label='Expected')
        sub3.scatter(tth_measured, eta_measured, s=5 , c='r', marker="v", label='Measured')
        # sub1.legend(bbox_to_anchor=(0.4, 1.1))
        plt.title(f'Expected {len(self.expected_gvectors)} (black), measured {len(self.measured_gvectors)} (red)')
        plt.xlabel('tth (deg)')
        plt.ylabel('eta (deg)')
        
        plt.show()
        fig.savefig(self.directory+self.log_file.replace('.log', '')+'_g'+str(self.grain_id)+'_scatter.png')
        self.add_to_log('Saved measured vs expected: '+self.log_file.replace('.log', '')+'_g'+str(self.grain_id)+'_scatter.png', True)
        return

#### TESTS    ###############
# from GrainSpotter import GrainSpotter
# list_of_grains = load_log_file(directory='/asap3/petra3/gpfs/p21.2/2021/data/11008399/processed/CS_1/load1/z_00/', 
#              log_file = 'V4_peaks_merged_grains.log')
# list_of_grains[0].print(True)