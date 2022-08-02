# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 16:00:49 2022

@author: shabalin

Class to run indexer, convert between different geometry formalisms.
"""

import os, yaml
import numpy as np
from datetime import datetime
# ImageD11 might need running "module load maxwell ImageD11" before starting python
from ImageD11 import indexing
from ImageD11 import transformer
from scipy.spatial.transform import Rotation

single_separator = "--------------------------------------------------------------\n"
double_separator = "==============================================================\n"


class Geometry:
    
    def __init__(self, directory = None):
        self.log = []
        self.directory = None
        self.par_file  = None
        self.yml_file  = None
        self.poni_file = None
        self.flt_file  = None
        self.gve_file  = None
        self.unitcell  = [None, None, None, None, None, None]
        self.symmetry  = None
        self.spacegroup= None
        self.chi       = None
        self.distance  = None
        self.buffer    = None
        self.saturation_level = None
        self.fit_tolerance = None
        self.min_bin_prob  = None
        self.no_bins   = None
        self.O         = [[None, None], [None, None]]
        self.omegasign = None
        self.t     = [0,0,0]
        self.tilt  = [None, None, None]
        self.wedge = None
        self.weight_hist_intensities = None
        self.wavelength = None
        self.y_center   = None
        self.y_size     = None
        self.dety_size  = None
        self.detz_size  = None
        self.z_center   = None
        self.z_size     = None
        self.add_to_log('Created Geometry object.', True)
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
        self.add_to_log(attr+': '+str(old)+' -> '+str(new))
        return
        
        
    def add_to_attr(self, attr, value):
        old_list = getattr(self, attr)
        setattr(self, attr, old_list+[value])
        new_list = getattr(self, attr)
        self.add_to_log(attr+': += '+str(new_list[-1]))
        return
    
    
    def print(self, also_log = False):
        print(double_separator+'Geometry object:')
        print('directory:', self.directory)
        print('par_file:' , self.par_file )
        print('yml_file:' , self.yml_file )
        print('poni_file:', self.poni_file)
        print('flt_file:' , self.flt_file )
        print('gve_file:' , self.gve_file )
        print('unitcell:' , self.unitcell )
        print('symmetry:' , self.symmetry )
        print('spacegroup:', self.spacegroup )
        print('chi:'      , self.chi)
        print('distance:' , self.distance)
        print('buffer:'   , self.buffer)
        print('saturation_level:', self.saturation_level)
        print('fit_tolerance:', self.fit_tolerance)
        print('min_bin_prob:' , self.min_bin_prob )
        print('no_bins:'  , self.no_bins)
        print('O: '       , self.O)
        print('omegasign:', self.omegasign)
        print('t:'        , self.t)
        print('tilt:'     , self.tilt )
        print('wedge:'    , self.wedge)
        print('weight_hist_intensities:', self.weight_hist_intensities)
        print('wavelength:', self.wavelength)
        print('y_center:'  , self.y_center)
        print('y_size:'    , self.y_size  )
        print('dety_size:' , self.dety_size  )
        print('z_center:'  , self.z_center)
        print('z_size:'    , self.z_size  )
        print('detz_size:' , self.detz_size  )
        if also_log:
            print(single_separator + 'Log:')
            for record in self.log: print(record)
        return
    

    def load_par(self, directory = None, par_file = None):
        if directory: self.set_attr('directory', directory)
        if par_file: self.set_attr('par_file', par_file)
        self.add_to_log(f'Reading file: {self.directory+self.par_file}', True)
        if not os.path.isfile(self.directory+self.par_file): raise FileNotFoundError  
        unitcell = [None, None, None, None, None, None]
        O = [[None, None], [None, None]]
        t = [None, None, None]
        tilt = [None, None, None]
        with open(self.directory+self.par_file, "r") as f:
            for line in f:
                if line[0] == '#' or len(line) < 2: continue
                words = line.split() # words in line
                if   ('cell__a'    == words[0]): unitcell[0] = float(words[1])
                elif ('cell__b'    == words[0]): unitcell[1] = float(words[1])
                elif ('cell__c'    == words[0]): unitcell[2] = float(words[1])
                elif ('cell_alpha' == words[0]): unitcell[3] = float(words[1])
                elif ('cell_beta'  == words[0]): unitcell[4] = float(words[1])
                elif ('cell_gamma' == words[0]): unitcell[5] = float(words[1])
                elif ('cell_lattice' in words[0]): self.set_attr('symmetry', words[1])
                elif ('chi' == words[0]): self.set_attr('chi', float(words[1]))
                elif ('distance' == words[0]): self.set_attr('distance', float(words[1]))
                elif ('fit_tolerance' == words[0]):
                    self.set_attr('fit_tolerance', float(words[1]))
                elif ('min_bin_prob' == words[0]):
                    self.set_attr('min_bin_prob', float(words[1]))
                elif ('no_bins' == words[0]): self.set_attr('no_bins', int(words[1]))
                elif ('o11' == words[0]): O[0][0] = int(words[1])
                elif ('o12' == words[0]): O[0][1] = int(words[1])
                elif ('o21' == words[0]): O[1][0] = int(words[1])
                elif ('o22' == words[0]): O[1][1] = int(words[1])
                elif ('omegasign' == words[0]):
                    self.set_attr('omegasign', int(float(words[1])))
                elif ('tilt_x' == words[0]): tilt[0] = float(words[1])
                elif ('tilt_y' == words[0]): tilt[1] = float(words[1])
                elif ('tilt_z' == words[0]): tilt[2] = float(words[1])
                elif ('t_x'    == words[0]): t[0]    = float(words[1])
                elif ('t_y'    == words[0]): t[1]    = float(words[1])
                elif ('t_z'    == words[0]): t[2]    = float(words[1])
                elif ('wedge' == words[0]): self.set_attr('wedge', float(words[1]))
                elif ('weight_hist_intensities' == words[0]):
                    if words[1] in ['False', 'None']: self.set_attr('weight_hist_intensities', False)
                    else: self.set_attr('weight_hist_intensities', int(words[1]))
                elif ('wavelength' == words[0]):
                    self.set_attr('wavelength', float(words[1]))
                elif ('y_center'  == words[0]): self.set_attr('y_center' , float(words[1]))
                elif ('y_size'    == words[0]): self.set_attr('y_size'   , float(words[1]))
                elif ('dety_size' == words[0]): self.set_attr('dety_size', int(words[1]))
                elif ('z_center'  == words[0]): self.set_attr('z_center' , float(words[1]))
                elif ('z_size'    == words[0]): self.set_attr('z_size'   , float(words[1]))
                elif ('detz_size' == words[0]): self.set_attr('detz_size', int(words[1]))
            self.set_attr('unitcell', unitcell)
            self.set_attr('O', O)
            self.set_attr('t', t)
            self.set_attr('tilt', tilt)
        f.close()
        self.add_to_log('File closed!', False)
        
        while None in np.asarray(self.O):
            print('Image orientation (O-matrix) is missing!')
            x = input('Set it as 4 spaced numbers (O11 O12 O21 O22):').split()
            self.O = [[int(v) for v in x[0:2]], [int(v) for v in x[2:]]]
            
        return
 
        
    def save_par(self, directory = None, par_file = None, overwrite = False):
        if directory: self.set_attr('directory', directory)
        if par_file: self.set_attr('par_file', par_file)
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            self.add_to_log('Created directory: '+self.directory, True)
        self.add_to_log(f'Writing file: {self.directory+self.par_file}', True)
        
        # Check if file exists, if yes then ask for the permission to overwrite
        while os.path.isfile(self.directory+self.par_file):
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
                    self.set_attr('par_file', x)
        
        f = open(self.directory+self.par_file ,"w")
        f.write( 'cell__a {}\n'.format(self.unitcell[0]) )
        f.write( 'cell__b {}\n'.format(self.unitcell[1]) )
        f.write( 'cell__c {}\n'.format(self.unitcell[2]) )
        f.write( 'cell_alpha {}\n'.format(self.unitcell[3]) )
        f.write( 'cell_beta {}\n'.format(self.unitcell[4]) )
        f.write( 'cell_gamma {}\n'.format(self.unitcell[5]) )
        f.write( 'cell_lattice_[P,A,B,C,I,F,R] {}\n'.format(self.symmetry) )
        f.write( 'chi {}\n'.format(self.chi) )
        f.write( 'distance {:0.3f}\n'.format(self.distance) )
        f.write( 'fit_tolerance {}\n'.format(self.fit_tolerance) )
        if self.min_bin_prob: f.write( 'min_bin_prob {}\n'.format(self.min_bin_prob) )
        if self.no_bins: f.write( 'no_bins {}\n'.format(self.no_bins) )
        f.write( 'o11 {}\n'.format(self.O[0][0]) )
        f.write( 'o12 {}\n'.format(self.O[0][1]) )
        f.write( 'o21 {}\n'.format(self.O[1][0]) )
        f.write( 'o22 {}\n'.format(self.O[1][1]) )
        f.write( 'omegasign {}\n'.format(self.omegasign) )
        f.write( 't_x {}\n'.format(self.t[0]) )
        f.write( 't_y {}\n'.format(self.t[1]) )
        f.write( 't_z {}\n'.format(self.t[2]) )
        f.write( 'tilt_x {:0.6f}\n'.format(self.tilt[0]) )
        f.write( 'tilt_y {:0.6f}\n'.format(self.tilt[1]) )
        f.write( 'tilt_z {:0.6f}\n'.format(self.tilt[2]) )
        f.write( 'wavelength {:0.6f}\n'.format(self.wavelength) )
        f.write( 'wedge {}\n'.format(self.wedge) )
        if self.weight_hist_intensities: f.write( 'weight_hist_intensities {}\n'.format(self.weight_hist_intensities) )
        f.write( 'y_center {:0.3f}\n'.format(self.y_center) )
        f.write( 'y_size {}\n'.format(self.y_size) )
        if self.dety_size: f.write( 'dety_size {}\n'.format(self.dety_size) )
        f.write( 'z_center {:0.3f}\n'.format(self.z_center) )
        f.write( 'z_size {}\n'.format(self.z_size) )
        if self.detz_size: f.write( 'detz_size {}\n'.format(self.detz_size) )
        f.close()
        self.add_to_log('File closed!', False)
        return
        
        
    def run_indexer(self, directory = None, flt_file = None, gve_file = None):
        if directory: self.set_attr('directory', directory)
        if flt_file: self.set_attr('flt_file', flt_file)
        if gve_file: self.set_attr('gve_file', gve_file)
        self.save_par(overwrite = True)
        self.add_to_log(f'Running indexer on: {self.directory}{self.flt_file}', True)
        self.add_to_log(f'Using par_file: {self.directory}{self.par_file}', True)
        obj = transformer.transformer()
        obj.loadfiltered( self.directory+self.flt_file )
        obj.loadfileparameters( self.directory+self.par_file )
        obj.compute_tth_eta( )
        obj.addcellpeaks( )
        obj.computegv( )
        obj.savegv( self.directory+self.gve_file )
        self.add_to_log(f'Resulted vectors saved in: {self.directory}{self.gve_file}', True)
        return

    
    def in_hexrd_definitions(self, det_num=1):
        while None in np.asarray(self.O):
            print('Image orientation (O-matrix) is missing!')
            x = input('Set it as 4 spaced numbers (O11 O12 O21 O22):').split()
            self.O = [[int(v) for v in x[0:2]], [int(v) for v in x[2:]]]
        ### HEXRD[x,y,z] = FABLE[-y,z,-x] ###
        r0 = Rotation.from_euler('z', -self.tilt[0])  # 1st is about FABLE:+x (radians) = about HEXRD:-z
        r1 = Rotation.from_euler('x', -self.tilt[1])  # 2nd is about FABLE:+y (radians) = about HEXRD:-x
        r2 = Rotation.from_euler('y',  self.tilt[2])  # 3rd is about FABLE:+z (radians) = about HEXRD:+y
        hexrd_tilt = (r0*r1*r2).as_euler('xyz') # radians
        
        det_imgsize_sf = np.asarray([self.detz_size  , self.dety_size  ])   # (det_slow, det_fast)
        det_pixsize_sf = np.asarray([self.z_size     , self.y_size     ])   # in um here
        det_beampos_sf = np.asarray([self.z_center   , self.y_center   ])   # in pixels
        det_centpos_sf = np.asarray([self.detz_size-1, self.dety_size-1])/2 # in pixels

        det_beampos_sf = (det_beampos_sf-det_centpos_sf)*det_pixsize_sf   # in um from img. center
        
        fable_imgsize_zy = np.matmul(np.asarray(self.O), det_imgsize_sf) # rot. about img. center
        fable_pixsize_zy = np.matmul(np.asarray(self.O), det_pixsize_sf) # rot. about img. center
        fable_beampos_zy = np.matmul(np.asarray(self.O), det_beampos_sf) # rot. about img. center

        hexrd_imgsize = [round(abs(fable_imgsize_zy[0])), round(abs(fable_imgsize_zy[1]))]
        hexrd_pixsize = [      abs(fable_pixsize_zy[0]) ,       abs(fable_pixsize_zy[1]) ]
        
        hexrd_centpos_befor_rot = np.asarray( [fable_beampos_zy[1], -fable_beampos_zy[0], 0] )
        hexrd_centpos_after_rot = (r0*r1*r2).apply(hexrd_centpos_befor_rot)
        hexrd_translt = hexrd_centpos_after_rot + np.asarray([0, 0, -self.distance])
        
        if   self.O == [[-1, 0], [ 0,-1]]: hexrd_orientation = 'none'
        elif self.O == [[ 0,-1], [-1, 0]]: hexrd_orientation = 't'
        elif self.O == [[ 1, 0], [ 0,-1]]: hexrd_orientation = 'v'
        elif self.O == [[-1, 0], [ 0, 1]]: hexrd_orientation = 'h'
        elif self.O == [[ 1, 0], [ 0, 1]]: hexrd_orientation = 'r180'
        elif self.O == [[ 0, 1], [-1, 0]]: hexrd_orientation = 'r90'
        elif self.O == [[ 0,-1], [ 1, 0]]: hexrd_orientation = 'r270'
        else                             : hexrd_orientation = 'unknown'
        
        energy = 12.39842/self.wavelength
        azimuth = 90
        if self.wedge: polar_angle = 90.0-np.degrees(self.wedge) # must be checked before using!
        else         : polar_angle = 90.0
        if self.chi: chi = np.degrees(self.chi) # must be checked before using!
        else       : chi = 0
        
        saturation_level = self.saturation_level
        buffer = self.buffer

        pars = {'id': 'instrument',
                'beam':{
                    'energy': energy,
                    'vector':{
                        'azimuth'    : 90.0,   # must be checked before using!
                        'polar_angle': polar_angle}}, # must be checked before using!
                'oscillation_stage':{
                    'chi': chi, 
                    'translation': [0, 0, 0]}, # in FABLE t is the translation of a grain(not of the stage)
                'detectors':{
                    f'detector_{det_num}':{
                        'transform':{
                            'tilt': [float(v) for v in hexrd_tilt],
                            'translation': [float(v/1000) for v in hexrd_translt],
                            'orientation': hexrd_orientation},
                        'pixels':{
                            'rows'   : int(hexrd_imgsize[0]),
                            'columns': int(hexrd_imgsize[1]),
                            'size'   : [float(v/1000) for v in hexrd_pixsize]},
                        'saturation_level': saturation_level,
                        'buffer': buffer}}}       
        return pars

    
    def load_yml(self, directory = None, yml_file = None, det_num = 1):
        if not self.O: raise ValueError('O-matix is missing!')
        if directory: self.set_attr('directory', directory)
        if yml_file: self.set_attr('yml_file', yml_file)
        self.add_to_log(f'Reading file: {self.directory+self.yml_file}', True)
        if not os.path.isfile(self.directory+self.yml_file): raise FileNotFoundError
        with open(self.directory+self.yml_file, "r") as f:
            pars = yaml.safe_load(f)
        det_name = 'detector_{:d}'.format(det_num)
        stage_t = pars['oscillation_stage']['translation']
        translation = pars['detectors'][det_name]['transform']['translation']*1000 # converted to um
        tilt        = pars['detectors'][det_name]['transform']['tilt'] # in radians
        pixels      = pars['detectors'][det_name]['pixels']
        
        if 'orientation' in pars['detectors'][det_name]['transform'].keys():
            hexrd_orientation = pars['detectors'][det_name]['transform']['orientation']
            if   hexrd_orientation == 'none': self.O == [[-1, 0], [ 0,-1]]
            elif hexrd_orientation == 't'   : self.O == [[ 0,-1], [-1, 0]]
            elif hexrd_orientation == 'v'   : self.O == [[ 1, 0], [ 0,-1]]
            elif hexrd_orientation == 'h'   : self.O == [[-1, 0], [ 0, 1]]
            elif hexrd_orientation == 'r180': self.O == [[ 1, 0], [ 0, 1]]
            elif hexrd_orientation == 'r90' : self.O == [[ 0, 1], [-1, 0]]
            elif hexrd_orientation == 'r270': self.O == [[ 0,-1], [ 1, 0]]
            print(pars['detectors'][det_name]['transform'])
        while None in np.asarray(self.O):
            print('Image orientation (O-matrix) is missing!')
            x = input('Set it as 4 spaced numbers (O11 O12 O21 O22):').split()
            self.O = [[int(v) for v in x[0:2]], [int(v) for v in x[2:]]]
                                     
        det_pixsize_sf = np.asarray([pixels['size'][0], pixels['size'][1]])*1000 # converted to um
        det_imgsize_sf = np.asarray([pixels['rows']   , pixels['columns']])

        fable_pixsize_zy = np.matmul(np.linalg.inv(self.O), det_pixsize_sf) # rot. about img. center
        fable_imgsize_zy = np.matmul(np.linalg.inv(self.O), det_imgsize_sf) # rot. about img. center
        fable_pixsize_sf =         abs(fable_pixsize_zy)
        fable_imgsize_sf = np.rint(abs(fable_imgsize_zy))

        ### HEXRD[x,y,z] = FABLE[-y,z,-x] ###
        r0 = Rotation.from_euler('y', -tilt[0])  # 1st is about HEXRD:+x (radians) = about FABLE:-y
        r1 = Rotation.from_euler('z',  tilt[1])  # 2nd is about HEXRD:+y (radians) = about FABLE:+z
        r2 = Rotation.from_euler('x', -tilt[2])  # 3rd is about HEXRD:+z (radians) = about FABLE:-x
        fable_tilt = (r2*r1*r0).as_euler('zyx')[::-1] # radians

        r0 = Rotation.from_euler('x', tilt[0])  # 1st is about HEXRD:+x (radians)
        r1 = Rotation.from_euler('y', tilt[1])  # 2nd is about HEXRD:+y (radians)
        r2 = Rotation.from_euler('z', tilt[2])  # 3rd is about HEXRD:+z (radians) 

        det_norm = (r2*r1*r0).apply( np.asarray([0,0,1]) ) # normal to the detector plane after tilts
        z_shift = (det_norm[0]*translation[0] + det_norm[1]*translation[1])/det_norm[2]
        fable_distance = - (translation[2] + z_shift) # FABLE +x is HEXRD -z

        hexrd_centpos_befor_rot = np.asarray([translation[0], translation[1], -z_shift])
        hexrd_centpos_after_rot = np.matmul( np.linalg.inv( (r2*r1*r0).as_matrix()), hexrd_centpos_befor_rot)
        fable_beampos_zy = np.asarray([-hexrd_centpos_after_rot[1], hexrd_centpos_after_rot[0]])
        fable_beampos_sf = np.matmul(np.linalg.inv(self.O), fable_beampos_zy)/fable_pixsize_sf # rot. about img. center
        fable_centpos_sf = fable_beampos_sf + (fable_imgsize_sf-1)/2
        fable_wavelength = 12.39842/pars['beam']['energy'] # keV
        fable_wedge = 90-pars['beam']['vector']['polar_angle'] # degrees, must be checked before using!
        fable_chi = pars['oscillation_stage']['chi'] # degrees, must be checked before using!

        self.set_attr('wavelength', fable_wavelength)
        self.set_attr('wedge'     , np.radians(fable_wedge))
        self.set_attr('chi'       , np.radians(fable_chi))
        self.set_attr('buffer'    , pars['detectors'][det_name]['buffer'])
        self.set_attr('saturation_level', pars['detectors'][det_name]['saturation_level'])                                             
        self.set_attr('tilt'     , fable_tilt)
        self.set_attr('distance' , fable_distance) # converted to um
        self.set_attr('y_size'   , fable_pixsize_sf[1]) # converted to um
        self.set_attr('z_size'   , fable_pixsize_sf[0]) # converted to um
        self.set_attr('dety_size', int(fable_imgsize_sf[1]))
        self.set_attr('detz_size', int(fable_imgsize_sf[0]))
        self.set_attr('y_center' , fable_centpos_sf[1])
        self.set_attr('z_center' , fable_centpos_sf[0])
        
        self.add_to_log('File closed!', False)
        return
    
    
    def save_yml(self, directory = None, yml_file = None, overwrite = False):
        if directory: self.set_attr('directory', directory)
        if yml_file: self.set_attr('yml_file', yml_file)
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            self.add_to_log('Created directory: '+self.directory, True)
        self.add_to_log(f'Writing file: {self.directory+self.yml_file}', True)
        
        # Check if file exists, if yes then ask for the permission to overwrite
        while os.path.isfile(self.directory+self.yml_file):
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
                    self.set_attr('yml_file', x)

        pars = self.in_hexrd_definitions(det_num = 1)
        with open(self.directory+self.yml_file, "w") as f:
            yaml.dump(pars, f)
        self.add_to_log('File closed!', False)
        return                  
  

    def load_poni(self, directory = None, poni_file = None):
        if directory: self.set_attr('directory', directory)
        if poni_file: self.set_attr('poni_file', poni_file)
        self.add_to_log(f'Reading file: {self.directory+self.poni_file}', True)
        if not os.path.isfile(self.directory+self.poni_file): raise FileNotFoundError

        with open(self.directory+self.poni_file, "r") as f:
            for line in f:
                if len(line) < 2: continue
                if line[0] == '#' and 'orientation' in line:
                    m = line.split('[[')[1].split(']]')[0]
                    O = [int(v) for v in m.replace('], [',' ').replace(', ',' ').split()]
                    self.O = [[O[0], O[1]], [O[2], O[3]]]
                words = line.split() # words in line
                if('Detector_config' in words[0]):
                    pixel1 = 1e6*float(words[2].replace(',', ''))
                    pixel2 = 1e6*float(words[4].replace(',', ''))
                    max_shape1 = int(words[6].replace('[', '').replace(',', ''))
                    max_shape2 = int(words[7].replace(']', '').replace('}', ''))
                if('Distance' in words[0]): distance = 1e6*float(words[1])
                if('Poni1' in words[0]): poni1 = 1e6*float(words[1])
                if('Poni2' in words[0]): poni2 = 1e6*float(words[1])
                if('Rot1' in words[0]): rot1 = float(words[1])
                if('Rot2' in words[0]): rot2 = float(words[1])
                if('Rot3' in words[0]): rot3 = float(words[1])
                if('Wavelength' in words[0]): wavelength = 1e10*float(words[1])
        f.close()
        self.add_to_log('File closed!', False)
        
        while None in np.asarray(self.O):
            print('Image orientation (O-matrix) is missing!')
            x = input('Set it as 4 spaced numbers (O11 O12 O21 O22):').split()
            self.O = [[int(v) for v in x[0:2]], [int(v) for v in x[2:]]]
        
        self.wavelength = wavelength
        self.dety_size = max_shape2
        self.detz_size = max_shape1
        self.y_size = pixel2
        self.z_size = pixel1
        fable_tiltrot_zy = np.matmul(np.asarray(self.O), np.asarray([rot1, rot2])) # rot. about img. center
        self.tilt = [rot3, fable_tiltrot_zy[1], fable_tiltrot_zy[0]]
        self.distance = distance/np.cos(rot1)/np.cos(rot2)
        self.y_center = -0.5 + (poni2-distance*np.tan(rot1))/pixel2
        self.z_center = -0.5 + (poni1+distance*np.tan(rot2)/np.cos(rot1))/pixel1
        return

    
    def save_poni(self, directory = None, poni_file = None, overwrite = False):
        if directory: self.set_attr('directory', directory)
        if poni_file: self.set_attr('poni_file', poni_file)
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            self.add_to_log('Created directory: '+self.directory, True)
        self.add_to_log(f'Writing file: {self.directory+self.poni_file}', True)

        # Check if file exists, if yes then ask for the permission to overwrite
        while os.path.isfile(self.directory+self.poni_file):
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
                    self.set_attr('poni_file', x)

        wavelength = 1e-10*self.wavelength
        rot1,rot2 = np.matmul(np.linalg.inv(self.O), np.asarray([self.tilt[2], self.tilt[1]])) # rot. about img. center
        rot3 = self.tilt[0]
        pixel1 = 1e-6*self.z_size
        pixel2 = 1e-6*self.y_size
        max_shape = [self.detz_size, self.dety_size]
        distance = 1e-6*self.distance*np.cos(rot1)*np.cos(rot2)
        poni1 = (self.z_center+0.5)*pixel1 - distance*np.tan(rot2)/np.cos(rot1)
        poni2 = (self.y_center+0.5)*pixel2 + distance*np.tan(rot1)
        
        pars = self.in_hexrd_definitions()
        hexrd_orientation = pars['detectors']['detector_1']['transform']['orientation']

        f = open(self.directory+self.poni_file ,"w")
        f.write('# Nota: C-Order, 1 refers to the Y axis, 2 to the X axis \n')
        f.write(f'# Converted from FABLE using image orientation: {str(self.O)} = {hexrd_orientation} in HEXRD\n')
        f.write('poni_version: 2\nDetector: Detector\n')

        Detector_config = {"pixel1": pixel1, "pixel2": pixel2, "max_shape": max_shape}
        f.write('Detector_config: ' + str(Detector_config).replace('\'', '\"') + '\n')
        f.write(f'Distance: {distance}\n')
        f.write(f'Poni1: {poni1}\n')
        f.write(f'Poni2: {poni2}\n')
        f.write(f'Rot1: {rot1}\n')
        f.write(f'Rot2: {rot2}\n')
        f.write(f'Rot3: {rot3}\n')
        f.write(f'Wavelength: {wavelength}\n')
        f.close()
        self.add_to_log('File closed!', False)
        return self    
    
    
    def average_geometries(list_of_geometries):
        from statistics import mean
        GM = list_of_geometries[0]
        n = len(list_of_geometries)
        GM.wavelength = mean([g.wavelength for g in list_of_geometries])
        GM.unitcell[0] = mean([g.unitcell[0] for g in list_of_geometries])
        GM.unitcell[1] = mean([g.unitcell[1] for g in list_of_geometries])
        GM.unitcell[2] = mean([g.unitcell[2] for g in list_of_geometries])
        GM.distance = mean([g.distance for g in list_of_geometries])
        GM.tilt[0] = mean([g.tilt[0] for g in list_of_geometries])
        GM.tilt[1] = mean([g.tilt[1] for g in list_of_geometries])
        GM.tilt[2] = mean([g.tilt[2] for g in list_of_geometries])
        GM.t[0] = mean([g.t[0] for g in list_of_geometries])
        GM.t[1] = mean([g.t[1] for g in list_of_geometries])
        GM.t[2] = mean([g.t[2] for g in list_of_geometries])
        GM.y_center = mean([g.y_center for g in list_of_geometries])
        GM.z_center = mean([g.z_center for g in list_of_geometries])
        for g in list_of_geometries:
    #         if g.directory != GM.directory: raise ValueError('Directories are not consistent!')
            if g.unitcell[3:] != GM.unitcell[3:]: raise ValueError('Unit cell anlges are not consistent!')
            if g.symmetry != GM.symmetry: raise ValueError('Symmetries are not consistent!')
            if g.spacegroup != GM.spacegroup: raise ValueError('Spacegroups are not consistent!')
    #         if g.chi != GM.chi: raise ValueError('chis are not consistent!')
            if g.fit_tolerance != GM.fit_tolerance: raise ValueError('fit_tolerances are not consistent!')
            if g.min_bin_prob != GM.min_bin_prob: raise ValueError('min_bin_probs are not consistent!')
            if g.no_bins != GM.no_bins: raise ValueError('no_bins are not consistent!')
            if g.O != GM.O: raise ValueError('O-matrices are not consistent!')
            if g.omegasign != GM.omegasign: raise ValueError('omegasigns are not consistent!')
    #         if g.wedge != GM.wedge: raise ValueError('wedges are not consistent!')
            if g.weight_hist_intensities != GM.weight_hist_intensities: raise ValueError('weight_hist_intensities are not consistent!')
            if g.y_size != GM.y_size: raise ValueError('y_sizes are not consistent!')
            if g.z_size != GM.z_size: raise ValueError('z_sizes are not consistent!')
            if g.dety_size != GM.dety_size: raise ValueError('dety_sizes are not consistent!')
            if g.detz_size != GM.detz_size: raise ValueError('detz_sizes are not consistent!')
        return GM

    def save_geometries_as_yml(list_of_geometries, directory, yml_file, overwrite = False):
        GM = list_of_geometries[0]
        pars = GM.in_hexrd_definitions(det_num = 1)
        for det_num, gm in enumerate(list_of_geometries):
            if gm.wavelength != GM.wavelength: print('Warning! Wavelengths are not consistent!')
            if gm.omegasign  != GM.omegasign:  print('Warning! Omegasigns are not consistent!')        
            if gm.wedge      != GM.wedge:      print('Warning! Wedges are not consistent!')
            if gm.chi        != GM.chi:        print('Warning! Chis are not consistent!')
            if gm.O          != GM.O  :        print('Warning! O-matrices are not consistent!')
            if gm.t          != GM.t  :        print('Warning! Stage translations are not consistent!')

            p = gm.in_hexrd_definitions(det_num = 1)
            if p['beam'] != pars['beam']: raise ValueError('Beam parameters are not consistent!')
            if p['id']   != pars['id']  : raise ValueError('Oscillation stages are not consistent!')
            pars['detectors'][f'detector_{det_num+1}'] = p['detectors'][f'detector_{1}']

        # Check if file exists, if yes then ask for the permission to overwrite    
        while os.path.isfile(directory+yml_file):       
            print('Warning! File already exist!')
            if overwrite:
                print('Overwriting...')
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
                    yml_file =  x

        with open(directory+yml_file, "w") as f:
            yaml.dump(pars, f)

        return directory+yml_file

#### TESTS    ###############
# G = Geometry(directory='/asap3/petra3/gpfs/p21.2/2021/data/11008399/processed/CS_1/load1/z_00/')
# G.load_par(par_file = 'V4_iconfig_Ti_fitted.par')
# G.save_par(par_file = 'V4_iconfig_Ti_fitted_test.par')
# G.run_indexer(flt_file = 'V4_peaks_t8000.flt', gve_file = 'V4_iconfig_Ti_fitted_test.gve')
# G.load_yml(yml_file = 'V4_iconfig_Ti_fitted.yml')
# G.save_yml(yml_file = 'V4_iconfig_Ti_fitted_test.yml')
# G.save_par(par_file = 'V4_iconfig_Ti_fitted_test.par')
# G.print(True)