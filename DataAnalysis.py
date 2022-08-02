# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:00:00 2022

@author: shabalin

Class to perform full 3D XRD analysis
"""

import sys, os, subprocess, pickle
import numpy as np
from datetime import datetime
from py3DXRD import Geometry

single_separator = "--------------------------------------------------------------\n"
double_separator = "==============================================================\n"

class DataAnalysis:
    
    def __init__(self, directory=None, name=None):
        self.log = []
        self.directory = None
        self.name = None
        self.material = None
        self.pressure = None
        self.temperature = None
        self.position = [0,0,0]
        self.rotation = [0,0,0]
        self.sweepscans = []
        self.fltfiles = []
        self.geometries = []
        self.yml_det_order = []
        self.gvefiles = []
        self.grainspotters = []
        self.grains = []
        self.add_to_log('Created DataAnalysis object.', True)
        if directory: self.set_attr('directory', directory)
        if name     : self.set_attr('name'     , name)
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
        print(double_separator+'DataAnalysis object:')
        print('directory:'    , self.directory)
        print('name:'         , self.name)
        print('material:'     , self.name)
        print('pressure:'     , self.pressure)
        print('temperature:'  , self.temperature)
        print('position:'     , self.position)
        print('rotation:'     , self.rotation)
        print('sweepscans:'   , len(self.sweepscans))
        print('fltfiles:'     , len(self.fltfiles))
        print('geometries:'   , len(self.geometries))
        print('yml_det_order:', yml_det_order)
        print('gvefiles:'     , len(self.gvefiles))
        print('grainspotters:', len(self.grainspotters))
        print('grains:'       , len(self.grains))
        if also_log:
            print(single_separator + 'Log:')
            for record in self.log: print(record)
        return
    
    def save_geometries_as_yml(self, yml_det_order = None):
        if yml_det_order: self.set_attr('yml_det_order', yml_det_order)
        gs = [self.geometries[ind] for ind in self.yml_det_order]
        yml_file = Geometry.save_geometries_as_yml(gs, self.directory, self.name+'.yml', overwrite = True)
        self.add_to_log('Exported all the geometries as .yml file: ' + yml_file, True)
        return
        
    def absorb(self, x):
        print('Not implemented yet!')
        return
#### TESTS    ###############