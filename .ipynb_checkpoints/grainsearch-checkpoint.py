# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 15:46:49 2021

@author: shabalin

Class to work with g-vectors.
"""


import os
import numpy as np
from statistics import mean
from .mod_360 import mod_360
from .merge_overlaps import merge_overlaps
from .group_to_chains import group_to_chains
from ImageD11 import indexing # this require running "module load maxwell ImageD11" before starting python
from ImageD11 import transformer # this require running "module load maxwell ImageD11" before starting python


class GrainSearch:
    
    
    def __init__(self, spacegroup):
        self.path = '/asap3/petra3/gpfs/p21.2/2021/data/11012254/Nb_5/100_sub/' # path to the directory with the .gve files
        self.gve_files = ['V1_peaks_merged.gve', 'V2_peaks_merged.gve'] # list of .gve files to load. Has to be a list!
        self.out_stem = 'V1V2_peaks_merged' # stem for output files  (like .ini and .log)
        self.spacegroup = 229 # [space group nr]
        self.tth_ranges = [] # [min max], multiple ranges can be specified
        self.cuts = [5, 0.3, 0.8] # [(min_measuments:grain is chosen if there are at least this amount of peaks per grain) (min_completeness: grain is chosen if there are at least this fraction of the expected peaks present) (min_uniqueness: no idea, just put any number)]
        self.eulerstep = 6 # [stepsize] : angle step size in Euler space
        self.uncertainties = [0.1, 0.2, 0.5] # [sigma_tth sigma_eta sigma_omega] in degrees
        self.nsigmas = 2 # [Nsig] : maximal deviation in sigmas
        self.Nhkls_in_indexing = None # [Nfamilies] : use first Nfamilies in indexing
        self.random = 10000 # random sampling of orientation space trying 10000 sample points
        self.positionfit = True # fit the position of the grain
        self.minfracg = False # stop search when minfracg (0..1) of the gvectors have been assigned to grains
        self.genhkl = False # generate list of reflections
        self.gve_files = []
        self.gves = [] # ndarray
        self.col_names_string = '#  gx  gy  gz  xc  yc  ds  eta  omega  spot3d_id  xl  yl  zl\n'
        self.header = []
        self.dshkls = []
        self.omega_step = 0
        self.eta_step = 0.5
        self.tth_ranges = []
        self.ds_ranges = []
        self.eta_ranges = []
        self.omega_ranges = []
        self.ds_tol = 0
        self.eta_tol = 0
        self.omega_tol = 0
        self.ds_gap = 0
        self.eta_gap = 0
        self.omega_gap = 0
        self.ring = []
        
        if spacegroup > 0:
            self.spacegroup = spacegroup
        
        return


    def add_gves_from_file(self, filename):
        if not os.path.exists(self.path+filename):
            raise ValueError(self.path + filename + ' - no such file or directory!')
        if '.gve' not in filename:
            raise ValueError(filename + ' - is not .gve file!')

        f = open(self.path+filename,"r")
        line = f.readline()
        while line and '# ds h k l' not in line:
            if line not in self.header:
                self.header.append(line)
            line = f.readline()
        while line and 'spot3d_id' not in line:
            if line not in self.dshkls:
                self.dshkls.append(line)
            line = f.readline()
            
        for line in f.readlines():
            try:
                g = [float(x) for x in line.split()]
                if len(g) == 12:
                    g[6] = mod_360(g[6],180) # comvert eta in [0, 360] range
                    g[7] = mod_360(g[7],0) # convert omega in [-180, 180] range
                    self.gves.append(g)
                    """
                    if len(self.gves) < 1:
                        self.gves = np.reshape(g, (1, -1))
                    else:
                        self.gves = np.append(self.gves, [g], 0)

                    """
                else:
                    print(line, ' -  empty or corrupted line!')
            except:
                raise ValueError(line, ' - could not read this line!')
        f.close()
        self.gve_files.append(filename)
        
        del line
        return
    
    
    def index(self):
        I = indexing.indexer()

        out_stem = self.out_stem
        self.write_gve_file(gve_file = 'temp')
        self.out_stem = out_stem

        I.readgvfile(self.path + "temp.gve") 
        os.remove(self.path + "temp.gve")
        I.assigntorings()
        in_rings = np.compress(np.greater(I.ra,-1),np.arange(I.gv.shape[0])) # list of indexed peaks
        self.gves = [self.gves[i] for i in in_rings]
        self.ring = [I.ra[i] for i in in_rings]
        
        del I, in_rings
        return

    
    def calculate_ranges(self, ds_gap, eta_gap, omega_gap):
        ds_P, eta_P, omega_P = 5,6,7
        ds_list    = [g[ds_P]    for g in self.gves]
        eta_list   = [g[eta_P]   for g in self.gves]
        omega_list = [g[omega_P] for g in self.gves]
        
        ds_chains    = group_to_chains(ds_list   , ds_gap)
        eta_chains   = group_to_chains(eta_list  , eta_gap)
        omega_chains = group_to_chains(omega_list, omega_gap)
        self.ds_ranges    = [ [ds_list[c[0]]-ds_gap      , ds_list[c[-1]]+ds_gap      ] for c in ds_chains    ]
        self.eta_ranges   = [ [eta_list[c[0]]-eta_gap    , eta_list[c[-1]]+eta_gap    ] for c in eta_chains   ]
        self.omega_ranges = [ [omega_list[c[0]]-omega_gap, omega_list[c[-1]]+omega_gap] for c in omega_chains ]

        del ds_list, eta_list, omega_list, ds_chains, eta_chains, omega_chains
        return

    """
    def calculate_ranges(self, ds_gap, eta_gap, omega_gap):
        if type(self.gves) == list: # check if the input is a list
            if type(self.gves[0]) == list: # check if it's a list of ranges
                pass
            elif len(self.gves) != 12:
                raise ValueError('Input g-vector has not 12 elements!')

        ds_rs = []
        et_rs = []
        om_rs = []
        for g in self.gves:
            ds_rs.append( [g[5]-ds_gap,g[5]+ds_gap] )
            et_rs.append( [g[6]-eta_gap,g[6]+eta_gap] )
            om_rs.append( [g[7]-omega_gap,g[7]+omega_gap] )
        self.ds_ranges = merge_overlaps(ds_rs, margin=ds_gap, target=None)
        self.eta_ranges = merge_overlaps(et_rs, margin=eta_gap, target=180)
        self.omega_ranges = merge_overlaps(om_rs, margin=omega_gap, target=0)
            
        self.ds_gap = ds_gap
        self.eta_gap = eta_gap
        self.omega_gap = omega_gap
        
        del ds_rs, et_rs, om_rs
        return
         
    def filter_gves(self, ds_tol, eta_tol, omega_tol):
        # Groups g-vectors that are close to each other takes, returns a list of such 'averaged' g-vectors.
        if ds_tol and eta_tol and omega_tol:
            g_sum = [self.gves[0]]
            n_sum = [1]
            for iP in range(1,len(self.gves)):
                to_add = True
                for iB in range(0,len(n_sum)):
                    if abs((g_sum[iB][5]/n_sum[iB]) - self.gves[iP][5]) < ds_tol:
                        d_eta = mod_360( (g_sum[iB][6]/n_sum[iB])-self.gves[iP][6], 0 )
                        if abs(d_eta) < eta_tol:
                            d_omega = mod_360( (g_sum[iB][7]/n_sum[iB])-self.gves[iP][7], 0 )
                            if abs(d_omega) < omega_tol:
                                g_sum[iB] = [g_sum[iB][i]+self.gves[iP][i] for i in range(0,len(g_sum[iB]))]
                                n_sum[iB] = n_sum[iB] +1
                                to_add = False
                if to_add:
                    g_sum.append(self.gves[iP])
                    n_sum.append(1)
            for iB in range(0,len(n_sum)):
                g_sum[iB] = [g/n_sum[iB] for g in g_sum[iB]]
            self.gves = g_sum
            self.ds_tol = ds_tol
            self.eta_tol = eta_tol
            self.omega_tol = omega_tol
            del g_sum, n_sum, d_eta, d_omega, to_add
        ind = np.argsort([g[5] for g in self.gves]) # Sort in asceding ds 
        self.gves = [self.gves[i] for i in ind]
        
        del ind
        return
    """
    
    def filter_gves(self, ds_tol, eta_tol, omega_tol):
        ds_P, eta_P, omega_P = 5,6,7
        # Groups g-vectors that are close to each other, returns a list of such 'averaged' g-vectors.
        print(f'Grouping {len(self.gves)} g-vectors using tolerance values:')
        print(f'ds_tol = {ds_tol} \t eta_tol = {eta_tol} \t omega_tol = {omega_tol}')
        if not ds_tol or not eta_tol or not omega_tol:
            raise ValueError('Tolerance values must be numeric!')
        if ds_tol<0 or eta_tol<0 or omega_tol<0:
            raise ValueError('Tolerance values must not be negative!')
        
        ds_list = [g[ds_P] for g in self.gves]
        chains_0 = group_to_chains(ds_list, ds_tol)
        
        chains_1 = []
        for C in chains_0:
            eta_list = [self.gves[i][eta_P] for i in C]
            subchains = group_to_chains(eta_list, eta_tol)
            tail_to_head = mod_360(eta_list[subchains[-1][-1]] - eta_list[subchains[0][0]], 0)
            if len(subchains)>1 and abs(tail_to_head)<eta_tol:
                subchains[0] = subchains[-1] + subchains[0]
                del subchains[-1]
            for c in subchains:
                chains_1 = chains_1 + [[C[i] for i in c]]
        
        chains_2 = []
        for C in chains_1:
            omega_list = [self.gves[i][omega_P] for i in C]
            subchains = group_to_chains(omega_list, omega_tol)
            tail_to_head = mod_360(omega_list[subchains[-1][-1]] - omega_list[subchains[0][0]], 0)
            if len(subchains)>1 and abs(tail_to_head)<omega_tol:
                subchains[0] = subchains[-1] + subchains[0]
                del subchains[-1]
            for c in subchains:
                chains_2 = chains_2 + [[C[i] for i in c]]
        
        chains_3 = []
        for C in chains_2:
            eta_list = [self.gves[i][eta_P] for i in C]
            subchains = group_to_chains(eta_list, eta_tol)               
            tail_to_head = mod_360(eta_list[subchains[-1][-1]] - eta_list[subchains[0][0]], 0)
            if len(subchains)>1 and abs(tail_to_head)<eta_tol:
                subchains[0] = subchains[-1] + subchains[0]
                del subchains[-1]
            for c in subchains:
                chains_3 = chains_3 + [[C[i] for i in c]]
 
        chains_4 = []
        for C in chains_3:
            ds_list = [self.gves[i][ds_P] for i in C]
            subchains = group_to_chains(ds_list, ds_tol)
            for c in subchains:
                chains_4 = chains_4 + [[C[i] for i in c]]
        
        av_gvectors = []
        for C in chains_4:
            g_list = [self.gves[i] for i in C]
            av_gvectors.append([*map(mean, zip(*g_list))])
        
        self.gves = [av_gvectors[i] for i in np.argsort([g[ds_P] for g in av_gvectors])] # Sort in asceding ds
        del ds_list, eta_list, omega_list, chains_0, chains_1, chains_2, chains_3, chains_4
        del subchains, tail_to_head, g_list, av_gvectors
        print(f'Resulted in {len(self.gves)} gvectors.')
        return
    
    
    def calc_histo(self):
        I = indexing.indexer()
        I.readgvfile(self.path+self.out_stem+'.gve')
        I.assigntorings()
        hkl_list = [I.unitcell.ringhkls[ds][-1] for ds in I.unitcell.ringds] # [-1] hkl with highest h (usualy last in the set)
        in_rings = np.compress(np.greater(I.ra,-1),np.arange(I.gv.shape[0])) # list of indexed peaks
        omgs = [I.omega[peak] for peak in in_rings]
        etas = [I.eta[peak] for peak in in_rings]
        self.omg_bins =  np.linspace(min(omgs),max(omgs),round(1+(max(omgs)-min(omgs))/self.omega_step), dtype=float)
        self.eta_bins = np.linspace(min(etas),max(etas),round(1+(max(etas)-min(etas))/self.eta_step), dtype=float)
        self.tth_omega = np.zeros( (len(hkl_list),len(self.omg_bins)), dtype=int )
        self.eta_omega = np.zeros( (len(hkl_list),len(self.eta_bins),len(self.omg_bins)), dtype=int )
        for peak in in_rings:
            t = abs(self.omg_bins-I.omega[peak])
            iO = np.where(t == np.amin(t))[0][0]
            t = abs(self.eta_bins-I.eta[peak])
            iE = np.where(t == np.amin(t))[0][0]
            iR = hkl_list.index(I.unitcell.ringhkls[I.unitcell.ringds[I.ra[peak]]][-1])
            if iR and iO:
                self.tth_omega[iR,iO] += 1
                if iE:
                    self.eta_omega[iR,iE,iO] += 1
                    
        filename = self.path+self.out_stem+'_omg_hkl_{:d}x{:d}.raw'.format(self.tth_omega.shape[1], self.tth_omega.shape[0])
        output_file = open(filename, 'wb')
        np.float32(self.tth_omega).tofile(output_file)
        output_file.close()
        filename = self.path+self.out_stem+'_omg_eta_hkl_{:d}x{:d}x{:d}.raw'.format(self.eta_omega.shape[2], self.eta_omega.shape[1], self.eta_omega.shape[0])
        output_file = open(filename, 'wb')
        np.float32(self.eta_omega).tofile(output_file)
        output_file.close()
        
        del I, hkl_list, in_rings, omgs, etas, t, iO, iE, iR, filename, output_file
        return
    
            
    def write_gve_file(self, gve_file = None):
        """ Save gvectors to .gve file. """
        if gve_file:
            self.out_stem = gve_file.replace('.gve','')
        f = open(self.path+self.out_stem+'.gve',"w")
        for line in self.header:
            f.write(line)
        for line in self.dshkls:
            f.write(line)
        f.write(self.col_names_string)                  
        ind = np.argsort([g[5] for g in self.gves]) # Sort in asceding ds 
        self.gves = [self.gves[i] for i in ind]
        for g in self.gves:
            line = ['%f' % x for x in g]
            line[8] = str(int(g[8]))
            f.write(' '.join(line)+' \n')
        f.close()
        
        del ind, line
        return
    
    
    def write_config(self, out_stem):
        """ Create Grainspotter configuration file(.ini). """
        if out_stem:
            self.out_stem = out_stem.replace('.ini','')
        f = open(self.path+self.out_stem+'.ini' ,"w")
        f.write( 'spacegroup {}            !# spacegroup [space group nr]\n'.format(self.spacegroup) )
        if self.ds_ranges:
            for v in self.ds_ranges:
                f.write( 'dsrange {:0.2f} {:0.2f}         !# dsrange [min max], reciprocal d-spacing range, multiple ranges can be specified\n'.format(v[0],v[1]) )
        else:
            f.write( '!dsrange 0.2 0.5          !# dsrange [min max], reciprocal d-spacing range, multiple ranges can be specified\n' )
        if self.tth_ranges:
            for v in self.tth_ranges:
                f.write( 'tthrange {} {}             !# tthrange [min max], multiple ranges can be specified\n'.format(v[0],v[1]) )
        else:
            f.write( '!tthrange 0 30            !# tthrange [min max], multiple ranges can be specified\n' )
        if self.eta_ranges:
            for v in self.eta_ranges:
                f.write( 'etarange {:0.1f} {:0.1f}      !# etarange [min max], full range is [0,360], multiple ranges can be specified\n'.format(v[0],v[1]) )
        else:
            f.write( '!etarange 0 360           !# etarange [min max], full range is [0,360], multiple ranges can be specified\n' )
        if self.omega_ranges:
            for v in self.omega_ranges:
                f.write('omegarange {:0.1f} {:0.1f}   !# omegarange [min max] degrees, full range is [-180,180], multiple ranges can be specified\n'.format(v[0],v[1]) )
        else:
            f.write( '!omegarange -180 180         !# omegarange [min max] degrees, full range is [-180,180], multiple ranges can be specified\n' )
        f.write( 'domega {}                !# domega [stepsize] in omega, degrees\n'.format(self.omega_step) )
        f.write( 'filespecs ./{}.gve ./{}_grains.log !# filespecs [gvecsfile grainsfile]\n'.format(self.out_stem, self.out_stem) )
        f.write( 'cuts {} {} {}            !# cuts [(min_measuments:grain is chosen if there are at least this amount of peaks per grain) (min_completeness: grain is chosen if there are at least this fraction of the expected peaks present) (min_uniqueness: no idea, just put any number)]\n'.format(self.cuts[0], self.cuts[1], self.cuts[2]) )
        f.write( 'eulerstep {}               !# eulerstep [stepsize] : angle step size in Euler space\n'.format(self.eulerstep) )
        f.write( 'uncertainties {} {} {} !# uncertainties [sigma_tth sigma_eta sigma_omega] in degrees\n'.format(self.uncertainties[0],self.uncertainties[1],self.uncertainties[2]) )
        f.write( 'nsigmas {}                 !# nsigmas [Nsig] : maximal deviation in sigmas\n'.format(self.nsigmas) )
        if self.Nhkls_in_indexing:
            f.write( 'Nhkls_in_indexing {}       !# Nhkls_in_indexing [Nfamilies] : use first Nfamilies in indexing\n'.format(self.Nhkls_in_indexing) )
        else:
            f.write( '!Nhkls_in_indexing 15     !# Nhkls_in_indexing [Nfamilies] : use first Nfamilies in indexing\n' )
        if self.random:
            f.write( 'random {}              !# random sampling of orientation space trying 10000 sample points\n'.format(self.random) )
        else:
            f.write( '!random 10000             !# random sampling of orientation space trying 10000 sample points\n' )
        if self.positionfit:
            f.write( 'positionfit               !# fit the position of the grain\n' )
        else:
            f.write( '!positionfit              !# fit the position of the grain\n' )
        if self.minfracg:
            f.write( 'minfracg {}              !# stop search when minfracg (0..1) of the gvectors have been assigned to grains\n'.format(self.minfracg) )
        else:
            f.write( '!minfracg 0.2             !# stop search when minfracg (0..1) of the gvectors have been assigned to grains\n' )
        if self.genhkl:
            f.write( 'genhkl                    !# generate list of reflections\n' )
        else:
            f.write( '!genhkl                   !# generate list of reflections\n' )
        f.close()
        
        return
