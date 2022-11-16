# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 11:30:49 2021

@author: shabalin

Class to load and process raw data in the sweep scans.
"""


import os, pdb
import numpy as np
import imageio, fabio
from numpy import float32
import pickle, yaml, subprocess
from skimage.transform import warp_polar
from skimage.util import img_as_float
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from hexrd import imageseries
from hexrd.imageseries.omega import OmegaWedges
#from .tth_eta_projections import tth_eta_projections
#from .calc_projections import calc_projections


class SweepScan:

    
    def __init__(self, path_out, det_code, omega_zero, omega_step):
        """ Initialze class for a single 3D XRD sweep scan. """
        self.det_code = det_code
        self.omega_zero = omega_zero
        self.omega_step = omega_step
        self.path_inp = ''
        self.path_out = path_out
        self.stem_inp = ''
        self.stem_out = ''
        self.frames = None
        self.file_names = None
        self.file_ext = None
        self.omega = None
        self.comment = ''
        self.imgs = []
        self.mask = None
        self.proj_max = None
        self.imgs_sum = None
        self.eta_tth = None
        self.omg_eta = None
        self.omg_tth = None
        self.tth_omega = None
        self.eta_omega = None
        self.bckg_image = None
        self.pars = None
        self.thr = None
        
        return
        

    def load_data(self, path_inp, stem_inp, file_ext, frames):
        """ Load data."""
        self.path_inp = path_inp
        self.stem_inp = stem_inp
        self.frames = frames
        self.file_ext = file_ext
        if 'cbf' in self.file_ext: # '.tif', '.cbf', '.edf', '.hf' for example
            import cbftiffmxrdfix
        
        all_files = sorted(os.listdir(self.path_inp))
        selected_files = [f for f in all_files if (self.stem_inp in f and self.file_ext in f)]
        self.file_names = [selected_files[ii] for ii in self.frames]
            
        self.omega = np.empty([len(self.frames),2], dtype=float)
        for ii in range(0,len(self.frames)):
            start_angle = self.omega_zero+self.omega_step*ii #self.frames[ii]
            self.omega[ii] = [start_angle, start_angle+self.omega_step]
        
        if self.omega_step > 0:
            imgs_omega = self.omega[::1,::1]
            imgs_frames = self.frames[::1]
        else: # CRUTCH for imageseries!
            imgs_omega = self.omega[::-1,::-1]
            imgs_frames = self.frames[::-1]
        #np.save("./omega.npy", imgs_omega)
        
        if self.det_code[0] in ['V', 'P']: # Varex or Pilatus detector.
            config_file = self.path_out+self.stem_out+self.det_code+"_param.yml"
            with open(config_file, 'w') as f:
                f.write('image-files:\n  directory: '+self.path_inp)
                f.write('\n  files: "')
                #self.imgs = []
                for fileName in [selected_files[ii] for ii in imgs_frames]:
                    f.write(fileName+'\n')
                    #self.imgs.append(imageio.imread(self.path_inp+fileName))
                f.write('\"\noptions:  \n  empty frames: 0\n  max-frames: 0')
                f.write('\nmeta:\n  omega: ')#'\"! load-numpy-array omega.npy\"')
            self.imgs = imageseries.open(config_file, 'image-files')
            self.imgs.metadata['omega'] = imgs_omega

        elif self.det_code[0]=='E': # E - Eiger detector
            config_file = self.path_out+self.stem_out+self.det_code+"_param.yml"
            with open(config_file, 'w') as f:
                #f.write('image-files:\n  directory: '+self.path_inp)
                f.write('hdf5:\n  directory: '+self.path_inp)
                f.write('\n  files: "')
                #self.imgs = []
                for fileName in [selected_files[ii] for ii in imgs_frames]:
                    f.write(fileName+'\n')
                    #self.imgs.append(imageio.imread(self.path_inp+fileName))
                f.write('\"\noptions:  \n  empty frames: 0\n  max-frames: 0')
                f.write('\nmeta:\n  omega: ')#'\"! load-numpy-array omega.npy\"')
            self.imgs = imageseries.open(config_file, 'image-files')
            self.imgs.metadata['omega'] = imgs_omega
            
        elif self.det_code[0]=='L': # L - Lambda detector
            raise TypeError(self.det_code + ' - detector is not set up yet!')
            
        else:
            raise TypeError(self.det_code + ' - detector is not set up yet!')
        
        del imgs_omega, imgs_frames, all_files, selected_files, config_file
        return
            
           
    def calculate_background(self, indices):
        """ Calculate background using a median of a subset of images. """
        try:
            sub_set = np.asarray( [self.imgs[i] for i in indices] )
            self.bckg_image = np.median(sub_set,axis=0).astype(float32)
            del sub_set
        except:
            raise ValueError(' - incorrect indices!')
            
        return

    
    def process_data(self, config_file_par, processing_options = None):
        # Process data.
        with open(config_file_par) as f:
            for line in f:
                if ('y_center' in line):
                    y_cen = float(line.split()[1])
                elif ('z_center' in line):
                    z_cen = float(line.split()[1])        
        f.close()
                
        try:
            if self.bckg_image.shape == self.imgs[0].shape:
                pass
        except:
            raise TypeError(type(self.bckg_image)+' - incorrect background!')

        ProcessedIS = imageseries.process.ProcessedImageSeries
        if processing_options==None:
            ops = [('dark', self.bckg_image)] #('flip', 'h')]
        else:
            ops = [('dark', self.bckg_image), processing_options]
        self.imgs = ProcessedIS(self.imgs, ops)
        self.proj_max = imageseries.stats.max(self.imgs, len(self.imgs))
#        self.tth_omega, self.eta_omega = tth_eta_projections(self.imgs, [y_cen, z_cen], self.imgs.shape[0]/2, 360)       
        
        cenpos = [z_cen, y_cen]
        maxrad = max( list(np.asarray(self.imgs[0].shape)-np.asarray(cenpos)) + cenpos )
        self.imgs_sum = self.imgs[0]
        img = warp_polar(img_as_float(self.imgs[0]), center=cenpos, radius=maxrad)
        self.eta_tth = img
        self.omg_eta = np.zeros([len(self.imgs), img.shape[0]], dtype=np.float32)
        self.omg_tth = np.zeros([len(self.imgs), img.shape[1]], dtype=np.float32)
        self.omg_eta[0,:] = self.eta_tth.sum(1)
        self.omg_tth[0,:] = self.eta_tth.sum(0)

        for i in range(1,len(self.imgs)):
            self.imgs_sum += self.imgs[i]
            img = warp_polar(img_as_float(self.imgs[i]), center=cenpos, radius=maxrad)
            self.eta_tth += img
            self.omg_eta[i,:] = img.sum(1)
            self.omg_tth[i,:] = img.sum(0)
    
        suggested_thr = np.percentile(self.proj_max.flatten(), 99.9)
        print(f"99.9 percentage threshold: {round(suggested_thr)}" )
        
        del y_cen, z_cen, ops, cenpos, maxrad, img, suggested_thr
        return

    
    def export_data(self, stem_out, thr=0):
        """ Export data using a certain threshold. """
        self.stem_out = stem_out
        self.thr = thr
        if not os.path.exists(self.path_out):
            os.makedirs(self.path_out)
        path = self.path_out+self.stem_out+self.det_code
        imageseries.write(self.imgs,"dummy","frame-cache",
                          cache_file=path+f"_t{self.thr}.npz",threshold=thr)
        imageio.imwrite(path+"_bckg.tif", self.bckg_image)
        imageio.imwrite(path+"_proj_max.tif", self.proj_max)
        np.save(path+"_omega.npy", self.omega)
        d = self.__dict__ # all the attributes of the class
        metadata = {k : d[k] for k in d if k not in ('omega','imgs',
                                                     'proj_max','bckg_image')}
        with open(path+"_metadata.yaml", 'w') as f:
            yaml.dump(metadata, f)
        
        del path, d, metadata
        return
    
    def save_tifs(self):
        """ Save images as tif files. """
        path = self.path_out+self.stem_out+self.det_code+"_tifs/"
        if not os.path.exists(path):
            os.makedirs(path)
        for i in range(0,len(self.imgs)):
            imageio.imwrite(path+"{:04d}.tif".format(i), self.imgs[i])
        
        del path
        return


    def write_peaksearch_yaml(self, thrs, min_peak_distance):
        """ Create configuration file for peaksearch: """
        genpath = self.path_out+self.stem_out+self.det_code
        with open(genpath+"_peaksearch.yaml", 'w') as f:
            f.write('# parameters for peaksearch')
            f.write(f'\ndet_code: {self.det_code}')
            f.write(f'\nimage_path: {self.path_inp}')
            f.write(f'\nimage_stem: {self.stem_inp} #Stem of images')
            filename = self.file_names[0].replace(self.file_ext,'')
            dig_part = filename.replace(self.stem_inp,'')
            f.write(f'\nfirst_image: {str(int(dig_part))}  # Index of first image')
            n = len(self.file_names)
            f.write(f'\nnbr_images: {n}  # Number of images to peaksearch')
            f.write(f'\nfiletype: {self.file_ext}')
            f.write(f'\noutput_dir: {self.path_out}')
            f.write(f'\nstem_out: {self.stem_out}')
            f.write(f'\ndark_image: {genpath}_bckg.tif')
            f.write(f'\nomegastep: {self.omega_step}')
            f.write(f'\nstartomega: {self.omega[0][0]}')
            f.write(f'\nthresholds: {str(thrs)}')
            f.write('\n#kwargs: \'--OmegaOverride\'#additional keyword aguments for peaksearch.py as one string')
            f.write(f'\nndigits: {len(dig_part)}')
            f.write('\n# parameters for merge_flt')
            f.write(f'\nmerged_name: \'{self.det_code}_peaks_merged.flt\'')
            f.write(f'\npixel_tol: {min_peak_distance} #minimum distance between peaks (pixels)')
            
        del genpath, filename, dig_part, n
        return
    
        
    def plot(self):        
        fig = plt.figure(figsize=(16, 10))

        sub = fig.add_subplot(221) # instead of plt.subplot(2, 2, 1)
        im = sub.imshow(np.log10(self.imgs_sum), interpolation='None')
        plt.title('Sum of images')
        divider = make_axes_locatable(sub)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

        sub = fig.add_subplot(222)
        im = sub.imshow(np.log10(self.eta_tth), interpolation='None')
        plt.title('eta_tth projection')
        divider = make_axes_locatable(sub)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

        sub = fig.add_subplot(223)
        im = sub.imshow(np.log10(self.omg_eta), interpolation='None')
        plt.title('omg_eta projection')
        divider = make_axes_locatable(sub)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

        sub = fig.add_subplot(224)
        im = sub.imshow(np.log10(self.omg_tth), interpolation='None')
        plt.title('omg_tth projection')
        divider = make_axes_locatable(sub)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')
        
        plt.show()
        fig.savefig(self.path_out+self.stem_out+self.det_code+'.png')
        
        del fig, sub, im, divider, cax
        return