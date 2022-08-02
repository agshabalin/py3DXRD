# -*- coding: utf-8 -*-
"""
Created on Wed Feb 02 11:00:00 2022

@author: shabalin

Class to work with a single sweep scan. Loads? process, and peaksearches files.
"""

import sys, os, subprocess, pdb, re
import numpy as np
from numpy import float32
from datetime import datetime
import imageio, fabio
import pickle, yaml
from skimage.transform import warp_polar
from skimage.util import img_as_float
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from hexrd import imageseries
from hexrd.imageseries.omega import OmegaWedges
from ImageD11 import labelimage
#sys.path.insert(1, '/home/shabalin/scripts/')

single_separator = "--------------------------------------------------------------\n"
double_separator = "==============================================================\n"

class SweepScan:
    
    def __init__(self, directory=None, name=None):
        self.log = []
#         self.absorbed =[]
        self.directory = None
        self.name = None
        self.position = []
        self.log_meta = None
        self.fio_meta = None
        self.sweep = {'fio_path':None, 'omega_start':None, 'omega_step':None,
                      'directory':None, 'stem':None, 'ndigits':None, 'ext':None}
        self.chunk = {'frames':[], 'filenames':[], 'omegas':[]}
        self.imgs = None
        self.processing = {'options':None, 'mask':None, 'bckg_indices':None, 'bckg':None}
        self.projs = {'imsum':None, 'immax':None,
                      'q0_pos':[], 'etatth':None, 'omgtth':None, 'omgeta':None}
        self.thresholds = []
        self.pix_tol    = None
        self.spline_file= None
        self.add_to_log('Created SweepScan object.', True)
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
        print(double_separator+'SweepScan object:')
        print('directory:' , self.directory)
        print('name:'      , self.name)
        print('position:'  , self.position)
        
        if self.log_meta:
            print('log_meta:', list(self.log_meta['entries'][0].keys()))
            print('total:', len(self.log_meta['entries']), 'entries with keys:')
            print([k for k in self.log_meta['entries'][0].keys()] )
            
        if self.fio_meta:
            print('fio_meta:', list(self.fio_meta['entries'][0].keys()))
            print('total:', len(self.fio_meta['entries']), 'entries with keys:')
            print([k for k in self.fio_meta['entries'][0].keys()] )
    
    
        print('sweep:'     , self.sweep, '\n')
        print('chunk:'     , self.chunk, '\n')
        
        try: print('imgs:'  , len(self.imgs), 'of images', '\n')
        except: print('imgs:', self.imgs, '\n')
        
        print('processing:' , self.processing, '\n')
        print('projs:'      , self.projs, '\n')
        print('thresholds:' , self.thresholds, '\n')
        print('pix_tol:'    , self.pix_tol, '\n')
        print('spline_file:', self.spline_file, '\n')
#         print('absorbed:', len(self.absorbed))
        if also_log:
            print(single_separator + 'Log:')
            for record in self.log: print(record)
        return
    
    
    def load_sweep(self, omega_start=None, omega_step=None,
                   directory=None, stem=None, ndigits=None, ext=None, frames = None):
        s = self.sweep
        c = self.chunk
        if omega_start is not None: s['omega_start'] = omega_start
        if omega_step  is not None: s['omega_step']  = omega_step
        if directory  : s['directory']   = directory
        if stem       : s['stem']        = stem
        if ndigits    : s['ndigits']     = ndigits
        if ext        : s['ext']         = ext
        if frames     : c['frames']      = frames
        
        print(single_separator+'Loading sweep:', s['directory']+s['stem']+'*'+s['ext'], '...')
        all_files = sorted(os.listdir(s['directory']))
        regex = re.compile(r'\d+')
        check_stem = True
        check_digits = True
        for fname in all_files:
            stem, ext = os.path.splitext(fname)
            if ext in ['.tif', '.cbf', '.edf', '.hdf5', '.h5']:
                if ext != s['ext']:
                    print('Detected file extension \''+ext+'\' vs provided \''+str(s['ext'])+'\'.')
                    x = input('Type: ! to overwrite, k - keep:')
                    if x in ['!']: s['ext'] = ext
                dig_part  = [x for x in regex.findall(stem)][-1]
                dig_len   = len(dig_part)
                stem_part = stem.replace(dig_part, '')
                if check_stem and stem_part != s['stem']:
                    print('Detected file stem \''+stem_part+'\' vs provided \''+str(s['stem'])+'\'.')
                    x = input('Type: ! to overwrite, k - keep, d - keep and don\'t ask anymore:')
                    if x in ['!']: s['stem'] = stem_part
                    if x in ['d']: check_stem = False
                if check_digits and dig_len != s['ndigits']:
                    print(f'Detected {dig_len} digits in file name vs provided '+str(s['ndigits']))
                    if s['ndigits'] == 'auto':
                        s['ndigits'] = dig_len
                    else:                    
                        x = input('Type: ! to overwrite, k - keep, d - keep and don\'t ask anymore:')
                        if x in ['!']: s['ndigits'] = dig_len
                        if x in ['d']: check_digits = False
        
#         if s['ext'] in ['.cbf', '.tif'] and 'p21.2' in self.directory:
#             import cbftiffmxrdfix

        matching_files = [f for f in all_files if (s['stem'] in f and s['ext'] in f)]
        n0 = int(matching_files[0].replace(s['stem'],'').replace(s['ext'],''))
        for i in range(c['frames'][-1]):
            expected_file = s['stem']+str(n0+i).zfill(s['ndigits'])+s['ext']
            if expected_file not in matching_files:
                raise FileNotFoundError(expected_file)
                                  
        c['filenames'] = [matching_files[ii] for ii in c['frames']]
        c['omegas'] = np.empty([len(c['frames']),2], dtype=float)
        for ii in range(len(c['frames'])):
            start_angle = s['omega_start']+s['omega_step']*c['frames'][ii]
            c['omegas'][ii] = [start_angle-0.5*s['omega_step'], start_angle+0.5*s['omega_step']]
        
        if s['omega_step'] > 0:
            imgs_omegas = c['omegas'][::1,::1]
            imgs_frames = c['frames'][::1]
        else: # CRUTCH for imageseries!
            imgs_omegas = c['omegas'][::-1,::-1]
            imgs_frames = c['frames'][::-1]
            
        if s != self.sweep: self.set_attr('sweep', s)
        if c != self.chunk: self.set_attr('chunk', c)
        self.add_to_log('Loading sweep: '+str(s), True)
        self.add_to_log('Selecting chunk: '+str(c), True)
        #np.save(s['directory']+self.name+"omega.npy", imgs_omegas)

        if s['ext'] in ['.tif', '.cbf', 'edf']: # Varex or Pilatus detector.
            config_file = self.directory+str(self.name)+"_sweepparam.yml"
            with open(config_file, 'w') as f:
                f.write('image-files:\n  directory: '+s['directory'])
                f.write('\n  files: "')
                #imgs = []
                for fname in c['filenames']:
                    f.write(fname+'\n')
                    #imgs.append(imageio.imread(s['directory']+fname))
                f.write('\"\noptions:  \n  empty frames: 0\n  max-frames: 0')
                f.write('\nmeta:\n  omega: ')#'\"! load-numpy-array omega.npy\"')
            f.close()
            imgs = imageseries.open(config_file, 'image-files')
            imgs.metadata['omega'] = imgs_omegas

        elif s['ext'] in ['.h5','.hdf5']: # Eiger detector
            config_file = self.directory+str(self.name)+"_sweepparam.yml"
            with open(config_file, 'w') as f:
                #f.write('image-files:\n  directory: '+s['directory'])
                f.write('hdf5:\n  directory: '+s['directory'])
                f.write('\n  files: "')
                #imgs = []
                for fname in c['filenames']:
                    f.write(fname+'\n')
                    #imgs.append(imageio.imread(s['directory']+fname))
                f.write('\"\noptions:  \n  empty frames: 0\n  max-frames: 0')
                f.write('\nmeta:\n  omega: ')#'\"! load-numpy-array omega.npy\"')
            f.close()
            imgs = imageseries.open(config_file, 'image-files')
            imgs.metadata['omega'] = imgs_omegas
    
        self.set_attr('imgs', imgs)
        del imgs_omegas, imgs_frames, all_files, matching_files, config_file
        self.add_to_log(f'{len(imgs)} images of size {imgs[0].shape} loaded.', True)
        return
    
    
    def process_imgs(self, options=None, mask=None, bckg=None):
        p = self.processing
        if options is not None: p['options'] = options
        if mask    is not None: p['mask'   ] = mask
        
        print('Processing sweep with options:', str(p['options']), '...')
        if bckg is None:
            p['bckg_indices'] = []
            p['bckg'] = np.zeros( (self.imgs[0].shape), float)
        elif 'auto' in bckg: # 'auto80' can be used for example.
            try:
                Nmax = min( len(self.imgs), int(bckg.replace('auto', '')) )
                if Nmax<1 or Nmax>200: raise ValueError(bckg+' - Nmax must be >1 and < 200')
            except:
                Nmax = min( len(self.imgs), 50 )
            p['bckg_indices'] = list(range(1,len(self.imgs)-1,round(len(self.imgs)/Nmax) )) # ~Nmax images
            p['bckg'] = calculate_bckg(self.imgs, p['bckg_indices'])
            self.add_to_log(f"Calculated background using {len(p['bckg_indices'])} images", True)
        else:
            p['bckg_indices'] = []
            p['bckg'] = bckg
            
        try   : p['bckg'].shape == self.imgs[0].shape
        except: raise TypeError(type(p['bckg'])+' - incorrect background!')

        ProcessedIS = imageseries.process.ProcessedImageSeries
        if p['options']: imgs = ProcessedIS(self.imgs, [('dark', p['bckg']), p['options']])
        else           : imgs = ProcessedIS(self.imgs, [('dark', p['bckg'])])
        self.set_attr('processing', p)
        self.set_attr('imgs', imgs)
        return
    

    def calculate_projections(self, q0_pos=None):
        if not q0_pos: q0_pos = self.projs['q0_pos']
        im_size = list(self.imgs[0].shape)
        min_r_0 = max(0 - q0_pos[0], q0_pos[0] - im_size[0], 0)
        min_r_1 = max(0 - q0_pos[1], q0_pos[1] - im_size[1], 0)
        max_r_0 = max(q0_pos[0] - 0, im_size[0] - q0_pos[0])
        max_r_1 = max(q0_pos[1] - 0, im_size[1] - q0_pos[1])
        min_rad = np.sqrt(min_r_0**2+min_r_1**2)
        max_rad = np.sqrt(max_r_0**2+max_r_1**2)
        
        print(f'Calculating projections for im_size={im_size}, q0_pos={q0_pos}, min_rad={min_rad}, max_rad={max_rad}\n...')
        immax = imageseries.stats.max(self.imgs, len(self.imgs))
        imsum = self.imgs[0]
        img = warp_polar(img_as_float(self.imgs[0]), center=q0_pos, radius=max_rad)
        etatth = img
        omgeta = np.zeros([len(self.imgs), img.shape[0]], dtype=np.float32)
        omgtth = np.zeros([len(self.imgs), img.shape[1]], dtype=np.float32)
        omgeta[0,:] = etatth.sum(1)
        omgtth[0,:] = etatth.sum(0)

        for i in range(1,len(self.imgs)):
            imsum += self.imgs[i]
            img = warp_polar(img_as_float(self.imgs[i]), center=q0_pos, radius=max_rad)
            etatth += img
            omgeta[i,:] = img.sum(1)
            omgtth[i,:] = img.sum(0)
    
        self.set_attr('projs', {'imsum':img_as_float(imsum), 'immax':immax, 'q0_pos':q0_pos,
                                'etatth':img_as_float(etatth), 'omgtth':omgtth, 'omgeta':omgeta})
        suggested_thr = np.percentile(immax.flatten(), 99.9)
        self.add_to_log(f"99.9 percentage threshold: {round(suggested_thr)}", True)
        del max_rad, img, etatth, omgeta, omgtth, suggested_thr
        return

    
    def plot(self):        
        fig = plt.figure(figsize=(16, 10))

        sub = fig.add_subplot(221) # instead of plt.subplot(2, 2, 1)
        im = sub.imshow(np.log10(self.projs['imsum']), interpolation='None')
        plt.title('Sum of images')
        divider = make_axes_locatable(sub)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

        sub = fig.add_subplot(222)
        im = sub.imshow(np.log10(self.projs['etatth']), interpolation='None')
        plt.title('eta_tth projection')
        divider = make_axes_locatable(sub)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

        sub = fig.add_subplot(223)
        im = sub.imshow(np.log10(self.projs['omgeta']), interpolation='None')
        plt.title('omg_eta projection')
        divider = make_axes_locatable(sub)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

        sub = fig.add_subplot(224)
        im = sub.imshow(np.log10(self.projs['omgtth']), interpolation='None')
        plt.title('omg_tth projection')
        divider = make_axes_locatable(sub)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')
        
        plt.show()
        self.add_to_log('Writing file: '+self.directory+self.name+'_projs.png', True)
        fig.savefig(self.directory+self.name+'_proj.png')
        return
    

    def calculate_thresholds(self):
        a = self.projs['immax'].flatten()
        t4000 = 0.5*np.percentile(a, 99.69)+0.5*np.percentile(a, 99.72)
        t2000 = 0.5*np.percentile(a, 99.37)+0.5*np.percentile(a, 99.38)
        t1000 = 0.5*np.percentile(a, 98.88)+0.5*np.percentile(a, 98.88)
        t0500 = 0.5*np.percentile(a, 98.17)+0.5*np.percentile(a, 98.16)
        thr = 100*np.round(0.01*(0.25*t4000 + 0.5*t2000 + 1*t1000 + 2*t0500)/4)
        return [int(t) for t in [thr, 2*thr, 4*thr, 8*thr]]

                
    def export_data(self, thr=None):
        # Export data using a certain threshold.
        self.add_to_log("Exporting data as: "+self.directory+self.name, True)
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            self.add_to_log('Created directory: '+self.directory, True)
        
        if not thr: thr = min(self.thresholds)
    
        path = self.directory+self.name
        
        try:
            self.add_to_log('Writing file: '+path+f"_t{thr}.npz", True)
            imageseries.write(self.imgs,"dummy","frame-cache",
                              cache_file=path+f"_t{thr}.npz", threshold=thr)
        except: self.add_to_log('Failed to write .npz file!', True)
        
        try:
            self.add_to_log('Writing file: '+path+"_bckg.tif", True)
            imageio.imwrite(path+"_bckg.tif"       , self.processing['bckg'])
        except: self.add_to_log('Failed to write *_bckg.tif file!', True)
        
        try:
            self.add_to_log('Writing file: '+path+"_proj_immax.tif", True)
            imageio.imwrite(path+"_proj_immax.tif" , self.projs['immax'])
        except: self.add_to_log('Failed to write *_proj_immax.tif file!', True)
        
        try:        
            self.add_to_log('Writing file: '+path+"_proj_imsum.tif", True)
            imageio.imwrite(path+"_proj_imsum.tif" , self.projs['imsum'])
        except: self.add_to_log('Failed to write *_proj_imsum.tif file!', True)

        try:            
            self.add_to_log('Writing file: '+path+"_proj_etatth.tif", True)
            imageio.imwrite(path+"_proj_etatth.tif", self.projs['etatth'])
        except: self.add_to_log('Failed to write *_proj_etatth.tif file!', True)

        try:            
            self.add_to_log('Writing file: '+path+"_proj_omgtth.tif", True)
            imageio.imwrite(path+"_proj_omgtth.tif", self.projs['omgtth'])
        except: self.add_to_log('Failed to write *_proj_omgtth.tif file!', True)

        try:            
            self.add_to_log('Writing file: '+path+"_proj_omgeta.tif", True)
            imageio.imwrite(path+"_proj_omgeta.tif", self.projs['omgeta'])
        except: self.add_to_log('Failed to write *_proj_omgeta.tif file!', True)

#         try:            
#             self.add_to_log('Writing file: '+path+"_omegas.npy", True)
#             np.save(path+"_omegas.npy", self.chunk['omegas'])
#         except: self.add_to_log('Failed to write *_omegas.npy file!', True)
        
        return


    def file_series_object(self):
        import fabio
        order = list(range(len(self.imgs)))
        def frm(i):
            omega = 0.5*self.chunk['omegas'][0][0] + 0.5*self.chunk['omegas'][0][1]
            f = fabio.fabioimage.fabioimage( self.imgs[i], {'Omega': omega} )
            f.filename = self.chunk['filenames'][i]
            f.currentframe = i
            return f

        yield( frm(order[0]) ) # first
        for i in order: yield frm(i)
    
    
    def peaksearch(self, thresholds = [], spline_path = None):
        if thresholds: self.set_attr('thresholds', thresholds)
        self.add_to_log(f'Running peaksearch on {len(self.imgs)} images using {len(self.thresholds)} thresholds...', True)
        
        if spline_path:
            os.system('cp '+spline_path+' '+self.directory+self.name+'.spline') # Copy spline file here
            self.set_attr('spline_file', self.name+'.spline')
        elif self.spline_file and not os.path.exists(self.directory+self.spline_file):
            os.system('cp '+self.spline_file+' '+self.directory+self.name+'.spline') # Probably it is some other directory so need to be copied here
            self.set_attr('spline_file', self.name+'.spline')
        else: self.set_attr('spline_file', None)
        
#         def one_thr(thr):
#             peaksfile = open(self.directory+self.name+f'_peaks_t{thr}.flt', 'w')
#             lio = labelimage.labelimage(self.imgs[0].shape, peaksfile)
#             for i_img, img in enumerate(self.imgs):
#                 omg = np.mean(self.chunk['omegas'][i_img])
#                 lio.peaksearch(img, thr, omg)
#                 lio.mergelast()
#             lio.finalise()
#             peaksfile.close()
#             return
        
#         from multiprocessing import Pool
#         with Pool(5) as p:
#             p.map(one_thr, self.thresholds)
#         for r in res:
#             self.add_to_log(r, True)
        for i_thr, thr in enumerate(self.thresholds):
            print(f'Threshold {thr}, saving results to:')
            peaksfile = open(self.directory+self.name+f'_peaks_t{thr}.flt', 'w')
            lio = labelimage.labelimage(self.imgs[0].shape, peaksfile)
            for i_img, img in enumerate(self.imgs):
                omg = np.mean(self.chunk['omegas'][i_img])
                lio.peaksearch(img, thr, omg)
                lio.mergelast()
            lio.finalise()
            peaksfile.close()
            self.add_to_log(self.directory+self.name+f'_peaks_t{thr}.flt', True)
        
        if self.spline_file:
            for thr in self.thresholds:
                flt_file = self.directory+self.name+f'_peaks_t{thr}.flt'
                apply_spline_to_fltfile(flt_file, flt_file, self.spline_file, self.projs['immax'].shape[0])
                self.add_to_log(f"For {flt_file} and threshold {thr} applied spline file:", self.spline)

        return


    def save_peaksearch_yaml(self, thresholds='auto', pix_tol=None, spline_path=None): # Configuration file for peaksearch
        
        if not thresholds:
            thresholds = self.thresholds
        elif 'auto' in thresholds:
            thresholds = self.calculate_thresholds()
        
        try:
            suggested_thr = np.percentile(self.projs['immax'].flatten(), 99.9)
            if thresholds[0] > suggested_thr:
                print('WARNING! Provided thresholds:', thresholds, f'are higher than 99.9 percentile: {suggested_thr}!')
        except:
            pass
    
        if not pix_tol: pix_tol = self.pix_tol
        
        if spline_path:
            os.system('cp '+spline_path+' '+self.directory+self.name+'.spline') # Copy spline file here
            self.set_attr('spline_file', self.name+'.spline')
        elif self.spline_file and not os.path.exists(self.directory+self.spline_file):
            os.system('cp '+self.spline_file+' '+self.directory+self.name+'.spline') # Probably it is some other directory so need to be copied here
            self.set_attr('spline_file', self.name+'.spline')
        else: self.set_attr('spline_file', None)
        
        self.add_to_log('Writing file: '+self.directory+self.name+'_peaksearch.yaml', True)
        with open(self.directory+self.name+'_peaksearch.yaml', 'w') as f:
            f.write('# parameters for peaksearch')
            f.write('\nimage_dir: {}'.format(self.sweep['directory']))
            f.write('\nimage_ext: {}'.format(self.sweep['ext']))
            f.write('\nimage_stem: {} #Stem of images'.format(self.sweep['stem']))
            f.write('\nndigits: {}'.format(self.sweep['ndigits']))
            dig_part = self.chunk['filenames'][0].replace(self.sweep['stem'],'').replace(self.sweep['ext'],'')
            f.write('\nfirst_image: {}  # Index of first image'.format( int(dig_part) ))
            f.write('\nnbr_images: {}  # Number of images to peaksearch'.format(len(self.chunk['omegas'])))
            f.write('\nomegastep: {}'.format(self.sweep['omega_step']))
            f.write('\nstartomega: {}'.format(0.5*self.chunk['omegas'][0][0]+0.5*self.chunk['omegas'][0][1])) # Use middle omega.
            f.write('\ndark_image: {}'.format(self.directory+self.name+'_bckg.tif'))
            f.write('\nthresholds: {}'.format(thresholds))
            if spline_path:
                f.write('\nspline: {}'.format(self.directory+self.spline_file))
            f.write('\noutput_dir: {}'.format(self.directory))
            f.write('\nstem_out: {}'.format(self.name))
            f.write('\n#kwargs: \'--OmegaOverride\'#additional keyword aguments for peaksearch.py as one string')
            f.write('\n# parameters for merge_flt')
            if pix_tol:
                f.write('\npixel_tol: {} #minimum distance between peaks (pixels)'.format(pix_tol))
            f.write('\nmerged_name: \'{}\''.format(self.name+'.flt'))
        f.close()
        return self.directory+self.name+'_peaksearch.yaml'

    
    def run_peaksearch(self, yaml_file='auto', use_imgs=False, use_temp_tifs=False, del_temp_tifs=True): # Wrapper for the ImageD11 peaksearch.py script.
        if yaml_file == 'auto': yaml_file = self.save_peaksearch_yaml()
        with open(yaml_file) as f: pars = yaml.safe_load(f)
        if pars['stem_out'] == None: pars['stem_out'] = ''
        
        first_im = int(pars['first_image'])
        last_im  = int(pars['first_image']) + int(pars['nbr_images']) - 1
        path_out = os.path.join(pars['output_dir'], pars['stem_out']+'_peaks')
        if use_temp_tifs:
            path_inp = os.path.join(self.directory+self.name+"_temp/", pars['image_stem'])
            pars['image_ext'] = '.tif'
            pars['dark_image'] = None
            
            if os.path.exists(self.directory+self.name+"_temp/"):
                subprocess.call('rm '+self.directory+self.name+"_temp/*.tif", shell=True) # Delete old tif files
            else:
                os.makedirs(self.directory+self.name+"_temp/")
            
            self.add_to_log(f"Saving temporary images in "+self.directory+self.name+"_temp/", True)
            for ind, img in enumerate(self.imgs):
                dig_part = str(first_im+ind).zfill(pars['ndigits'])
                imageio.imwrite(path_inp+dig_part+pars['image_ext'], img)
        else:
            path_inp = os.path.join(pars['image_dir'], pars['image_stem'])

        # construct the command for peaksearch.py
        command = 'peaksearch.py -o {} -n {} '.format(path_out, path_inp)
        command+= '-F {} --ndigits {:d} '.format(pars['image_ext'], pars['ndigits'])
        command+= '-f {:d} -l {:d} '.format(first_im, last_im)
        command+= '-S {:.3f} -T {:.3f} -p Y'.format(pars['omegastep'], pars['startomega'])
        if pars['dark_image']: command+= ' -d {}'.format(pars['dark_image']) 
        for t in pars['thresholds']: command += ' -t {:d}'.format(t)
        if 'kwargs' in pars: command += ' {}'.format(pars['kwargs'])
        
        self.add_to_log('Running peaksearch in: '+self.directory+'\n'+command, True)
        if use_imgs:
            import time
            reallystart = time.time()
            try:
                from argparse import ArgumentParser
                from py3DXRD import peaksearcher
                parser = ArgumentParser()
                myparser = peaksearcher.get_options(parser)
                options , args = myparser.parse_known_args(command.split())
                options.file_series_object =  self.file_series_object()
                peaksearcher.peaksearch_driver(options, args)
            except:
                if myparser is not None:
                    myparser.print_help()
                print("\n\n And here is the problem:\n")
                raise

            end = time.time()
            t = end-reallystart
            print("Total time = %f /s" % ( t ))
        else:
            process = subprocess.run(command.split(), check=True,
                                     stdout=subprocess.PIPE, universal_newlines=True)
            self.add_to_log('Output:'+process.stdout, False)
            print('Last line in the output:'+process.stdout.splitlines()[-1])
        self.set_attr('thresholds', pars['thresholds'])
        
        if 'spline' in pars.keys():
            for t in pars['thresholds']:
                flt_file = pars['output_dir']+pars['stem_out']+f'_peaks_t{t}.flt'
                apply_spline_to_fltfile(flt_file, flt_file, pars['spline'], self.projs['immax'].shape[0])
                self.add_to_log(f"For {flt_file} and threshold {t} applied spline file:", pars['spline'])
                     
        if use_temp_tifs and del_temp_tifs:
            subprocess.call('rm -r '+self.directory+self.name+"_temp/", shell=True)
            self.add_to_log(f"Deleted temporary images.", True)
        return

    
    def merge_peaks(self, par_file, yaml_file='auto'): # Wrapper for ImageD11 merge_flt.py
        if yaml_file == 'auto': yaml_file = self.save_peaksearch_yaml()
        with open(yaml_file) as f: pars = yaml.safe_load(f)
        if pars['stem_out'] == None: pars['stem_out'] = ''

        if 'merged_name' in pars:
            file_out = os.path.join(pars['output_dir'],pars['merged_name'])
        else:
            file_out = os.path.join(pars['output_dir'],pars['stem_out']+'_peaks_merged.flt')
        inp = os.path.join(pars['output_dir'], pars['stem_out']+'_peaks')
        # construct the command for merge_flt.py
        command = 'merge_flt.py {} {} {} {:d} '.format(par_file,inp,file_out,pars['pixel_tol'])
        command+= ('{:d} '*len(pars['thresholds'])).format(*pars['thresholds'])
        self.add_to_log('Merging flt files matching: '+inp+'\n'+command, True)
        process=subprocess.run(command.split(), check=True,
                               stdout=subprocess.PIPE, universal_newlines=True)
        self.add_to_log('Output:'+process.stdout, False)
        print('Last line in the output:'+process.stdout.splitlines()[-1])
        self.set_attr('thresholds', pars['thresholds'])
        self.set_attr('pix_tol', pars['pixel_tol'])
        return


    def save_tifs(self): # Save images as tif files.
        path = self.directory+self.name+"_tifs/"
        self.add_to_log(f"Saving tifs to: "+path, True)
        if not os.path.exists(path):
            os.makedirs(path)
            self.add_to_log('Created directory: '+path, True)
        for i in range(0,len(self.imgs)):
            imageio.imwrite(path+"{:04d}.tif".format(i), np.float32(self.imgs[i]))
        self.add_to_log(f"Saved {len(self.imgs)} .tif files", True)        
        return
    
    
    def delete_tifs(self): # Save images as tif files.
        command = 'rm -r '+self.directory+self.name+"_tifs/"
        subprocess.call(command, shell=True)
        self.add_to_log(f"Deleted temporary .tif files", True)
        return
    
                        
def calculate_bckg(imgs, indices): # Uses median of a subset of images
    try: sub_set = np.asarray( [imgs[i] for i in indices] )
    except: raise ValueError(' - incorrect indices!')
    return np.median(sub_set,axis=0).astype(float32)


def apply_spline_to_fltfile(flt_file_in, flt_file_out, spline, sc_dim):
    from ImageD11 import columnfile, blobcorrector
    if spline == 'perfect':
        cor = blobcorrector.perfect()
    else:
        cor = blobcorrector.correctorclass( spline )

    inc = columnfile.columnfile( flt_file_in )
    inc.s_raw = sc_dim - inc.s_raw
    for i in range( inc.nrows ):
        inc.sc[i], inc.fc[i] = cor.correct( inc.s_raw[i], inc.f_raw[i] )
    inc.sc = sc_dim - inc.sc
    inc.writefile( flt_file_out )

#### TESTS    ###############
# path = '/asap3/petra3/gpfs/p21.2/2021/data/11008399/'
# omega_step = 0.25
# omega_start = 65+0.5*omega_step
# iL = 1
# det_num = 4
# path_inp = path+'raw/CS_1/load{:d}/supersweep/Varex_{:d}/'.format(iL, det_num) # path to raw data
# stem = 'sweep_'
# iC = 0
# frames = list(range(iC+1,iC+559))

# G = SweepScan('/home/shabalin/11008399/test/','V4_test')
# G.load_sweep(omega_start=omega_start, omega_step=omega_step,
#                     directory=path_inp, stem='sweep_', ndigits=5, ext='.tif',
#                     frames=frames)
# G.process_imgs(bckg_indices=list(range(10,550,10)) )
# G.calculate_projections(q0_pos=[1441.2,1444.3])
# G.save_peaksearch_yaml(thresholds=[4000,8000,16000,32000], pix_tol=10)
# G.export_data()
# G.plot()
# G.run_peaksearch(use_temp_files =True)
# G.merge_peaks(par_file = '/home/shabalin/11008399/test/V4_iconfig_Ti_fitted.par')
# G.print(True)