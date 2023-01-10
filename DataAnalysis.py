# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:00:00 2022

@author: shabalin

Class to perform full 3D XRD analysis
"""

import sys
import os
import subprocess
import pickle
import tifffile
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import py3DXRD
from py3DXRD.p212_tools import load_p212_log, load_p212_fio, parse_p212_command
from py3DXRD.angles_and_ranges import convert_ranges, mod_360, disorientation

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
        self.position = [0, 0, 0]
        self.rotation = [0, 0, 0]
        self.sweepProcessors = []
        self.peakIndexers = []
        self.yml_det_order = []
        self.gvectorEvaluator = None
        self.grainSpotter = None
        self.grains = []
        self.absorbed = []
        self.pos_tol = None
        self.ang_tol = None
        self.add_to_log('Initialized DataAnalysis object.', True)
        if directory:
            self.set_attr('directory', directory)
        if name:
            self.set_attr('name', name)
        return

    def add_to_log(self, str_to_add, also_print=False):
        self.log.append(str(datetime.now()) + '> ' + str_to_add)
        if also_print:
            print(str_to_add)
        return

    def set_attr(self, attr, value):
        old = getattr(self, attr)
        setattr(self, attr, value)
        new = getattr(self, attr)
        self.add_to_log(attr + ': ' + str(old) + ' -> ' + str(new))
        return

    def add_to_attr(self, attr, value):
        old_list = getattr(self, attr)
        setattr(self, attr, old_list + [value])
        new_list = getattr(self, attr)
        self.add_to_log(attr + ': += ' + str(new_list[-1]))
        return

    def print(self, also_log=False):
        print(double_separator + 'DataAnalysis object:')
        print('directory:', self.directory)
        print('name:', self.name)
        print('absorbed:', len(self.absorbed))
        print('material:', self.name)
        print('pressure:', self.pressure)
        print('temperature:', self.temperature)
        print('position:', self.position)
        print('rotation:', self.rotation)
        print('sweepProcessors:', len(self.sweepProcessors))
        print('peakIndexers:', len(self.peakIndexers))
        print('yml_det_order:', self.yml_det_order)
        print('gvectorEvaluator:', self.gvectorEvaluator)
        print('grainSpotter:', self.grainSpotter)
        print('grains:', len(self.grains))
        print('pos_tol:', self.pos_tol)
        print('ang_tol:', self.ang_tol)
        if also_log:
            print(single_separator + 'Log:')
            for record in self.log:
                print(record)
        return

    def save_geometries_as_yml(self, yml_det_order=None):
        if yml_det_order:
            self.set_attr('yml_det_order', yml_det_order)
        gs = [self.sweepProcessors[ind - 1].geometry for ind in self.yml_det_order]
        yml_file = py3DXRD.Geometry.save_geometries_as_yml(
            gs, self.directory, self.name + '.yml', overwrite=True)
        self.add_to_log(
            'Exported all the geometries as .yml file: ' +
            yml_file,
            True)
        return

    def process_images(
            self,
            frames=None,
            save_tifs=False,
            q0_pos=None,
            rad_ranges=None,
            thr=None):
        for SP in self.sweepProcessors:
            print(single_separator + f'\nPROCESSING for {SP.name}:')
            if not os.path.exists(SP.directory):
                # Create the directory if does not exist.
                os.makedirs(SP.directory)

            # Loading and processing images:
            if type(frames) == list:
                SP.chunk['frames'] = frames  # If substack is needed.
            # Loading and processing images usually takes SOME time.
            SP.load_sweep()
            if save_tifs:
                SP.save_tifs()  # Usually not needed.
            # bckg = 'auto50' - median over equally spaced images (not more
            # than ~50).
            SP.process_imgs(bckg='auto50')

            # Usually takes LONG time.
            SP.calculate_projections(q0_pos=q0_pos, rad_ranges=rad_ranges)
            SP.plot()
            if not thr or thr == 'auto':
                SP.calculate_thresholds()
                thr = min(SP.thresholds)
            SP.export_data(thr=thr)

            imgs = SP.imgs
            SP.set_attr('imgs', None)  # As images are too heavy to save
            pickle.dump(
                SP,
                open(
                    SP.directory +
                    SP.name +
                    "_SweepProcessor.p",
                    "wb"))
            SP.set_attr('imgs', imgs)

    def peaksearch(
            self,
            peaksearch_thresholds,
            peakmerge_thresholds,
            min_peak_dist,
            use_temp_tifs=False,
            del_temp_tifs=False):
        for SP in self.sweepProcessors:
            print(single_separator + f'\nPEAKSEARCHING for {SP.name}:')
            if not peaksearch_thresholds or peaksearch_thresholds == 'auto':
                SP.calculate_thresholds()
                peaksearch_thresholds = SP.thresholds
            if not peakmerge_thresholds or peakmerge_thresholds == 'auto':
                peakmerge_thresholds = peaksearch_thresholds
            yaml_file = SP.save_peaksearch_yaml(
                thresholds=peaksearch_thresholds,
                spline_path=SP.geometry.spline_file)
            SP.run_peaksearch(
                yaml_file=yaml_file,
                use_imgs=False,
                use_temp_tifs=use_temp_tifs,
                del_temp_tifs=del_temp_tifs)  # Usually takes LONG time.
            for thr in peaksearch_thresholds:
                PI = py3DXRD.PeakIndexer(directory=SP.directory)
                PI.load_flt(flt_file=SP.name + f'_peaks_t{thr}.flt')
                PI.save_flt(
                    flt_file=SP.name +
                    f'_peaks_t{thr}_backup.flt',
                    overwrite=True)
                if thr < max(peaksearch_thresholds):
                    PI.remove_not_in_range(int_range=[4 * thr, 8 * thr])
                else:
                    PI.remove_not_in_range(int_range=[4 * thr, 8000 * thr])
                PI.save_flt(
                    flt_file=SP.name +
                    f'_peaks_t{thr}.flt',
                    overwrite=True)
                if len(PI.peaks) < 1 and thr in peakmerge_thresholds:
                    peakmerge_thresholds.remove(thr)

            yaml_file = SP.save_peaksearch_yaml(
                thresholds=peakmerge_thresholds, pix_tol=min_peak_dist)
            SP.merge_peaks(yaml_file=yaml_file)
            # As images are too heavy and not needed futhermore.
            SP.set_attr('imgs', None)
            pickle.dump(
                SP,
                open(
                    SP.directory +
                    SP.name +
                    "_SweepProcessor.p",
                    "wb"))
            SP.projs = None  # To save memory
            SP.processing['bckg'] = None  # To save memory
            SP.processing['mask'] = None  # To save memory

    def index(self, move_det_xyz_mm=[0, 0, 0]):
        if not move_det_xyz_mm:
            move_det_xyz_mm = [0, 0, 0]
        self.set_attr('peakIndexers', [])
        for iP, SP in enumerate(
                self.sweepProcessors):    # (e.g. Varex 1,2,3,4).
            print(single_separator + f'\nINDEXING of {SP.name}:')
            # Evaluating peaks:
            PI = SP.generate_PeakIndexer()
            if move_det_xyz_mm != [0, 0, 0]:
                PI.geometry.move_detector(move_det_xyz_mm)
            PI.run_indexer(
                gve_file=PI.name +
                '_' +
                self.material['name'] +
                '.gve')
            PI.set_attr('name', PI.name + '_' + self.material['name'])
            if iP == 0:
                GE = PI.generate_GvectorEvaluator()
                GE.set_attr('name', self.name)
            else:
                GE.absorb(PI.generate_GvectorEvaluator())

            pickle.dump(
                PI,
                open(
                    PI.directory +
                    PI.name +
                    "_PeakIndexer.p",
                    "wb"))
            self.add_to_attr('peakIndexers', PI)
        self.set_attr('gvectorEvaluator', GE)
        self.save_geometries_as_yml()

    def evaluateGvectors(
            self,
            ds_tol=None,
            tth_gap=0.5,
            ds_gap=0.1,
            eta_gap=1,
            omega_gap=None,
            to_plot=True,
            save_arrays=False):
        print(single_separator + f'\nEVALUATING GVECTORS for {self.name}:')
        if ds_tol == 'auto':
            ds_tol = None
        self.gvectorEvaluator.remove_not_inrings(ds_tol=ds_tol)
        if not omega_gap:
            omega_gap = 2 * abs(self.sweepProcessors[0].sweep['omega_step'])
        self.gvectorEvaluator.calculate_ranges(
            tth_gap=tth_gap,
            ds_gap=ds_gap,
            eta_gap=eta_gap,
            omega_gap=omega_gap)
#             GE.group_gvectors(0.1, 1, 1)
        self.gvectorEvaluator.calc_histo(
            0.5, 0.5, plot=to_plot, save_arrays=save_arrays)
        self.gvectorEvaluator.save_gve(overwrite=True)
        pickle.dump(
            self.gvectorEvaluator,
            open(
                self.directory +
                self.name +
                "_GvectorEvaluator.p",
                "wb"))
        self.gvectorEvaluator.set_attr(
            'ds_eta_omega', np.zeros(
                (1, 1, 1)))  # To save memory

    def searchGrains(self, grainSpotter=None):
        print(single_separator + f'\nRUNNING GRAINSPOTTER for {self.name}')
        if grainSpotter:
            self.set_attr('grainSpotter', grainSpotter)
        self.grainSpotter.set_attr('spacegroup', self.material['spacegroup'])
        self.grainSpotter.set_attr('directory', self.directory)
        self.grainSpotter.set_attr('ini_file', self.name + '.ini')
        self.grainSpotter.set_attr('domega', abs(
            self.sweepProcessors[0].sweep['omega_step']))

        # Usually takes SOME time.
        self.grainSpotter.run_grainspotter(
            gve_file=self.name + '.gve',
            log_file=self.name + '.log')
        # Below we load .log file (more parameters but lower precision)
        # and .gff (fewer parameters but higher precision) and merge them
        grains_log = py3DXRD.Grain.load_log(self.directory, self.name+'.log')
        grains_gff = py3DXRD.Grain.load_gff(self.directory, self.name+'.gff')
        for gl in grains_log:
            gf = [g for g in grains_gff if g.grain_id == gl.grain_id][0]
            gl.set_attr('mean_IA', gf.mean_IA)
            gl.set_attr('pos_chisq', gf.pos_chisq)
            gl.set_attr('position', gf.position)
            gl.set_attr('u', gf.u)
            gl.set_attr('ubi', gf.ubi)
        self.set_attr('grains', grains_log)
        
        pickle.dump(
            self.grainSpotter,
            open(
                self.directory +
                self.name +
                "_GrainSpotter.p",
                "wb"))
        for g in self.grains:
            g.identify_measured(self.gvectorEvaluator.gvectors)

    def runPolyXSim(self, polyxsim):
        print(single_separator + f'\nRUNNING POLYXSIM for {self.name}')
        GS = self.grainSpotter
        PS = polyxsim
        sim_dir = self.directory + 'sim_' + self.material['name'] + '/'
        if not os.path.exists(sim_dir):
            os.makedirs(sim_dir)  # Create the directory if does not exist
        subprocess.call('rm -r ' + sim_dir, shell=True)
        PS.set_attr('directory', sim_dir)

        for g in self.grains:
            print(single_separator + f'\nGRAIN {g.grain_id}')
            if not g.size:
                g.set_attr('size', 0.05)
#             g.identify_measured(self.gvectorEvaluator.gvectors)
            GE_sim = py3DXRD.GvectorEvaluator(directory=self.directory)
            for SP in self.sweepProcessors:
                det_num = int(SP.name.split('_d')[-1])
                PS.set_attr('inp_file', f'g{g.grain_id:03d}_d{det_num}.inp')
                PS.set_attr('stem', f'g{g.grain_id:03d}_d{det_num}')
                PS.set_attr('geometry', SP.geometry)
                PS.set_attr('grains', [g])
                PS.save_inp(overwrite=True)
                d = min(
                    0.2,
                    np.sqrt(
                        g.position[0]**2 +
                        g.position[1]**2 +
                        g.position[2]**2))
                PS.set_attr('sample_cyl', [2 * d, 2 * d])
                PS.run_PolyXsim()  # Usually takes SOME time.
                GE_this_det = py3DXRD.GvectorEvaluator(directory=PS.direc)
                GE_this_det.set_attr('geometries', [PS.geometry])
                GE_this_det.load_gve(gve_file=PS.stem + '.gve')
                try:
                    GE_sim.absorb(GE_this_det)
                except BaseException:
                    GE_sim = GE_this_det    # if the very first sweep
            GE_sim.calc_histo(0.5, 0.5, plot=False, save_arrays=False)
            GE_sim.remove_not_ranges(
                GS.ds_ranges,
                GS.tth_ranges,
                GS.omega_ranges,
                convert_ranges(
                    GS.eta_ranges))
            # Sort in asceding ds
            ind = np.argsort([g['ds'] for g in GE_sim.gvectors])
            GE_sim.save_gve(gve_file=f'g{g.grain_id:03d}.gve', overwrite=True)
            g.set_attr('expected_gvectors', [GE_sim.gvectors[i] for i in ind])
            g.plot_measured_vs_expected()  # Usually takes LONG time.

    def mark_peaks(self):
        print(single_separator + f'\nMARKING GRAIN PEAKS for {self.name}')
        gvectors = []
        for G in self.grains:
            gvectors += G.measured_gvectors

        ind = np.argsort([g['omega'] for g in gvectors]
                         )  # Sort in asceding omega
        gvectors = [gvectors[i] for i in ind]

        peaks = []

        for g in gvectors:
            spot3d_id_reg = self.gvectorEvaluator.spot3d_id_reg
            GE_id = int(g['spot3d_id'] / self.gvectorEvaluator.spot3d_id_reg)
            spot3d_id = g['spot3d_id'] - \
                self.gvectorEvaluator.spot3d_id_reg * GE_id
            if GE_id == 0:
                GE = self.gvectorEvaluator
            else:
                GE = self.gvectorEvaluator.absorbed[GE_id - 1]
            PI = pickle.load(
                open(
                    GE.directory +
                    GE.name +
                    '_PeakIndexer.p',
                    "rb"))
            PI_id = int(spot3d_id / PI.spot3d_id_reg)
            if PI_id > 0:
                PI = PI.absorbed[PI_id - 1]
                spot3d_id = spot3d_id - int(spot3d_id / PI.spot3d_id_reg)
            p_list = [p for p in PI.peaks if p['spot3d_id'] == spot3d_id]

            if len(p_list) > 0:
                peaks.append(p_list[0])
            else:
                raise ValueError('peak for gvector ' +
                                 str(g['spot3d_id']) + ' not found!')

        print(f'Peaks to mark: {len(peaks)}')
        for SP in self.sweepProcessors:
            SP_FULL = pickle.load(
                open(
                    SP.directory +
                    SP.name +
                    "_SweepProcessor.p",
                    "rb"))
            omegas = [(omg[0] + omg[1]) / 2 for omg in SP_FULL.chunk['omegas']]
            omgeta_rngs = SP_FULL.projs['omgeta_rngs']
            immax_rngs = SP_FULL.projs['immax_rngs']
            ranges = SP_FULL.projs['ranges']

            for p in peaks:
                d_omg = abs(np.asarray(omegas) - p['omega'])
                omg_ind = np.where(d_omg == min(d_omg))[0][0]
#                 eta_ind = round(g['eta']-180)
                r = [
                    p['fc'] -
                    SP.geometry.y_center,
                    p['sc'] -
                    SP.geometry.z_center]
                thh_ind = round(np.sqrt(r[0]**2 + r[1]**2))
                eta_ind = round(
                    mod_360(
                        np.arctan2(
                            r[1],
                            r[0]) *
                        180 /
                        np.pi,
                        180))
                this_rng = [rng for rng in ranges if rng[0]
                            < thh_ind and rng[1] > thh_ind][0]
                rng_ind = ranges.index(this_rng)
                omgeta_rngs[rng_ind][omg_ind, eta_ind] = -1
                immax_rngs[round(p['fc']), round(p['sc'])] = -1

            path = SP.directory + 'projs_ranges/' + SP.name
            SP_FULL.add_to_log(
                'Writing file: ' +
                path +
                "_immax_rngs_marked.tif",
                True)
            SP.add_to_log(
                'Writing file: ' +
                path +
                "_immax_rngs_marked.tif",
                True)
            tifffile.imsave(path + "_immax_rngs_marked.tif", immax_rngs)
            for ir, rng in enumerate(ranges):
                SP_FULL.add_to_log(
                    'Writing file: ' +
                    path +
                    "_omgeta_rngs_marked.tif",
                    True)
                SP.add_to_log(
                    'Writing file: ' +
                    path +
                    "_omgeta_rngs_marked.tif",
                    True)
                tifffile.imsave(
                    path + f"_omgeta_rngs_marked-{rng[0]}-{rng[1]}.tif",
                    omgeta_rngs[ir])

            pickle.dump(
                SP_FULL,
                open(
                    SP.directory +
                    SP.name +
                    "_SweepProcessor.p",
                    "wb"))

        return

    
    def remove_duplicates(self, ang_tol, pos_tol):
        # Finds grains that are close to each other in positions and orientations
        # In such pairs duplicate is the one with worse metrics. This function removes such grains.
        print(double_separator+f'\nRemoving duplicates in {len(self.grains)} grains')
        print('Tolerances [ang, pos]: ', [ang_tol, pos_tol])
        if not ang_tol   >= 0: raise ValueError('ang_tol must be non-negative!')
        if not pos_tol   >= 0: raise ValueError('pos_tol must be non-negative!')
        self.set_attr('ang_tol', ang_tol)
        self.set_attr('pos_tol', pos_tol)
        duplicates = []
        
        for i1, g1 in enumerate(self.grains):
            q1 = len(g1.measured_gvectors)/g1.mean_IA # Measure of quality.
            for i2, g2 in enumerate(self.grains):
                if g2.grain_id <= g1.grain_id: continue # Skip those that have been already checked.
                d_pos = np.linalg.norm( np.asarray(g1.position) - np.asarray(g2.position) )
                if d_pos > pos_tol: continue # Skip if the positions are not close.
                d_ang = disorientation(g1.u, g2.u, sgno = self.material['spacegroup'])
                if d_ang > ang_tol: continue # Skip if the orinetations are not close.

                q2 = len(g2.measured_gvectors)/g2.mean_IA # Measure of quality.
                print("Grain #{:03d} and Grain #{:03d} differences:".format(g1.grain_id, g2.grain_id),
                      "Orientations = {0:.2f} deg.".format(d_ang),
                      "Positions = {0:.2f} mm.".format(d_pos),
                      "Quality = {0:.2f}".format(q1), "vs {0:.2f}".format(q2))
                duplicates.append( g1.grain_id if q1 < q2 else g2.grain_id ) # Duplicate is the one with worse metrics.
                
        duplicates = np.sort(list(set(duplicates)))
        print("List of grain_ids to remove:", duplicates)
        list_of_grains = [g for g in self.grains if g.grain_id not in duplicates]
        self.set_attr('grains', list_of_grains)
        return duplicates

    
    def absorb(self, x):
        print('Not implemented yet!')
        return


def plot_sinogram(name, list_DATApaths):
    import tifffile
    omegas = []
    ypos_s = []
    for path in list_DATApaths:
        print("PATH:", path)
        DATA = pickle.load(open(path, "rb"))
        SP = DATA.sweepProcessors[0]
        GE = DATA.gvectorEvaluator
        GS = DATA.grainSpotter
        GE.remove_not_ranges(
            GS.ds_ranges,
            GS.tth_ranges,
            GS.omega_ranges,
            GS.eta_ranges)
        # -1*0.5*SP.sweep['omega_step'] # useful for checking positions
        omg_sys_err = 0
        omegas = omegas + [g['omega'] + omg_sys_err for g in GE.gvectors]
        ypos_s = ypos_s + [DATA.position[1] for g in GE.gvectors]

#         for SP in DATA.sweepProcessors:
#             PI = SP.generate_PeakIndexer()
#             omg_sys_err = -0*0.5*SP.sweep['omega_step']
#             omegas = omegas + [p['omega'] + omg_sys_err for p in PI.peaks]
#             ypos_s = ypos_s + [DATA.position[1] for p in PI.peaks]
    fig = plt.figure(figsize=(30, 10))
    plt.scatter(omegas, ypos_s, s=1)

    omgy = np.zeros([len(SP.chunk['omegas']) + 2,
                    len(list_DATApaths) + 2], dtype=np.float32)
    omg_min = min(omegas)
    omg_stp = abs(SP.sweep['omega_step'])
    y_min = min(ypos_s)
    y_stp = (max(ypos_s) - min(ypos_s)) / (len(list_DATApaths) - 1)
    for omg, y in zip(omegas, ypos_s):
        omg_ind = (omg - omg_min) / omg_stp
        y_ind = (y - y_min) / y_stp
        omgy[round(omg_ind), round(y_ind)] += 1
    directory = DATA.directory.replace('/' + DATA.directory.split('/')[-2], '')
    name = DATA.name.replace('_'.join(DATA.name.split('_')[0:2]), name)
    tifffile.imsave(directory + name + '_sinogram.tif', omgy)


def set_MultiDATA(directory, name, list_DATApaths, num_0):
    def new_name(old_name):
        return old_name.replace('_'.join(old_name.split('_')[0:2]), name)

    path = list_DATApaths[num_0]
    print("PATH:", path)
    DATA = pickle.load(open(path, "rb"))
    DATA.directory = directory
    DATA.name = new_name(DATA.name)

    for ii, SP in enumerate(DATA.sweepProcessors):
        DATA.sweepProcessors[ii].directory = directory
        DATA.sweepProcessors[ii].name = new_name(DATA.sweepProcessors[ii].name)
        DATA.sweepProcessors[ii].geometry.directory = directory
        DATA.sweepProcessors[ii].geometry.par_file = new_name(
            DATA.sweepProcessors[ii].geometry.par_file)

    for ii, PI in enumerate(DATA.peakIndexers):
        DATA.peakIndexers[ii].directory = directory
        DATA.peakIndexers[ii].name = new_name(DATA.peakIndexers[ii].name)
        DATA.peakIndexers[ii].geometry.directory = directory
        DATA.peakIndexers[ii].geometry.par_file = new_name(
            DATA.peakIndexers[ii].geometry.par_file)

    DATA.gvectorEvaluator.directory = directory
    DATA.gvectorEvaluator.name = new_name(DATA.gvectorEvaluator.name)
    for ii, GM in enumerate(DATA.gvectorEvaluator.geometries):
        DATA.gvectorEvaluator.geometries[ii].directory = directory
        DATA.gvectorEvaluator.geometries[ii].par_file = new_name(
            DATA.gvectorEvaluator.geometries[ii].par_file)

    max_ds_tol = 0
    max_detshift = 0
    min_distance = 999999999999999999
    for path in list_DATApaths:
        if path == list_DATApaths[num_0]:
            continue
        DATA_ = pickle.load(open(path, "rb"))
        detshift = list(np.asarray(DATA.position) - np.asarray(DATA_.position))
        DATA_.index(move_det_xyz_mm=detshift)
        DATA.gvectorEvaluator.absorb(DATA_.gvectorEvaluator)
        max_ds_tol = max(max_ds_tol, DATA_.gvectorEvaluator.ds_tol)
        max_detshift = max(max_detshift, np.linalg.norm(detshift))
        min_distance = min(min_distance, min(
            [gm.distance for gm in DATA_.gvectorEvaluator.geometries]))

    print(max_ds_tol)
    # 1000 to convert detshift from mm to microns
    extra_tth = np.degrees(np.arctan(1000 * max_detshift / min_distance))
    extra_ds_tol = DATA.gvectorEvaluator.geometries[0].ds_from_tth(extra_tth)
    DATA.gvectorEvaluator.set_attr('ds_tol', max_ds_tol + extra_ds_tol)

    return DATA

#### TESTS    ###############
