import os, sys, subprocess, pickle
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
# sys.path.insert(0, '/home/shabalin/')
sys.path.insert(0, '/asap3/petra3/gpfs/common/p21.2/scripts/')
import py3DXRD
import experiment_settings
single_separator = "--------------------------------------------------------------"
double_separator = "=============================================================="


### SETTING THE SWEEPS AND OTHER OBJECTS:
def set_Geometry(det_num, material):
    GM = py3DXRD.Geometry() # initialization
    GM.load_par(directory = experiment_settings.path_gen + 'processed/calibration/',
                par_file = f'CeO2_avg_final.par') # load calibrated parameters
    GM.spline_file = None # Perfect detector or spline needed
    GM.set_attr('distance'   , 714498.555)
    GM.set_attr('spacegroup' , material['spacegroup'])
    GM.set_attr('symmetry'   , material['symmetry'  ])
    GM.set_attr('unitcell'   , material['unitcell'  ])
    return GM


def set_SweepProcessor(i_load, i_slow, i_fast, det_num):
    default_xyz = [0, 0, 0]
    y_points = np.linspace(-3.5, 3.5, 71) # y-motor positions if needed
    
#     meta_key = path_gen + f'raw/fio/eh3scan_00029.fio'
#     meta_key = path_gen + f'raw/newMgAl102/{i_load}/{i_load}.log'
#     meta_key = path_gen + f'raw/macrotest/{i_load}/{i_fast:03d}_y1_{y_points[i_fast]:.4f}.fio'
    meta_key = experiment_settings.path_gen + f'raw/polypd002/{i_load}/{i_fast:03d}_y1_{y_points[i_fast]:.4f}.fio' # metadata of the sweep
    
    SP = experiment_settings.set_p212_sweep(i_slow, i_fast, det_num, default_xyz, meta_key) # this function usually don't need modifications
    
    if SP.sweep['stem'][-1] != '_': SP.sweep['stem'] += '_' # correct for possible mismatch between file stems in the sweep command and actual file names
    SP.processing['options'] = None # ('flip', 'r270') # detector flips if needed to be applied before those in geometry file
#     SP.directory = SP.directory.replace('newMgAl102', 'newMgAl102_4') # if needed to modify the output directory (ie when the default one must not be overwritten)
    return SP


def set_GrainSpotter(material):
    GS = experiment_settings.set_grainspotter(material, domega=None)
#     GS.load_ini(ini_file = 'example.ini')

    GS.set_attr('tth_ranges'   , [ [8.0, 17.0] ] ) # 12.7]])
#     GS.set_attr('ds_ranges'    , [ [0.45, 1.2] ] ) # [0.5, 1.0]]) # GV.ds_ranges)
    GS.set_attr('eta_ranges'   , [ [5, 85], [275, 355]] ) # GV.eta_ranges)
    GS.set_attr('omega_ranges' , [ [-179.5,  179.5]] ) # GV.omega_ranges)
    GS.set_attr('cuts'         , [ 12, 0.6,  0.6] )
    GS.set_attr('uncertainties', [0.2, 1.5,  1.0] ) # [sigma_tth sigma_eta sigma_omega] in degrees
    GS.set_attr('nsigmas'      , 1)
    GS.set_attr('eulerstep'    , 6)
    GS.set_attr('Nhkls_in_indexing', None)
    GS.set_attr('random', 10000)
    GS.set_attr('positionfit', True)
    return GS


def set_PolySim(grainspotter = None, material = None):
#     PS = py3DXRD.PolySim(directory = path_gen+"")
#     PS.load_inp(inp_file = "")
    PS = experiment_settings.set_polyxsim(grainspotter, material)
    PS.set_attr('inp_file', grainspotter.log_file.strip('.log'))
    PS.set_attr('beamflux', 1e12)
    PS.set_attr('beampol_factor', 1)
    PS.set_attr('beampol_direct', 0)
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


def set_DATA(i_load, i_slow, i_fast, detectors, material):
    DATA = py3DXRD.DataAnalysis()
    DATA.set_attr('material', material)
    for det_num in detectors:
        SP = set_SweepProcessor(i_load, i_slow, i_fast, det_num)
        GM = set_Geometry(det_num, material)
        SP.set_attr('geometry', GM)
        DATA.add_to_attr('sweepProcessors', SP)
    DATA.set_attr('yml_det_order', [1])
    DATA.set_attr('directory'    , SP.directory)
    DATA.set_attr('name'         , f's{i_slow:03d}_f{i_fast:03d}_'+material['name'])
    DATA.set_attr('position'     , SP.position)
    DATA.set_attr('rotation'     , [0,0,0])
    if SP.log_meta:
        load_values = [ent['load'] for ent in SP.log_meta['entries']]
        DATA.set_attr('pressure', sum(load_values) / len(load_values))
    return DATA


### SELECTING WHICH SWEEPS TO ANALYZE:
# loads_str = "1293  1305  1317 1419  1431  1443  1455  1467"
# print(single_separator+'\nLOADS:', [int(v) for v in loads_str.split()].sort())
load_states       = [1] # Indices of     loads to analyze
slow_translations = [0]                   #,1,2] # Indices of positions to analyze
fast_translations = list(range(0,71))     #,1,2,4,5,6,7,8] # Indices of positions to analyze
detectors         = [3]                   # Detector code to analyze (4 for Varex-4 etc)
material = experiment_settings.materials_table['Pd']

### LOADING RAW DATA. PEAKSERCHING, MERGING, APPLYING INSTRUMENT CONFIGURATION.
for i_load in load_states[:]: # 
    for i_slow in slow_translations[:]:     # (e.g. idty1).
        for i_fast in fast_translations[0:1]: # (e.g. idtz2).
            print(double_separator + f'\nLOAD = {i_load}, TRANSLATIONS: slow = {i_slow}, fast = {i_fast}')
            
            DATA = set_DATA(i_load, i_slow, i_fast, detectors, material)
            
            DATA.process_images(frames = 'all', thr = 'auto')
            DATA.peaksearch(peaksearch_thresholds = 'auto', peakmerge_thresholds = 'auto', min_peak_dist = 10)
            pickle.dump(DATA, open(DATA.directory+DATA.name+"_DATA.p","wb") )
            
#             DATA = pickle.load(open(DATA.directory+DATA.name+"_DATA.p","rb") )
            DATA.index(move_det_xyz_mm = [0, -1*DATA.position[1], 0])
            DATA.evaluateGvectors(tth_gap=0.5, ds_gap=0.1, eta_gap=1)
            DATA.searchGrains(grainSpotter = set_GrainSpotter(material))
            pickle.dump(DATA, open(DATA.directory+DATA.name+"_DATA.p","wb") )
            
#             DATA = pickle.load(open(DATA.directory+DATA.name+"_DATA.p","rb") )
            DATA.runPolyXSim(polyxsim = set_PolySim(DATA.grainSpotter))
print('DONE!')

list_DATApaths = []
for i_load in load_states[:]: # 
    for i_slow in slow_translations[:]:     # (e.g. idty1).
        for i_fast in fast_translations[:]: # (e.g. idtz2).
            print(double_separator + f'\nLOAD = {i_load}, TRANSLATIONS: slow = {i_slow}, fast = {i_fast}')
            DATA = set_DATA(i_load, i_slow, i_fast, detectors, material)
            list_DATApaths.append(DATA.directory+DATA.name+"_DATA.p")

list_DATApaths = [p for n, p in enumerate(list_DATApaths) if p not in list_DATApaths[:n]]
print(list_DATApaths)

from py3DXRD.DataAnalysis import plot_sinogram, set_MultiDATA
plot_sinogram(list_DATApaths)
DATA_ALL = set_MultiDATA(list_DATApaths, 35)

DATA_ALL.evaluateGvectors(tth_gap=0.5, ds_gap=0.1, eta_gap=1)
DATA_ALL.searchGrains(grainSpotter = set_grainspotter())
pickle.dump(DATA_ALL, open(DATA_ALL.directory+DATA_ALL.name+"_DATA.p","wb") )
# DATA.runPolyXSim(polyxsim = set_polyxsim(DATA_ALL.grainSpotter, material = DATA_ALL.material))