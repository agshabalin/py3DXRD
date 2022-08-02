import os, sys
import numpy as np
sys.path.insert(0, '/home/shabalin/py3DXRD/')
from Geometry import Geometry
from PolySim  import PolySim

if __name__=="__main__":
    if len(sys.argv) != 3:
        raise ValueError('Provide 2 arguments: input_geometry_file, output_geometry_file!')

    input_geometry_file  = sys.argv[1]
    output_geometry_file = sys.argv[2]

    if   input_geometry_file[0] == '/':
        pass
    elif input_geometry_file[0:2] == './':
        input_geometry_file  = os.getcwd() +  input_geometry_file[1:]
        output_geometry_file = os.getcwd() + output_geometry_file[1:]
    else:
        input_geometry_file  = os.getcwd() + '/' + input_geometry_file
        output_geometry_file = os.getcwd() + '/' + output_geometry_file    



    inp_filename = input_geometry_file.split('/')[-1]
    inp_directory = input_geometry_file.replace(inp_filename, '')
    out_filename = output_geometry_file.split('/')[-1]
    out_directory = output_geometry_file.replace(out_filename, '')
    print('Input  directory:',inp_directory)
    print('Input  filename:', inp_filename)
    print('Output  directory:',out_directory)
    print('Output  filename:', out_filename)
    print('Output geometry file:', output_geometry_file)

    PS = PolySim()
    PS.geometry = Geometry()

    if   inp_filename.split('.')[1] == 'inp':
        PS.load_inp(directory = inp_directory, inp_file = inp_filename)
    elif inp_filename.split('.')[1] == 'par':
        PS.geometry.load_par( directory = inp_directory, par_file  = inp_filename)
    elif inp_filename.split('.')[1] == 'yml':
        PS.geometry.load_yml( directory = inp_directory, yml_file  = inp_filename)
    elif inp_filename.split('.')[1] == 'poni':
        PS.geometry.load_poni(directory = inp_directory, poni_file = inp_filename)

    if   out_filename.split('.')[1] == 'inp':
        PS.save_inp(directory = out_directory, inp_file = out_filename)
    elif out_filename.split('.')[1] == 'par':
        PS.geometry.save_par( directory = out_directory, par_file  = out_filename)
    elif out_filename.split('.')[1] == 'yml':
        PS.geometry.save_yml( directory = out_directory, yml_file  = out_filename)
    elif out_filename.split('.')[1] == 'poni':
        PS.geometry.save_poni(directory = out_directory, poni_file = out_filename)