single_separator = "--------------------------------------------------------------\n"
double_separator = "==============================================================\n"


class peak:
    
    
    def __init__(self, titles, values, sourcefile):
        self.source = sourcefile
        self.titles = titles
        self.values = [float(v) if '.' in v else int(v) for v in values]
        self.dict = dict( zip(self.titles, self.values) )
        if 'spot3d_id' in self.dict.keys():
            self.spot3d_id = self.dict['spot3d_id']
        else:
            self.spot3d_id = None
        
    def show(self):
        print(single_separator + f'Spot3d_id:{self.spot3d_id}')
        print('Source file:', self.source)
        print(self.dict, '\n')
        

    def read_flt_file(filepath):
        if filepath[-4:] != '.flt': raise ValueError(filepath + ' - is not a .flt file!')
        print(double_separator + 'Reading file:' + filepath)
        list_of_peaks = []
        titles = "sc  fc  omega  Number_of_pixels  avg_intensity  s_raw  f_raw  sigs  sigf  covsf  sigo  covso  covfo  sum_intensity  sum_intensity^2  IMax_int  IMax_s  IMax_f  IMax_o  Min_s  Max_s  Min_f  Max_f  Min_o  Max_o  dety  detz  onfirst  onlast  spot3d_id".split()
        print(f'Titles (default):\n#  '+'  '.join(titles))
        f = open(filepath,"r")
        line = f.readline()
        while line:
            if '#  sc  fc  omega' in line:
                titles = line[2:-1].split()
                print(f'Titles (updated) :\n#  '+'  '.join(titles))
            else:
                words = line.split()
                if len(words) == len(titles):
                    list_of_peaks.append(peak(titles, words, filepath))
            line = f.readline()
        f.close()
        print(f'{len(list_of_peaks)} peaks loaded.\n')
        return list_of_peaks
    
    
    def write_flt_file(list_of_peaks, filepath):
        print(double_separator + 'Writing file: ' + filepath)
        titles = list(list_of_peaks[0].dict.keys())
        print(f'Titles:\n#  '+'  '.join(titles[0:30]))
        f = open(filepath,"w")
        f.write('#  ' + '  '.join(titles[0:30]) + '\n' )
        for p in list_of_peaks:
            values =[f'{v:d}' if type(v)==type(1) else f'{v:.4f}' for v in list(p.dict.values())]
            f.write('  ' + '  '.join([v for v in values[0:30]]) + '\n' )
        f.close()            
        print(f'{len(list_of_peaks)} peaks saved.\n')
        return