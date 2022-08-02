def load_p212_log(path):
    p212_log = {'file': path.split('/')[-1]}
    p212_log['directory'] = path.replace(p212_log['file'], '').replace('../', '')
    with open(path, "r") as f:
        p212_log['entries'] = []
        p212_log['data_and_time'] = f.readline()[:-1]
        for line in f:
            if line[0] == '#':
                if 'crosshead' in line and 'position' in line:
                    p212_log['crosshead_position'] = float(line.split()[2])
                elif 'sweep' in line:
                    cmd = line[1:-1]
                    for s in ['\'', ',', ']','[', '--']: cmd = cmd.replace(s,'')
                    p212_log['command'] = cmd
                elif 'Detector' in line:
                    p212_log['det_num'] = (line.split()[1])
                elif 'timestamp' in line:
                    titles = line[1:].split()
            else:
                words = line.split()
                if len(words) == len(titles):
                    v = [float(v) for v in words]
                    p212_log['entries'].append( dict(zip(titles,v)) )
    f.close()
    return p212_log