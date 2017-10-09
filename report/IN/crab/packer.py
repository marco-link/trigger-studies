#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import numpy
import pandas

dataset = 'Sync(nofilter)'

for path in os.listdir('.'):
    if not '.' in path:
        data = []

        for f in os.listdir(path):
            if('.root' in f):
                print(f)
                data.append(pandas.read_csv('{}/{}'.format(path, f), sep=',', header=0, comment='#'))

        data = pandas.concat(data)

        header = ','.join(data.keys().values)
        print(header)

        out = []
        for k in data.keys():
            if(k in ['ID', 'run', 'lumi'] or 'HLT_' in k):
                out.append('{:.0f}')
            else:
                out.append('{:.3f}')
        mask = ', '.join(out) + '\n'


        f = open('{}-{}.csv'.format(dataset, path), 'w')
        f.write(header + '\n')
        for x in data.values:
            f.write(mask.format(*x))
