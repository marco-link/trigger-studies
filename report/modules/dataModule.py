#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import numpy
import pandas

def loadData(dataset, limits=None):
    merge(dataset)
    return getData(dataset, limits)


def getData(dataset, limits=None):
    data = []

    infiles = os.listdir('data/{}'.format(dataset))
    infiles = [infiles] if not isinstance(infiles,  list) else infiles # handling singles
    if not infiles == []:
        if not (limits==None):
            x = []
            for f in infiles:
                i = int(f.replace('.csv', ''))
                if(i>=limits[0] and i<=limits[1]):
                    x.append(f)
            infiles = x

        for f in infiles:
            print('loading file: {}/{}'.format(dataset, f))
            temp = pandas.read_csv('data/{}/{}'.format(dataset, f), sep=',', header=0, comment='#')
            temp['dataset'] = dataset
            data.append(temp)

    if(data==[]):
        print('NO FILES FOUND')
        return pandas.DataFrame(columns=['ID', 'run'])
    else:
        data = pandas.concat(data)
        for k in data.keys():
            if(k in ['ID', 'run', 'lumi'] or 'HLT_' in k):
                data[k] = data[k].astype(int)
        return data


def merge(dataset):
    infiles = os.listdir('IN')

    if(infiles == []):
        print('nothing to merge!')
    else:
        print('merging new files: {}'.format(infiles))
        events = 0
        confls = 0

        infiles = [infiles] if not isinstance(infiles,  list) else infiles # handling singles
        for f in infiles:
            if dataset in f:
                newdata = pandas.read_csv('IN/{}'.format(f), sep=',', header=0, comment='#')
                header = ','.join(newdata.keys().values)
                out = []
                for k in newdata.keys():
                    if(k in ['ID', 'run', 'lumi'] or 'HLT_' in k):
                        out.append('{:.0f}')
                    else:
                        out.append('{:.3f}')
                mask = ', '.join(out) + '\n'

                data = getData(dataset, limits=[numpy.amin(newdata['run']), numpy.amax(newdata['run'])])
                move = True

                print('\n\nmerging {:d} events from {}'.format(len(newdata), f))

                for idx, event in newdata.iterrows():
                    events = events + 1
                    if not (event['ID'] in data['ID'][data['run']==event['run']].values):
                        save(dataset, event, header, mask)
                    else:
                        print('conflict merging run: {:d} event: {:d} from {}'.format(int(event['run']), int(event['ID']), f))
                        confls = confls + 1
                        move = False

                    if(idx%1000 == 1):
                        print('[{:4.1f}%] processing {:d}/{:d} from {} ({:d} total conflicts)'.format(100*idx/len(newdata), idx, len(newdata), f, confls))

                if move:
                    print('no conflicts! moving {} to data/raw'.format(f))
                    os.rename('IN/{}'.format(f), 'data/raw/{}'.format(f))

        print('merged {} events, found {} conflicts'.format(events, confls))


def save(dtset, event, header, mask):
    path = 'data/{}/{:d}.csv'.format(dtset, int(event['run']))

    if not os.path.isfile(path):
        f = open(path, 'w')
        f.write('{}\n'.format(header))
    else:
        f = open(path, 'a')

    f.write(mask.format(*event.values))
    f.close()


def getRunlist(data):
    runlist = pandas.Series(data['run'].values.ravel()).unique()
    return numpy.sort(runlist)


def printRunlist(data):
    print('runlist:\n', getRunlist(data))


def getLumi(data, path='data/lumi.csv'):
    runlist = getRunlist(data)
    lumilist = pandas.read_csv(path, sep=',', header=0, comment='#')
    lumilist['run'], lumilist['fill'] = lumilist['run:fill'].str.split(':', 1).str
    lumi = 0

    for run in runlist:
        l = lumilist['recorded(/ub)'][lumilist['run'].astype(int) == run].values

        if len(l)==1:
            lumi = lumi + l[0]/1e9
        else:
            print('ERROR: run {} not found in lumilist, luminosity calculation aborted'.format(run))
            return -1

    print('luminosity in 1/fb: {:.3f}'.format(lumi))
    return lumi