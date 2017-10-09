#!/usr/bin/python3
# -*- coding: utf-8 -*-

import modules.triggerplotModule
import modules.dataModule
import numpy
import pandas

def JetCuts(data):
    return (data > 65)


# quantities settings
HT = {'key': 'HT', 'label': 'HT in GeV', 'limits': [500, 2000, 30]}
Mjj = {'key': 'Mjj', 'label': 'invariant dijetmass $M_{jj}$ in GeV', 'limits': [500, 2000, 30]}

pt1 = {'key': 'pt1', 'label': 'leading jet $p_t$', 'limits': [0, 1000, 20]}
pt2 = {'key': 'pt2', 'label': 'subleading jet $p_t$', 'limits': [0, 1000, 20]}
softdrop1 = {'key': 'softdrop1', 'label': 'leading jet softdropmass in GeV', 'limits': [0, 200, 5]}
softdrop2 = {'key': 'softdrop2', 'label': 'subleading jet softdropmass in GeV', 'limits': [0, 200, 5]}

trimMasstrigger = ['HLT_AK8PFHT750_TrimMass50_v', 'HLT_AK8PFHT800_TrimMass50_v', 'HLT_AK8PFHT850_TrimMass50_v', 'HLT_AK8PFHT900_TrimMass50_v']

# load data
data = modules.dataModule.loadData('Sync')
data = data[numpy.logical_and(abs(data['eta1']) < 2.4, abs(data['eta1']) < 2.4)] # etacut for synchronization (cut on data is 2.5)
data = data[JetCuts(data['softdrop1']) | JetCuts(data['softdrop2'])]

modules.dataModule.printRunlist(data)


# define tags
data['1tag'] = JetCuts(data['softdrop1']) | JetCuts(data['softdrop2']) # at least one jet passed jetcut
data['2tag'] = JetCuts(data['softdrop1']) & JetCuts(data['softdrop2']) # both jets passed jetcut



### define combinations

# denominator
data['PFHT430_OR_PFHT780_OR_PFJet500'] = data['HLT_PFHT430_v'] + data['HLT_PFHT780_v'] + data['HLT_PFJet500_v'] > 0
data['Mu50_OR_IsoMu27'] = data['HLT_Mu50_v'] + data['HLT_IsoMu27_v'] > 0

# PFHT1050_OR_AK8PFJet500
data['PFHT1050_OR_AK8PFJet500'] = data['HLT_PFHT1050_v'] + data['HLT_AK8PFJet500_v'] > 0

# OR_substructure
data['OR_substructure'] = data['HLT_AK8PFHT750_TrimMass50_v'] + data['HLT_AK8PFHT800_TrimMass50_v'] + data['HLT_AK8PFHT850_TrimMass50_v'] + data['HLT_AK8PFHT900_TrimMass50_v'] > 0

# all_trigger
data['OR_all_trigger'] = numpy.logical_or(data['OR_substructure'], data['PFHT1050_OR_AK8PFJet500'])


# define datasets
full = {'data': data, 'label': 'SingleMuon (no filter)', 'key': 'dataset', 'sets': ['SingleMuon-postfix'], 'denom': 'HLT_IsoMu27_v'}
data['dataset'] = data['dataset'].mask(data['run'].values > 299504, other='SingleMuon-postfix')

menu = {'data': data, 'label': 'SingleMuon (no filter)', 'key': 'menu', 'sets': ['menu_2017_v1', 'menu_2017_v2', 'menu_2017_v3'], 'denom': 'HLT_IsoMu27_v'}
data['menu'] = data['dataset']
data['menu'] = data['menu'].mask(numpy.logical_and(data['run'].values <= 299329, data['run'].values >= 296070), other='menu_2017_v1')
data['menu'] = data['menu'].mask(numpy.logical_and(data['run'].values <= 302019, data['run'].values >= 299368), other='menu_2017_v2')
data['menu'] = data['menu'].mask(numpy.logical_and(data['run'].values <= 302479, data['run'].values >= 302026), other='menu_2017_v3')


runs = {'data': data, 'label': 'SingleMuon (no filter)', 'key': 'runs', 'sets': ['RunB', 'RunC', 'RunD'], 'denom': 'HLT_IsoMu27_v'}
data['runs'] = data['dataset']
data['runs'] = data['runs'].mask(numpy.logical_and(data['run'].values <= 299329, data['run'].values >= 297046), other='RunB')
data['runs'] = data['runs'].mask(numpy.logical_and(data['run'].values <= 302029, data['run'].values >= 299368), other='RunC')
data['runs'] = data['runs'].mask(numpy.logical_and(data['run'].values <= 303434, data['run'].values >= 302030), other='RunD')




# define file
texfile = 'Sync.tex'
modules.triggerplotModule.clearfile(texfile)
modules.triggerplotModule.write2tex('\section{Synchronization}\n', texpath=texfile)


modules.triggerplotModule.doEffPlot(full, trigger=['HLT_PFHT1050_v', 'HLT_AK8PFJet500_v', 'OR_substructure', 'OR_all_trigger'], quant=HT, texpath=texfile, fit=False)
modules.triggerplotModule.doEffPlot(full, trigger=['HLT_PFHT1050_v', 'HLT_AK8PFJet500_v', 'OR_substructure', 'OR_all_trigger'], quant=Mjj, texpath=texfile, fit=False)

#modules.triggerplotModule.doEffPlot(runs, trigger=['OR_all_trigger'], quant=Mjj, texpath=texfile, fit=False)
#modules.triggerplotModule.doEffPlot(menu, trigger=['OR_all_trigger'], quant=Mjj, texpath=texfile, fit=False)

#modules.triggerplotModule.doEffPlot(runs, trigger=['HLT_PFHT1050_v'], quant=Mjj, texpath=texfile, fit=False)
#modules.triggerplotModule.doEffPlot(menu, trigger=['HLT_PFHT1050_v'], quant=Mjj, texpath=texfile, fit=False)

modules.triggerplotModule.do2DPlot(full, 'PFHT1050_OR_AK8PFJet500', Mjj, softdrop1, texpath=texfile, cuts=[[-1, 10000],[65, 105]])
modules.triggerplotModule.do2DPlot(full, 'OR_substructure', Mjj, softdrop1, texpath=texfile, cuts=[[-1, 10000],[65, 105]])
modules.triggerplotModule.do2DPlot(full, 'OR_all_trigger', Mjj, softdrop1, texpath=texfile, cuts=[[-1, 10000],[65, 105]])