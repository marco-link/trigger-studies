#!/usr/bin/python3
# -*- coding: utf-8 -*-

import modules.triggerplotModule
import modules.dataModule
import numpy
import pandas

def JetCuts(data):
    return (data > 65) & (data < 105)


# quantities settings
HT = {'key': 'HT', 'label': 'HT in GeV', 'limits': [500, 2000, 30]}
Mjj = {'key': 'Mjj', 'label': 'invariant dijetmass $M_{jj}$ in GeV', 'limits': [500, 2000, 30]}

pt1 = {'key': 'pt1', 'label': 'leading jet $p_t$', 'limits': [0, 1000, 20]}
pt2 = {'key': 'pt2', 'label': 'subleading jet $p_t$', 'limits': [0, 1000, 20]}
softdrop1 = {'key': 'softdrop1', 'label': 'leading jet softdropmass in GeV', 'limits': [0, 200, 5]}
softdrop2 = {'key': 'softdrop2', 'label': 'subleading jet softdropmass in GeV', 'limits': [0, 200, 5]}

eta1 = {'key': 'eta1', 'label': 'leading jet $\\eta$', 'limits': [-2.5, 2.5, 0.1]}
eta2 = {'key': 'eta2', 'label': 'subleading jet $\\eta$', 'limits': [-2.5, 2.5, 0.1]}

phi1 = {'key': 'phi1', 'label': 'leading jet $\\varphi$', 'limits': [-numpy.pi, numpy.pi, 0.1]}
phi2 = {'key': 'phi2', 'label': 'subleading jet $\\varphi$', 'limits': [-numpy.pi, numpy.pi, 0.1]}

trimMasstrigger = ['HLT_AK8PFHT750_TrimMass50_v', 'HLT_AK8PFHT800_TrimMass50_v', 'HLT_AK8PFHT850_TrimMass50_v', 'HLT_AK8PFHT900_TrimMass50_v']

# load data
data = modules.dataModule.loadData('SingleMuon')

modules.dataModule.printRunlist(data)


# define tags
data['1tag'] = JetCuts(data['softdrop1']) | JetCuts(data['softdrop2']) # at least one jet passed jetcut
data['2tag'] = JetCuts(data['softdrop1']) & JetCuts(data['softdrop2']) # both jets passed jetcut

data['Mjj>1200GeV'] = data['Mjj'] > 1200


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
full = {'data': data, 'label': 'SingleMuon', 'key': 'dataset', 'sets': ['SingleMuon-postfix'], 'denom': 'Mu50_OR_IsoMu27'}
data['dataset'] = data['dataset'].mask(data['run'].values > 299504, other='SingleMuon-postfix')

menu = {'data': data, 'label': 'SingleMuon', 'key': 'menu', 'sets': ['menu_2017_v1', 'menu_2017_v2', 'menu_2017_v3'], 'denom': 'Mu50_OR_IsoMu27'}
data['menu'] = data['dataset']
data['menu'] = data['menu'].mask(numpy.logical_and(data['run'].values <= 299329, data['run'].values >= 296070), other='menu_2017_v1')
data['menu'] = data['menu'].mask(numpy.logical_and(data['run'].values <= 302019, data['run'].values >= 299368), other='menu_2017_v2')
data['menu'] = data['menu'].mask(numpy.logical_and(data['run'].values <= 302479, data['run'].values >= 302026), other='menu_2017_v3')


runs = {'data': data, 'label': 'SingleMuon', 'key': 'runs', 'sets': ['RunB', 'RunC', 'RunD', 'RunE', 'RunF'], 'denom': 'Mu50_OR_IsoMu27'}
data['runs'] = data['dataset']
data['runs'] = data['runs'].mask(numpy.logical_and(data['run'].values <= 299329, data['run'].values >= 297046), other='RunB')
data['runs'] = data['runs'].mask(numpy.logical_and(data['run'].values <= 302029, data['run'].values >= 299368), other='RunC')
data['runs'] = data['runs'].mask(numpy.logical_and(data['run'].values <= 302663, data['run'].values >= 302031), other='RunD')
data['runs'] = data['runs'].mask(numpy.logical_and(data['run'].values <= 304797, data['run'].values >= 303572), other='RunE')
data['runs'] = data['runs'].mask(numpy.logical_and(data['run'].values <= 306138, data['run'].values >= 305040), other='RunF')



# define file
texfile = 'SingleMuon.tex'
modules.triggerplotModule.clearfile(texfile)
modules.triggerplotModule.write2tex('\section{SingleMuon}\n', texpath=texfile)


modules.triggerplotModule.doEffPlot(full, trigger=['HLT_PFHT1050_v', 'PFHT1050_OR_AK8PFJet500', 'OR_substructure', 'OR_all_trigger'], quant=Mjj, texpath=texfile, fit=False)

modules.triggerplotModule.doEffPlot(full, trigger=['HLT_PFHT1050_v', 'PFHT1050_OR_AK8PFJet500', 'OR_substructure', 'OR_all_trigger'], quant=softdrop1, texpath=texfile, fit=False, mask='Mjj>1200GeV')

modules.triggerplotModule.doEffPlot(runs, trigger=['OR_all_trigger'], quant=Mjj, texpath=texfile, fit=False)
modules.triggerplotModule.doEffPlot(menu, trigger=['OR_all_trigger'], quant=Mjj, texpath=texfile, fit=False)

# define plots
for tag in ['1tag', '2tag']:
    modules.triggerplotModule.doEffPlot(full, trigger=['HLT_PFHT1050_v', 'PFHT1050_OR_AK8PFJet500', 'OR_substructure', 'OR_all_trigger'], quant=Mjj, texpath=texfile, fit=False, mask=tag)


modules.triggerplotModule.do2DPlot(full, 'OR_substructure', Mjj, softdrop1, texpath=texfile, cuts=[[-1, 10000],[65, 105]])


texfile = 'SingleMuon_appendix.tex'
modules.triggerplotModule.clearfile(texfile)



for tag in ['1tag', '2tag']:
    modules.triggerplotModule.doEffPlot(runs, trigger=['OR_all_trigger'], quant=Mjj, texpath=texfile, fit=False, mask=tag)
    modules.triggerplotModule.doEffPlot(menu, trigger=['OR_all_trigger'], quant=Mjj, texpath=texfile, fit=False, mask=tag)

    modules.triggerplotModule.doEffPlot(menu, trigger=['HLT_PFHT1050_v'], quant=Mjj, texpath=texfile, fit=False, mask=tag)
    modules.triggerplotModule.doEffPlot(menu, trigger=['HLT_AK8PFJet500_v'], quant=Mjj, texpath=texfile, fit=False, mask=tag)

    for trig in trimMasstrigger:
            modules.triggerplotModule.doEffPlot(menu, trigger=[trig], quant=Mjj, texpath=texfile, fit=False, mask=tag)

modules.triggerplotModule.do2DPlot(full, 'PFHT1050_OR_AK8PFJet500', Mjj, softdrop1, texpath=texfile, cuts=[[-1, 10000],[65, 105]])
modules.triggerplotModule.do2DPlot(full, 'OR_all_trigger', Mjj, softdrop1, texpath=texfile, cuts=[[-1, 10000],[65, 105]])