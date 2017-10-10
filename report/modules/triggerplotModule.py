#!/#usr/bin/python3
# -*- coding: utf-8 -*-

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import scipy.optimize
import scipy.interpolate
import scipy.stats
import scipy.special

import modules.dataModule

def clearfile(path):
    # delete file at beginning
    f = open(path, 'w')
    f.close()


def makeTable(data, labels, section, texpath):
    f = open(texpath, 'a')
    f.write('\\section*{{{}}}\n'.format(section))
    f.write('\\begin{frame}{\insertsection}\n')
    f.write('\\begin{center}\n')
    f.write('\\begin{tabular}{lc}\n')
    f.write('\\toprule\n')
    f.write('trigger & passed\\\\ \n')
    f.write('\midrule\n')
    for l in labels:
        f.write('{} & {} \\\\\n'.format(l.replace('_v', '').replace('_', '\_'), numpy.sum(data[l])))
    f.write('\\bottomrule\n')
    f.write('\\end{tabular}')
    f.write('\\end{center}\n')
    f.write('\\end{frame}\n\n')
    f.close()


def write2tex(txt, texpath):
    f = open(texpath, 'a')
    f.write('{}\n'.format(txt))
    f.close()


def makeSlide(name, texpath, caption=''):
    f = open(texpath, 'a')
    f.write('\\begin{frame}\n')
    f.write('\\begin{figure}\n')
    f.write('\\begin{center}\n')
    f.write('\includegraphics[width=\linewidth]{{{}}}\n'.format(name))
    f.write('\\end{center}\n')
    if not caption == '':
        f.write('\caption{{{}}}\n'.format(caption.replace('_', '\_')))
    f.write('\\end{figure}\n')
    f.write('\\end{frame}\n\n')
    f.close()


def getError(k, n, gamma=0.682):
    #clopper pearson
    # https://de.wikipedia.org/wiki/Konfidenzintervall_f%C3%BCr_die_Erfolgswahrscheinlichkeit_der_Binomialverteilung
    alpha = 0.5 * (1-gamma)
    lower = k/n - scipy.stats.beta.ppf(alpha, k, n-k+1)
    upper = scipy.stats.beta.ppf(1-alpha, k+1, n-k) - k/n

    lower[k==0] = 0
    upper[k==n] = 0

    return [lower, upper]


def getEfficiency(data, trigger, quant, denominator):
    # define bins
    width = quant['limits'][2]
    center = numpy.arange(quant['limits'][0], quant['limits'][1], width)
    bins = numpy.append(center - 0.5 * width, center[-1] + 0.5 * width)

    # apply denominator
    data = data[data[denominator]==1]

    denominatorentries = numpy.histogram(data[quant['key']], bins)[0]
    entries = numpy.histogram(data[quant['key']][data[trigger] == 1], bins)[0]

    eff = entries/denominatorentries
    mask = numpy.isfinite(eff)

    sigma = getError(entries[mask], denominatorentries[mask])

    return center[mask], eff[mask], sigma


def doEffPlot(dataset, trigger, quant, texpath, fit=False, x0=[0.9, 0.005, 1000], mask=''):
    data = dataset['data']
    if not mask=='':
        data = data[data[mask] == 1]

    fig = matplotlib.pyplot.figure(figsize=(10, 6))
    p1 = fig.add_subplot(111)
    for trig in trigger:
        for dtset in dataset['sets']:
            x, y, sigma = getEfficiency(data[data[dataset['key']] == dtset], trig, quant, dataset['denom'])
            lumi = modules.dataModule.getLumi(data[data[dataset['key']] == dtset])

            if not(numpy.sum(y)==0):
                col = p1.errorbar(x, y, yerr=sigma, fmt='.', label='{}\n{} ({:.2f} fb$^{{-1}}$)'.format(trig, dtset, lumi), capsize=2)[0].get_color()
                if(fit):
                    xfit, yfit, label, success, p, perr, c = doFit(x, y, sigma, x0=x0)
                    if (success):
                        p1.plot(xfit, yfit, ':', color=col, label=label)

    p1.text(0.01, 0.01, 'denominator: {}'.format(dataset['denom'].replace('_', ' ')), verticalalignment='bottom', horizontalalignment='left', transform=p1.transAxes)
    if (mask==''):
        label = ', '.join(trigger).replace('_v', '')
    else:
        label = '{} ({})'.format(', '.join(trigger).replace('_v', ''), mask)
    p1.set_title(label)
    p1.set_xlabel(quant['label'])
    p1.set_ylabel('efficiency')

    p1.set_xlim(quant['limits'][0], quant['limits'][1])
    p1.set_ylim(-0.05, 1.05)

    p1.set_position([0.1, 0.1, 0.6, 0.8])
    p1.legend(loc='upper left', bbox_to_anchor = (1, 1))
    p1.grid(alpha=0.5)

    #fig.tight_layout()
    out = 'res/{}_{}-{}-{}-{}-{}.pdf'.format(mask, dataset['label'], ','.join(trigger).replace('_v', ''), dataset['denom'], dataset['key'], quant['key']).replace(' ', '_')
    makeSlide(out, texpath, caption='')#'{} on {}'.format(', '.join(trigger).replace('_v', ''), ', '.join(dataset['sets'])))
    fig.savefig(out, dpi=300, transparent=False)
    #matplotlib.pyplot.show()
    matplotlib.pyplot.close()


def doFit(xdata, ydata, sigma, x0, cuteff=0.99):
    # only fit on efficiency > 0.6 on last entries connected
    mask = []
    y = True
    for x in numpy.flip(ydata > 0.6, 0):
        if not x:
            y = False
        mask.append(y)
    mask = numpy.flip(mask, 0)

    xdata = xdata[mask]
    ydata = ydata[mask]

    def getChi2(para, xdata, ydata, error, func):
        low, up = error
        fit = func(xdata, para)
        delta = ydata - fit

        sigma = []
        for i in range(len(fit)):
            if(delta[i]>0):
                sigma.append(low[i])
            else:
                sigma.append(up[i])
        x = (delta/sigma)**2

        return numpy.sum(x)

    def fitfunc(x, para):
        a, b, c = para
        return 0.5 * a * ( 1 + scipy.special.erf(b * (x - c)))

    res = scipy.optimize.minimize(lambda x: getChi2(x, xdata, ydata, sigma, fitfunc), x0=x0, tol=1e-4, method='Nelder-Mead')


    print(res)

    para = res.x
    Chi2 = res.fun
    dof = len(xdata)-len(para)
    pvalue = scipy.stats.chi2.sf(Chi2, dof)

    a, b, c = para
    cut = c - numpy.log(1/cuteff - 1)/b
    cut = scipy.special.erfinv(2*cuteff - 1) / b + c

    # estimate error on parameters
    def getParaError(func, value, esterr=1):
        res = scipy.optimize.minimize(lambda x: (func(x)-func(value)-1)**2, x0=[value+esterr], tol=1e-6, method='Nelder-Mead')
        up = res.x[0]

        res = scipy.optimize.minimize(lambda x: (func(x)-func(value)-1)**2, x0=[value-esterr], tol=1e-6, method='Nelder-Mead')
        down = res.x[0]

        return [abs(up-value), abs(down-value)]

    paraerr = [getParaError(lambda x: getChi2([x, para[1], para[2]], xdata, ydata, sigma, fitfunc), para[0], esterr=0.01), getParaError(lambda x: getChi2([para[0], x, para[2]], xdata, ydata, sigma, fitfunc), para[1], esterr=0.01), getParaError(lambda x: getChi2([para[0], para[1], x], xdata, ydata, sigma, fitfunc), para[2], esterr=100)]

    label = '$y = \\frac{{a}}{{2}} [1 + erf(b (x - c))]$\n$a={:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$\n$b={:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$\n$c={:.0f}^{{+{:.1f}}}_{{-{:.1f}}}$\n$\chi^2 / dof = {:.3f}/{}$\npvalue: {:.2f}\n{:.1f}%-plateau: {:.1f}'.format(para[0], *paraerr[0], para[1], *paraerr[1], para[2], *paraerr[2], Chi2, dof, pvalue, 100*cuteff, cut)
    if(len(xdata)>0):
        xfit = numpy.linspace(numpy.amin(xdata), numpy.amax(xdata), 1000)
    else:
        xfit = -999

    return xfit, fitfunc(xfit, para), label, res.success, para, paraerr, cut


def do2DPlot(dataset,  trigger, quant1, quant2, texpath, cuts=None, mask=''):
    data = dataset['data']
    if not mask=='':
        data = data[data[mask] == 1]
    datasets = dataset['sets']
    denominator = dataset['denom']

    # define bins
    width1 = quant1['limits'][2]
    center1 = numpy.arange(quant1['limits'][0], quant1['limits'][1] + width1, width1)
    bins1 = numpy.append(center1 - 0.5 * width1, center1[-1] + 0.5 * width1)

    width2 = quant2['limits'][2]
    center2 = numpy.arange(quant2['limits'][0], quant2['limits'][1] + width2, width2)
    bins2 = numpy.append(center2 - 0.5 * width2, center2[-1] + 0.5 * width2)

    # apply denominator & filter datasets
    data = data[data[denominator] == 1]
    data = data[data[dataset['key']].isin(datasets).values]
    lumi = modules.dataModule.getLumi(data[data[dataset['key']].isin(datasets).values])
    denominatorentries = numpy.histogram2d(data[quant1['key']], data[quant2['key']], bins=[bins1, bins2])[0]
    entries = numpy.histogram2d(data[quant1['key']][data[trigger]==1], data[quant2['key']][data[trigger]==1], bins=[bins1, bins2])[0]

    # plotting
    fig = matplotlib.pyplot.figure(figsize=(10,6))

    p1 = fig.add_subplot(111)

    eff = entries/denominatorentries
    eff[denominatorentries < 3] = numpy.nan
    eff = numpy.ma.masked_invalid(eff)

    cmap = matplotlib.pyplot.get_cmap('viridis')
    cmap.set_bad(color='w', alpha=0.8)

    fig.colorbar(p1.pcolormesh(bins1, bins2, eff.T, vmin=0, vmax=1, cmap=cmap))

    if mask == '':
        p1.set_title('{}'.format(trigger.replace('_v', '')))
    else:
        p1.set_title('{} ({})'.format(trigger.replace('_v', ''), mask))
    p1.set_xlabel(quant1['label'])
    p1.set_ylabel(quant2['label'])
    p1.text(0.01, 0.99, 'denominator: {}\ndatasets: {} ({:.2f} fb$^{{-1}}$)'.format(denominator.replace('_', ' '), ', '.join(datasets), lumi), verticalalignment='top', horizontalalignment='left', transform=p1.transAxes, color='r')

    p1.set_xlim(quant1['limits'][0], quant1['limits'][1])
    p1.set_ylim(quant2['limits'][0], quant2['limits'][1])

    if not (cuts == None):
        x = [cuts[0][0], cuts[0][0], cuts[0][1], cuts[0][1], cuts[0][0]]
        y = [cuts[1][0], cuts[1][1], cuts[1][1], cuts[1][0], cuts[1][0]]
        p1.plot(x, y, 'r-')


    fig.tight_layout()
    out = 'res/2D-{}-{}-{}_vs_{}_{}-{}.pdf'.format(dataset['label'], trigger, quant1['key'], quant2['key'], datasets, mask).replace(' ', '_')
    makeSlide(out, texpath, caption='2D: {} on {}'.format(trigger.replace('_v', ''), ', '.join(datasets)))
    fig.savefig(out, dpi=300, transparent=False)
    matplotlib.pyplot.close()


def makeRunPlot(data, dataset, trigger, quant, denominator, cuts, texpath):
    data = data[data[denominator]]

    fig = matplotlib.pyplot.figure(figsize=(10, 6))
    p1 = fig.add_subplot(111)
    x = []
    for run in data['run']:
        if not run in x:
            x.append(run)
    runs = numpy.sort(x)
    x = numpy.arange(len(runs))

    for cut in cuts:
        m = []
        n = []
        for r in runs:
            n.append(numpy.sum(data[trigger][numpy.logical_and(data[quant['key']] > cut, data['run'] == r)]))
            m.append(numpy.sum(data[denominator][numpy.logical_and(data[quant['key']] > cut, data['run'] == r)]))

        n = numpy.array(n)
        m = numpy.array(m)
        sigma = getError(n, m)
        p1.errorbar(x, n/m, yerr=sigma, fmt='.', label='{} > {}'.format(quant['label'], cut) , capsize=2)


    p1.set_xticks(x)
    p1.set_xticklabels(runs, rotation=90, size=6)
    p1.text(0.01, 0.01, 'dataset: {}\ndenominator: {}'.format(dataset, denominator.replace('_', ' ')), verticalalignment='bottom', horizontalalignment='left', transform=p1.transAxes)

    p1.set_title('{}'.format(trigger.replace('_v', '')))
    p1.set_xlabel('runnumber')
    p1.set_ylabel('efficiency')

    #p1.set_ylim(-0.05, 1.05)

    p1.legend(loc=4)
    p1.grid(alpha=0.5)

    #fig.tight_layout()
    out = 'res/runcomp_{}-{}-{}-{}-{}.pdf'.format(trigger.replace('_v', ''), denominator, dataset, quant['key'])
    makeSlide(out, texpath, caption='runcomparison for {} on {}'.format(trigger.replace('_v', ''), dataset))
    fig.savefig(out, dpi=300, transparent=False)
    #matplotlib.pyplot.show()
    matplotlib.pyplot.close()


#def makeEffPointPlot(data, dataset, triggers, quant, denominator, texpath, x0=[0.95, 0.005, 1000]):
    ## define bins
    #width = quant['limits'][2]
    #center = numpy.arange(quant['limits'][0], quant['limits'][1], width)
    #bins = numpy.append(center - 0.5 * width, center[-1] + 0.5 * width)

    #data = data[data[denominator]]

    #fig = matplotlib.pyplot.figure(figsize=(10, 6))
    #p1 = fig.add_subplot(211)
    #p2 = fig.add_subplot(212)
    #x = []
    #for run in data['run']:
        #if not run in x:
            #x.append(run)
    #runs = numpy.sort(x)
    #x = numpy.arange(len(runs))
    #y = []
    #z = []

    #for trig in triggers:
        #for r in runs:
            #denominatorentries = numpy.histogram(data[quant['key']][(data['run'] == r)], bins)[0]
            #entries = numpy.histogram(data[quant['key']][(data['run'] == r) & (data[trig]==1)], bins)[0]
            #xfit, yfit, label, success, para, paraerr, cut = doFit(center, entries, denominatorentries, x0=x0)
            #if success:
                #y.append(cut)
                #z.append(para[0])
            #else:
                #y.append(numpy.nan)
                #z.append(numpy.nan)

        #p1.errorbar(x, y, yerr=0, fmt='.', label='{}'.format(trig) , capsize=2)
        #p2.errorbar(x, z, yerr=0, fmt='.', label='{}'.format(trig) , capsize=2)


    #p1.set_xticks(x)
    #p1.set_xticklabels(runs, rotation=90, size=6)
    #p1.text(0.01, 0.01, 'dataset: {}\ndenominator: {}'.format(dataset, denominator.replace('_', ' ')), verticalalignment='bottom', horizontalalignment='left', transform=p1.transAxes)

    #p1.set_title('efficiency point')
    #p1.set_xlabel('runnumber')
    #p1.set_ylabel('99%-plateau')


    #p1.legend(loc=4)
    #p1.grid(alpha=0.5)

    #p1.set_xlabel('runnumber')
    #p1.set_ylabel('plateau')

    #p2.set_xticks(x)
    #p2.set_xticklabels(runs, rotation=90, size=6)

    #p2.legend(loc=4)
    #p2.grid(alpha=0.5)

    ##fig.tight_layout()
    #out = 'res/effpoint_{}-{}-{}.pdf'.format(denominator, dataset, quant['key'])
    #makeSlide(out, texpath, caption='runeffpoint on {}'.format(dataset))
    #fig.savefig(out, dpi=300, transparent=False)
    ##matplotlib.pyplot.show()
    #matplotlib.pyplot.close()


#def doEffBarPlot(data, triggers, runs, quantity, label, limits=[0, 2500], width=30, splines=False, showruns=True): #FIXME
    #center = numpy.arange(limits[0], limits[1], width)
    #bins = numpy.append(center - 0.5 * width, center[-1] + 0.5 * width)


    #trig_label = (', '.join(triggers)).replace('_v', '')
    #run_label = 'runs:'
    #for run in runs:
        #run_label = run_label + '\n{}'.format(run)

    #baselineentries = 0
    #entries = 0
    #for trig in triggers:
        #for run in numpy.sort(runs):
            #baselineentries = baselineentries + numpy.histogram(data[quantity][data['run'] == run], bins)[0]
            #entries = entries + numpy.histogram(data[quantity][(data['run'] == run) & (data[trig]==1)], bins)[0]

    ## begin plot
    #fig = matplotlib.pyplot.figure(figsize=(10, 2 * 6))

    #p1 = fig.add_subplot(211)

    #p1.bar(center, baselineentries, width=width, label='baseline\nentries: {}'.format(numpy.sum(baselineentries)))[0]
    #p1.bar(center, entries, width=width, label='{}\nentries: {}'.format(trig_label, numpy.sum(entries)))[0]


    #p1.set_title(trig_label)
    #p1.set_xlabel(label)
    #p1.set_ylabel('entries')

    #p1.legend(loc=0)


    #p2 = fig.add_subplot(212, sharex=p1)

    #eff = entries/baselineentries
    #mask = numpy.isfinite(eff)
    #center = center[mask]
    #entries = entries[mask]
    #baselineentries = baselineentries[mask]
    #eff = eff[mask]

    #sigma = get_error(entries, baselineentries)


    #p2.errorbar(center, eff, yerr=sigma, fmt='.', label=trig_label)

    ##splines
    #if(splines):
        #xfit = numpy.linspace(min(center), max(center), 1000)
        #fitmask = (center >= numpy.min(xfit)) & (center <= numpy.max(xfit))
        #tck = scipy.interpolate.splrep(center, eff)
        #yfit = scipy.interpolate.splev(xfit, tck, der=0)
        #p2.plot(xfit, yfit, label='cubic splines')

    #if showruns:
        #p1.text(0.01, 0.99, run_label, verticalalignment='top', horizontalalignment='left', transform=p1.transAxes, size=6)


    #p2.set_xlim(limits[0], limits[1])

    #p2.set_xlabel(label)
    #p2.set_ylabel('efficiency')

    #p2.legend(loc=0)

    #fig.tight_layout()
    #fig.savefig('eff/{}_{}.pdf'.format(trig_label, quantity), dpi=300, transparent=False)
    #matplotlib.pyplot.close()