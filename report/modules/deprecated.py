# DEPRECATED code

print('using deprecated code!\nbe carefull\n\n\n')

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

def makeRunPlot(data, dataset, trigger, quant, denominator, cuts, texpath): #FIXME
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

    p1.set_title('{}'.format(trigger.replace('_v', '')), fontweight="bold", size=14)
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


def makeEffPointPlot(data, dataset, triggers, quant, denominator, texpath, x0=[0.95, 0.005, 1000]): #FIXME
    # define bins
    width = quant['limits'][2]
    center = numpy.arange(quant['limits'][0], quant['limits'][1], width)
    bins = numpy.append(center - 0.5 * width, center[-1] + 0.5 * width)

    data = data[data[denominator]]

    fig = matplotlib.pyplot.figure(figsize=(10, 6))
    p1 = fig.add_subplot(211)
    p2 = fig.add_subplot(212)
    x = []
    for run in data['run']:
        if not run in x:
            x.append(run)
    runs = numpy.sort(x)
    x = numpy.arange(len(runs))
    y = []
    z = []

    for trig in triggers:
        for r in runs:
            denominatorentries = numpy.histogram(data[quant['key']][(data['run'] == r)], bins)[0]
            entries = numpy.histogram(data[quant['key']][(data['run'] == r) & (data[trig]==1)], bins)[0]
            xfit, yfit, label, success, para, paraerr, cut = doFit(center, entries, denominatorentries, x0=x0)
            if success:
                y.append(cut)
                z.append(para[0])
            else:
                y.append(numpy.nan)
                z.append(numpy.nan)

        p1.errorbar(x, y, yerr=0, fmt='.', label='{}'.format(trig) , capsize=2)
        p2.errorbar(x, z, yerr=0, fmt='.', label='{}'.format(trig) , capsize=2)


    p1.set_xticks(x)
    p1.set_xticklabels(runs, rotation=90, size=6)
    p1.text(0.01, 0.01, 'dataset: {}\ndenominator: {}'.format(dataset, denominator.replace('_', ' ')), verticalalignment='bottom', horizontalalignment='left', transform=p1.transAxes)

    p1.set_title('efficiency point')
    p1.set_xlabel('runnumber')
    p1.set_ylabel('99%-plateau')


    p1.legend(loc=4)
    p1.grid(alpha=0.5)

    p1.set_xlabel('runnumber')
    p1.set_ylabel('plateau')

    p2.set_xticks(x)
    p2.set_xticklabels(runs, rotation=90, size=6)

    p2.legend(loc=4)
    p2.grid(alpha=0.5)

    #fig.tight_layout()
    out = 'res/effpoint_{}-{}-{}.pdf'.format(denominator, dataset, quant['key'])
    makeSlide(out, texpath, caption='runeffpoint on {}'.format(dataset))
    fig.savefig(out, dpi=300, transparent=False)
    #matplotlib.pyplot.show()
    matplotlib.pyplot.close()


def doEffBarPlot(data, triggers, runs, quantity, label, limits=[0, 2500], width=30, splines=False, showruns=True): #FIXME
    center = numpy.arange(limits[0], limits[1], width)
    bins = numpy.append(center - 0.5 * width, center[-1] + 0.5 * width)


    trig_label = (', '.join(triggers)).replace('_v', '')
    run_label = 'runs:'
    for run in runs:
        run_label = run_label + '\n{}'.format(run)

    baselineentries = 0
    entries = 0
    for trig in triggers:
        for run in numpy.sort(runs):
            baselineentries = baselineentries + numpy.histogram(data[quantity][data['run'] == run], bins)[0]
            entries = entries + numpy.histogram(data[quantity][(data['run'] == run) & (data[trig]==1)], bins)[0]

    # begin plot
    fig = matplotlib.pyplot.figure(figsize=(10, 2 * 6))

    p1 = fig.add_subplot(211)

    p1.bar(center, baselineentries, width=width, label='baseline\nentries: {}'.format(numpy.sum(baselineentries)))[0]
    p1.bar(center, entries, width=width, label='{}\nentries: {}'.format(trig_label, numpy.sum(entries)))[0]


    p1.set_title(trig_label)
    p1.set_xlabel(label)
    p1.set_ylabel('entries')

    p1.legend(loc=0)


    p2 = fig.add_subplot(212, sharex=p1)

    eff = entries/baselineentries
    mask = numpy.isfinite(eff)
    center = center[mask]
    entries = entries[mask]
    baselineentries = baselineentries[mask]
    eff = eff[mask]

    sigma = get_error(entries, baselineentries)


    p2.errorbar(center, eff, yerr=sigma, fmt='.', label=trig_label)

    #splines
    if(splines):
        xfit = numpy.linspace(min(center), max(center), 1000)
        fitmask = (center >= numpy.min(xfit)) & (center <= numpy.max(xfit))
        tck = scipy.interpolate.splrep(center, eff)
        yfit = scipy.interpolate.splev(xfit, tck, der=0)
        p2.plot(xfit, yfit, label='cubic splines')

    if showruns:
        p1.text(0.01, 0.99, run_label, verticalalignment='top', horizontalalignment='left', transform=p1.transAxes, size=6)


    p2.set_xlim(limits[0], limits[1])

    p2.set_xlabel(label)
    p2.set_ylabel('efficiency')

    p2.legend(loc=0)

    fig.tight_layout()
    fig.savefig('eff/{}_{}.pdf'.format(trig_label, quantity), dpi=300, transparent=False)
    matplotlib.pyplot.close()