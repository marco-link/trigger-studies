## @package triggerplotModule
# module to generate trigger efficiency plots
# @note requires numpy, scipy and matplotlib


## @var fitthresh
# y-threshhold for datapoints to be considered in fitting
fitthresh = 0.8

## @var worklabel
# label shown in the top left corner of generated plots
worklabel = 'CMS private work'



import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import scipy.optimize
import scipy.interpolate
import scipy.stats
import scipy.special

import modules.dataModule



## clears a file; usually called before generating plots.
# @param path   filepath of the file to clear
def clearfile(path):
    # delete file at beginning
    f = open(path, 'w')
    f.close()


## writes text to textfile
# @param txt        text to write in file
# @param textpath   path of textfile
def write2tex(txt, texpath):
    f = open(texpath, 'a')
    f.write('{}\n'.format(txt))
    f.close()


## writes LaTeX formated slide with graphic to textfile
# @param name       filename of the graphic
# @param texpath    filepath of textfile
# @param caption    (optional) caption of the slide
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


## function to calculate the asymmetric error for the trigger efficiency 
# using Clopper-Pearson interval like defined in
# https://de.wikipedia.org/wiki/Konfidenzintervall_f%C3%BCr_die_Erfolgswahrscheinlichkeit_der_Binomialverteilung
# @param k      number of hits
# @param n      number of experiments
# @param gamma  (optional) confidence level for interval; default is one sigma (68,2%)
#
# @retval list  containing the lower and upper error
def getError(k, n, gamma=0.682):
    alpha = 0.5 * (1-gamma)
    lower = k/n - scipy.stats.beta.ppf(alpha, k, n-k+1)
    upper = scipy.stats.beta.ppf(1-alpha, k+1, n-k) - k/n

    lower[k==0] = 0
    upper[k==n] = 0

    return [lower, upper]


## calculates efficiency
# @param data           TODO
# @param trigger        trigger TODO
# @param quant          quantity used for x-axis
# @param denominator    denominator used for filtering
#
# @retval numpy_array   center TODO
# @retval numpy_array   eff TODO
# @retval numpy_array   sigma TODO
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


## generates a trigger efficiency plot
# @param dataset    TODO
# @param trigger    TODO
# @param quant      TODO
# @param texpath    TODO
# @param fit        (optional) TODO
# @param x0       (optional) TODO
# @param mask       (optional) TODO
def doEffPlot(dataset, trigger, quant, texpath, fit=False, x0=[0.9, 0.005, 1000], mask=''):
    data = dataset['data']
    lumis = []
    for dtset in dataset['sets']:
        lumis.append(modules.dataModule.getLumi(data[data[dataset['key']] == dtset]))

    if not mask=='':
        data = data[data[mask] == 1]

    fig = matplotlib.pyplot.figure(figsize=(10, 6))
    p1 = fig.add_subplot(111)

    p1.text(0.01, 0.99, worklabel, verticalalignment='top', horizontalalignment='left', transform=p1.transAxes, size = 12, fontweight="bold")

    for trig in trigger:
        for dtset, lumi in zip(dataset['sets'], lumis):
            x, y, sigma = getEfficiency(data[data[dataset['key']] == dtset], trig, quant, dataset['denom'])

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
    p1.set_title(label, size=12)
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


## fits a modified errof function to data
# @param xdata      TODO
# @param ydata      TODO
# @param sigma      TODO
# @param x0         TODO
# @param cuteff     (optional) TODO
#
# @retval           TODO bla
# @retval           TODO bla
# @retval           TODO bla
# @retval           TODO bla
# @retval           TODO bla
# @retval           TODO bla
# @retval           TODO bla
def doFit(xdata, ydata, sigma, x0, cuteff=0.99):
    # only fit on efficiency > fitthresh on last entries connected
    mask = []
    y = True
    for x in numpy.flip(ydata > fitthresh, 0):
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

    # error function
    def fitfunc(x, para):
        a, b, c = para
        return 0.5 * a * ( 1 + scipy.special.erf(b * (x - c)))

    # logistic function
    #def fitfunc(x, para):
        #a, b, c = para
        #return a / (1 + numpy.exp(-b * (x - c)))

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

    label = '$y = \\frac{{a}}{{2}} [1 + erf(b (x - c))]$\n$a={:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$\n$b={:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$\n$c={:.0f}^{{+{:.1f}}}_{{-{:.1f}}}$\n$\chi^2 / dof = {:.3f}/{}$\np-value: {:.2g}\n$y>${:.2f}$a$: $x>${:.1f}\n'.format(para[0], *paraerr[0], para[1], *paraerr[1], para[2], *paraerr[2], Chi2, dof, pvalue, cuteff, cut)
    if(len(xdata)>0):
        xfit = numpy.linspace(numpy.amin(xdata), numpy.amax(xdata), 1000)
    else:
        xfit = -999

    return xfit, fitfunc(xfit, para), label, res.success, para, paraerr, cut


## generates 2D trigger efficiency plots
# @param dataset    TODO
# @param trigger    TODO
# @param quant1     TODO
# @param quant2     TODO
# @param texpath    TODO
# @param cuts       (optional) TODO
# @param mask       (optional) TODO
def do2DPlot(dataset,  trigger, quant1, quant2, texpath, cuts=None, mask=''):
    data = dataset['data']
    datasets = dataset['sets']
    denominator = dataset['denom']
    lumi = modules.dataModule.getLumi(data[data[dataset['key']].isin(datasets).values])
    if not mask=='':
        data = data[data[mask] == 1]


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
    denominatorentries = numpy.histogram2d(data[quant1['key']], data[quant2['key']], bins=[bins1, bins2])[0]
    entries = numpy.histogram2d(data[quant1['key']][data[trigger]==1], data[quant2['key']][data[trigger]==1], bins=[bins1, bins2])[0]

    # plotting
    fig = matplotlib.pyplot.figure(figsize=(10,6))

    p1 = fig.add_subplot(111)

    p1.text(0.01, 0.99, worklabel, verticalalignment='top', horizontalalignment='left', transform=p1.transAxes, size = 12, fontweight="bold", backgroundcolor= (1., 1., 1., 0.6))

    eff = entries/denominatorentries
    eff[denominatorentries < 3] = numpy.nan
    eff = numpy.ma.masked_invalid(eff)

    cmap = matplotlib.pyplot.get_cmap('viridis')
    cmap.set_bad(color='w', alpha=0.8)

    fig.colorbar(p1.pcolormesh(bins1, bins2, eff.T, vmin=0, cmap=cmap))

    if mask == '':
        p1.set_title('{}'.format(trigger.replace('_v', '')), fontweight="bold", size=14)
    else:
        p1.set_title('{} ({})'.format(trigger.replace('_v', ''), mask), fontweight="bold", size=14)
    p1.set_xlabel(quant1['label'])
    p1.set_ylabel(quant2['label'])
    p1.text(0.01, 0.01, 'denominator: {}\ndatasets: {} ({:.2f} fb$^{{-1}}$)'.format(denominator.replace('_', ' '), ', '.join(datasets), lumi), verticalalignment='bottom', horizontalalignment='left', transform=p1.transAxes, backgroundcolor= (1., 1., 1., 0.6))

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