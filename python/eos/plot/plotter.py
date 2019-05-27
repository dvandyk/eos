# Copyright (c) 2018 Frederik Beaujean
# Copyright (c) 2017, 2018 Danny van Dyk
#
# This file is part of the EOS project. EOS is free software;
# you can redistribute it and/or modify it under the terms of the GNU General
# Public License version 2, as published by the Free Software Foundation.
#
# EOS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA  02111-1307  USA

import eos
from logging import debug, info, warn
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import gaussian_kde
import sys

""" Plotter it used to produce publication quality plots with EOS.

Plotter uses matplotlib to produce EOS plots within PDF files. It is the backend
for the eos-plot-* scripts.
"""
class Plotter:
    def __init__(self, instructions, output):
        """
        Parameters
        ----------

        instructions : Dictionary containing the instructions on what to plot in which manner.
        output       : Name of the output PDF file.
        """
        self.instructions = instructions
        self.output = output
        self.fig = None
        self.ax = None
        self.xrange = None
        self.yrange = None


    """ Setting up the plot based on the provided instruction. """
    def setup_plot(self):
        if not 'plot' in self.instructions:
            raise KeyError('no plot metadata specified')

        myplot = self.instructions['plot']

        self.fig, self.ax = plt.subplots()

        mytitle = ''
        myylabel = ''
        myxlabel = ''
        if 'title' in myplot:
            mytitle = myplot['title']

        if 'size' in myplot:
            xwidth, ywidth = myplot['size']
            # size is specified in cm, matplotlib expects inches
            # convert from cm to inches
            xwidth /= 2.54 # cm / inch
            ywidth /= 2.54 # cm / inch
            plt.gcf().set_size_inches((xwidth, ywidth))

        plt.locator_params(axis='x', nbins=5)
        plt.locator_params(axis='y', nbins=5)

        if 'x' in myplot:
            myx = myplot['x']

            if 'label' in myx:
                myxlabel = myx['label']

            if 'unit' in myx:
                myxlabel += r'\,[' + myx['unit'] + r']'

            if 'range' in myx:
                self.xrange = myx['range']
                self.ax.set_xlim(tuple(self.xrange))

            if 'scale' in myx:
                self.xscale = float(myx['scale'])
                self.xticks = matplotlib.ticker.FuncFormatter(lambda x, pos, xscale=self.xscale: '${0:.2f}$'.format(x / xscale))
                self.ax.xaxis.set_major_formatter(self.xticks)

        if 'y' in myplot:
            myy = myplot['y']

            if 'label' in myy:
                myylabel = myy['label']

            if 'unit' in myy:
                myylabel += r'\,[' + myy['unit'] + r']'

            if 'range' in myy:
                self.yrange = myy['range']
                self.ax.set_ylim(tuple(self.yrange))

            if 'scale' in myy:
                self.yscale = float(myy['scale'])
                self.yticks = matplotlib.ticker.FuncFormatter(lambda y, pos, yscale=self.yscale: '${0:.2f}$'.format(y / yscale))
                self.ax.yaxis.set_major_formatter(self.yticks)

        self.ax.set(xlabel=myxlabel, ylabel=myylabel, title=mytitle)


    """ Plots a single EOS observable w/o uncertainties as a function of one kinemtic variable or one parameter. """
    class ObservablePlot:
        def __init__(self, plotter, item):
            self.oname = item['observable']
            info('   plotting EOS observable "{}"'.format(oname))

            # create parameters
            self.parameters = eos.Parameters.Defaults()
            if 'parameters' in item and 'parameters-from-file' in item:
                warn('    overriding values read from \'parameters-from-file\' with explicit values in \'parameters\'')

            if 'parameters-from-file' in item and type(item['parameters-from-file']) is str:
                self.parameters.override_from_file(item['parameters-from-file'])

            if 'parameters' in item and type(item['parameters']) is dict:
                for key, value in item['parameters'].items():
                    self.parameters.set(key, value)

            # create kinematics
            self.kinematics = eos.Kinematics()
            if not 'kinematic' in item and not 'parameter' in item:
                raise KeyError('neither kinematic nor parameter found; do not know how to map x to a variable')
            if 'kinematic' in item and 'parameter' in item:
                raise KeyError('both kinematic and parameter found; do not know how to map x to a variable')
            if 'kinematic' in item:
                self.var = self.kinematics.declare(item['kinematic'], np.nan)
            else:
                self.var = self.parameters.declare(item['parameter'], np.nan)
                if 'kinematics' in item:
                    for k, v in item['kinematics'].items():
                        self.kinematics.declare(k, v)

            # create (empty) options
            self.options = eos.Options()

            # create observable
            observable = eos.Observable.make(self.oname, self.parameters, self.kinematics, self.options)

            # determine plot settings
            self.color         = item['color']   if 'color'   in item else 'black'
            self.samples       = item['samples'] if 'samples' in item else 100
            self.xlo, self.xhi = item['range']   if 'range'   in item else self.xrange


        def plot(self):
            xvalues = np.linspace(self.xlo, self.xhi, self.samples + 1)
            ovalues = np.array([])
            for xvalue in xvalues:
                self.var.set(xvalue)
                ovalues = np.append(ovalues, self.observable.evaluate())

            plt.plot(xvalues, ovalues, color=self.color)

        @staticmethod
        def make(plot, item):
            return(PlotObservable(plotter, item))

    """ Plots an uncertainty band as a function of one kinematic variable.

    This routine expects the uncertainty propagation to have produces an HDF5 file.
    """
    def plot_uncertainty(self, item):
        if 'hdf5-file' not in item:
            raise KeyError('no hdf5-file specified')
            return

        h5fname = item['hdf5-file']
        info('   plotting uncertainty propagation from file "{}"'.format(h5fname))

        uncfile = eos.data.UncertaintyDataFile(h5fname)
        _xvalues = []
        for o in uncfile.parameters:
            kin = o[1].split(b',')
            if len(kin) > 1:
                raise ValueError('more than one kinematic variable specified')

            name,value = kin[0].split(b'=')
            value = float(value)
            _xvalues.append(float(value))

        _xvalues = np.array(_xvalues)
        if 'range' in item:
            xmin,xmax = item['range']
            _xvalues = np.ma.masked_outside(_xvalues, float(xmin), float(xmax))

        data = uncfile.data()
        _ovalues_lower   = []
        _ovalues_central = []
        _ovalues_higher  = []
        for i in range(len(uncfile.parameters)):
            lower   = np.percentile(data[:, i], q=15.865)
            central = np.percentile(data[:, i], q=50.000)
            higher  = np.percentile(data[:, i], q=84.135)
            _ovalues_lower.append(lower)
            _ovalues_central.append(central)
            _ovalues_higher.append(higher)

        color = 'black'
        if 'color' in item:
            color = item['color']

        alpha = 1.0
        if 'opacity' in item:
            alpha = item['opacity']

        label = item['label'] if 'label' in item else None

        xvalues = np.linspace(np.min(_xvalues),np.max(_xvalues),100)
        # work around CubicSpline missing in SciPy version < 0.18
        if scipy.__version__ >= '0.18':
            from scipy.interpolate import CubicSpline
            interpolate = lambda x, y, xv: CubicSpline(x, y)(xv)
        else:
            from scipy.interpolate import spline
            interpolate = spline

        ovalues_lower   = interpolate(_xvalues, _ovalues_lower,   xvalues)
        ovalues_central = interpolate(_xvalues, _ovalues_central, xvalues)
        ovalues_higher  = interpolate(_xvalues, _ovalues_higher,  xvalues)

        plt.fill_between(xvalues, ovalues_lower, ovalues_higher, lw=0, color=color, alpha=alpha, label=label)
        plt.plot(xvalues, ovalues_lower,   color=color, alpha=alpha)
        plt.plot(xvalues, ovalues_central, color=color, alpha=alpha)
        plt.plot(xvalues, ovalues_higher,  color=color, alpha=alpha)


    """ Plots constraints from the EOS library of experimental and theoretical likelihoods. """
    def plot_constraint(self, item):
        import yaml

        constraints = eos.Constraints()

        if 'constraints' not in item:
            raise KeyError('no constraints specified')

        names = item['constraints']
        if type(names) == str:
            names = [names]

        for name in names:
            entry = constraints[name]
            if not entry:
                raise ValueError('unknown constraint {}'.format(name))

            constraint = yaml.load(entry.serialize())

            xvalues = None
            xerrors = None
            yvalues = None
            yerrors = None

            if constraint['type'] == 'Gaussian':
                kinematics = constraint['kinematics']
                width = 1
                if item['variable'] in kinematics:
                    xvalues = [kinematics[item['variable']]]
                    xerrors = None
                elif (item['variable'] + '_min' in kinematics) and (item['variable'] + '_max' in kinematics):
                    xvalues = [(kinematics[item['variable'] + '_max'] + kinematics[item['variable'] + '_min']) / 2]
                    xerrors = [(kinematics[item['variable'] + '_max'] - kinematics[item['variable'] + '_min']) / 2]
                    if item['rescale-by-width']:
                        width = (kinematics[item['variable'] + '_max'] - kinematics[item['variable'] + '_min'])

                yvalues = [constraint['mean'] / width]
                sigma_hi = np.sqrt(float(constraint['sigma-stat']['hi'])**2 + float(constraint['sigma-sys']['hi'])**2) / width
                sigma_lo = np.sqrt(float(constraint['sigma-stat']['lo'])**2 + float(constraint['sigma-sys']['lo'])**2) / width
                yerrors = [[sigma_hi, sigma_lo]]
            elif constraint['type'] == 'MultivariateGaussian(Covariance)':
                if 'observable' not in item:
                    raise KeyError('observable needs to be specified for MultivariateGaussian(Covariance) constraints')
                dim = constraint['dim']
                covariance = np.array(constraint['covariance'])
                observables = constraint['observables']
                means = constraint['means']
                kinematics = constraint['kinematics']

                xvalues = []
                xerrors = []
                yvalues = []
                yerrors = []
                for i in range(0, dim):
                    width = 1

                    if not observables[i] == item['observable']:
                        continue
                    _kinematics = kinematics[i]
                    if item['variable'] in _kinematics:
                        xvalues.append(_kinematics[item['variable']])
                        xerrors.append(0)
                    elif (item['variable'] + '_min' in _kinematics) and (item['variable'] + '_max' in _kinematics):
                        xvalues.append((_kinematics[item['variable'] + '_max'] + _kinematics[item['variable'] + '_min']) / 2)
                        xerrors.append((_kinematics[item['variable'] + '_max'] - _kinematics[item['variable'] + '_min']) / 2)
                        if item['rescale-by-width']:
                            width = (_kinematics[item['variable'] + '_max'] - _kinematics[item['variable'] + '_min'])

                    yvalues.append(means[i] / width)
                    yerrors.append(np.sqrt(covariance[i, i]) / width)
            elif constraint['type'] == 'MultivariateGaussian':
                if 'observable' not in item:
                    raise KeyError('observable needs to be specified for MultivariateGaussian constraints')
                dim = constraint['dim']
                sigma_stat_hi = np.array(constraint['sigma-stat-hi'])
                sigma_stat_lo = np.array(constraint['sigma-stat-lo'])
                sigma_sys = np.array(constraint['sigma-sys'])
                sigma = np.sqrt(np.power(sigma_sys, 2) + 0.25 * np.power(sigma_stat_hi + sigma_stat_lo, 2))
                observables = constraint['observables']
                means = constraint['means']
                kinematics = constraint['kinematics']

                xvalues = []
                xerrors = []
                yvalues = []
                yerrors = []
                for i in range(0, dim):
                    width = 1

                    if not observables[i] == item['observable']:
                        continue
                    _kinematics = kinematics[i]
                    if item['variable'] in _kinematics:
                        xvalues.append(_kinematics[item['variable']])
                        xerrors.append(0)
                    elif (item['variable'] + '_min' in _kinematics) and (item['variable'] + '_max' in _kinematics):
                        xvalues.append((_kinematics[item['variable'] + '_max'] + _kinematics[item['variable'] + '_min']) / 2)
                        xerrors.append((_kinematics[item['variable'] + '_max'] - _kinematics[item['variable'] + '_min']) / 2)
                        if item['rescale-by-width']:
                            width = (_kinematics[item['variable'] + '_max'] - _kinematics[item['variable'] + '_min'])

                    yvalues.append(means[i] / width)
                    yerrors.append(sigma[i] / width)
            else:
                raise ValueError('type of constraint presently not supported')

            xvalues = np.array(xvalues)
            if xerrors:
                xerrors = np.array(xerrors)
            yvalues = np.array(yvalues)
            yerrors = np.array(yerrors)

            color = item['color'] if 'color' in item else 'black'
            label = item['label'] if 'label'   in item else None

            plt.errorbar(x=xvalues, y=yvalues, xerr=xerrors, yerr=yerrors.T,
                color=color, elinewidth=1.0, fmt='_', linestyle='none', label=label)


    """ Plots 2D contours of a pair of parameters based on pre-existing random samples. """
    def plot_contours2d(self, item):
        if 'hdf5-file' not in item:
            raise KeyError('no hdf5-file specified')

        h5fname = item['hdf5-file']
        info('   plotting 2D contours from file "{}"'.format(h5fname))
        datafile = eos.data.load_data_file(h5fname)

        if 'variables' not in item:
            raise KeyError('no variables specificed')

        xvariable, yvariable = item['variables']

        data   = datafile.data()
        xindex = datafile.variable_indices[xvariable]
        xdata  = data[:, xindex]
        yindex = datafile.variable_indices[yvariable]
        ydata  = data[:, yindex]

        if not np.array(self.xrange).any():
            self.xrange = [np.amin(xdata), np.amax(xdata)]
            self.ax.set_xlim(tuple(self.xrange))
        if not np.array(self.yrange).any():
            self.yrange = [np.amin(ydata), np.amax(ydata)]
            self.ax.set_ylim(tuple(self.yrange))
        plt.show()

        xbins = 100
        ybins = 100

        H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(xbins, ybins), normed=True)
        x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,xbins))
        y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((ybins,1))
        pdf = (H*(x_bin_sizes*y_bin_sizes))

        # find the PDF value corresponding to a given cummulative probability
        plevel = lambda x, pdf, P: pdf[pdf > x].sum() - P
        pone_sigma = scipy.optimize.brentq(plevel, 0., 1., args=(pdf, 0.68))
        ptwo_sigma = scipy.optimize.brentq(plevel, 0., 1., args=(pdf, 0.95))
        pthree_sigma = scipy.optimize.brentq(plevel, 0., 1., args=(pdf, 0.99))
        levels = [pone_sigma, ptwo_sigma, pthree_sigma]
        labels = ['68%', '95%', '99%']

        CS = plt.contour(pdf.transpose(),
                         colors='OrangeRed',
                         extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                         levels=levels[::-1])

        fmt = {}
        for level, label in zip(CS.levels, labels[::-1]):
            fmt[level] = label

        plt.clabel(CS, inline=1, fmt=fmt, fontsize=10)


    """ Plots a 1D Kernel Density Estimate (KDE) of pre-existing random samples. """
    def plot_kde(self, item):
        if 'hdf5-file' not in item:
            raise KeyError('no hdf5-file specified')

        h5fname = item['hdf5-file']
        info('   plotting histogram from file "{}"'.format(h5fname))
        datafile = eos.data.load_data_file(h5fname)

        if 'variable' not in item:
            raise KeyError('no variable specificed')
        variable = item['variable']

        print(datafile.variable_indices)
        if variable not in datafile.variable_indices:
            raise ValueError('variable {} not contained in data file'.format(variable))

        index = datafile.variable_indices[variable]
        data  = datafile.data()[:, index]
        alpha = item['opacity']   if 'opacity'   in item else 0.3
        color = item['color']     if 'color'     in item else 'blue'
        bw    = item['bandwidth'] if 'bandwidth' in item else None

        kde = gaussian_kde(data)
        kde.set_bandwidth(bw_method='silverman')
        if 'bandwidth' in item:
            kde.set_bandwidth(bw_method=kde.factor * item['bandwidth'])

        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 1000)

        plt.plot(x, kde(x), color=color)


    """ Plots a 1D histogram of pre-existing random samples. """
    def plot_histogram(self, item):
        if 'hdf5-file' not in item:
            raise KeyError('no hdf5-file specified')

        h5fname = item['hdf5-file']
        info('   plotting histogram from file "{}"'.format(h5fname))
        datafile = eos.data.load_data_file(h5fname)

        if 'variable' not in item:
            raise KeyError('no variable specificed')
        variable = item['variable']

        print(datafile.variable_indices)
        if variable not in datafile.variable_indices:
            raise ValueError('variable {} not contained in data file'.format(variable))

        index = datafile.variable_indices[variable]
        data  = datafile.data()[:, index]
        alpha = item['opacity'] if 'opacity' in item else 0.3
        color = item['color']   if 'color'   in item else 'blue'
        bins  = item['bins']    if 'bins'    in item else 100
        plt.hist(data, alpha=alpha, bins=bins, color=color, density=1)


    """ Plots a 2D histogram of pre-existing random samples. """
    def plot_histogram2d(self, item):
        if 'hdf5-file' not in item:
            raise KeyError('no hdf5-file specified')

        h5fname = item['hdf5-file']
        info('   plotting 2D histogram from file "{}"'.format(h5fname))
        datafile = eos.data.load_data_file(h5fname)

        if 'variables' not in item:
            raise KeyError('no variables specificed')

        xvariable, yvariable = item['variables']

        xindex = datafile.variable_indices[xvariable]
        yindex = datafile.variable_indices[yvariable]
        data = datafile.data()

        cmap = plt.get_cmap('viridis')
        cmap.set_under('w', 1)

        xdata = data[:, xindex]
        ydata = data[:, yindex]
        bins  = item['bins']    if 'bins'    in item else 100
        plt.hist2d(xdata, ydata, bins=bins, cmin=1)


    def plot_function(self, item):
        if 'f' not in item:
            raise KeyError('no function specificied')
        f = item['f']
        alpha  = item['opacity'] if 'opacity' in item else 1.0
        color  = item['color']   if 'color'   in item else 'black'
        style  = item['style']   if 'style'   in item else '-'
        points = item['points']  if 'points'  in item else 100
        label  = item['label']   if 'label'   in item else None

        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, points)
        y = []

        for xvalue in x:
            y.append(eval(f, {}, {'x': xvalue}))

        plt.plot(x, y, color=color, alpha=alpha, linestyle=style, label=label)


    """ Inserts an EOS watermark into the plots. """
    class WatermarkPlot:
        def __init__(self, plotter, item):
            xdelta, ydelta = (0.04, 0.04)

            hpos, vpos = item['position'] if 'position' in item else ['right', 'top']

            if hpos == 'right':
                self.x = 1 - xdelta
            elif hpos == 'left':
                self.x = xdelta
            elif hpos == 'center':
                self.x = 0.5
            else:
                raise ValueError('invalid horizontal position \'{}\''.format(hpos))

            if vpos == 'bottom':
                self.y = 0 + ydelta
            elif vpos == 'top':
                self.y = 1 - ydelta
            elif vpos == 'center':
                self.y = 0.5
            else:
                raise ValueError('invalid vertical position \'{}\''.format(hpos))

            if 'preliminary' in item and item['preliminary']:
                self.color = 'red'
                self.version = 'Preliminary'
            else:
                self.color = 'OrangeRed'
                self.version = 'v{version}'.format(version=eos.version())

        def plot(self):
            logofont = matplotlib.font_manager.FontProperties(family='sans-serif', size='15')
            ax = plt.gca()
            ax.text(self.x, self.y, r'\textsf{{\textbf{{EOS {version}}}}}'.format(version=self.version),
                    transform=ax.transAxes, fontproperties=logofont,
                    color=color, alpha=0.5, bbox=dict(facecolor='white', alpha=0.5, lw=0),
                    horizontalalignment=hpos, verticalalignment=vpos, zorder=+5)

        @staticmethod
        def make(plotter, item):
            return(WatermarkPlot(plotter, item))


    """ Plots the contents specified in the instructions provided to Plotter. """
    def plot_contents(self):
        if not 'contents' in self.instructions:
            return

        plot_handlers = {
            #'constraint':  Plotter.plot_constraint,
            #'contours2D':  Plotter.plot_contours2d,
            #'function':    Plotter.plot_function,
            #'histogram':   Plotter.plot_histogram,
            #'histogram2D': Plotter.plot_histogram2d,
            #'kde':         Plotter.plot_kde,
            'observable':  Plotter.plot_observable,
            #'uncertainty': Plotter.plot_uncertainty,
            'watermark':   Plotter.plot_eos_watermark,
        }

        anonymous_types = [
            'watermark',
        ]

        contents = self.instructions['contents']

        for item in contents:
            if not type(item) is dict:
                TypeError('wrong data type for content item {}'.format(str(item)))

            name = item['name'] if 'name' in item else 'unnamed'
            if not 'type' in item:
                raise KeyError('plot content "{}" has no type'.format(name))
            item_type = item['type']

            if item_type not in anonymous_types and 'name' not in item:
                raise KeyError('unnamed plot content')
            elif 'name' not in item:
                name = None
                debug('plotting anonymous contents of type \'{}\''.format(item_type))
            else:
                name = item['name']
                info('plotting "{}"'.format(name))

            if item_type not in plot_functions:
                KeyError('unknown content type: "{}"'.format(item_type))

            plot_functions[item_type](self, item)

        if 'legend' in self.instructions['plot']:
            if 'location' in self.instructions['plot']['legend']:
                self.ax.legend(loc=self.instructions['plot']['legend']['location'])


    """ Produces the plot. """
    def plot(self):
        self.setup_plot()
        self.plot_contents()

        plt.show()
        plt.savefig(self.output)


def variable_to_latex(variable):
    p = eos.Parameters.Defaults()
    if variable in [pp.name() for pp in p]:
        return p[variable].latex()
    else:
        return r'\verb{' + variable + '}'


