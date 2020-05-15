#!/usr/bin/env python

import UnbinnedAnalysis

import argparse
from GtBurst import IRFS
from GtBurst.Configuration import Configuration
from GtBurst.commands import gtdocountsmap
from GtBurst.commands import gtbuildxmlmodel
from GtBurst.commands import gtdolike
from GtBurst import dataHandling
from GtBurst.fast_ts_map import FastTSMap
from GtBurst.my_fits_io import pyfits
from GtBurst.likelihood_profile_writer import LikelihoodProfiler

import os
import subprocess
import glob
import numpy
import xml.etree.ElementTree as ET
import numpy as np

# Use non-interactive matplotlib backend
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def clean_intervals(tstarts, tstops, gti_starts, gti_stops):
    gtis = []

    for t1, t2 in zip(gti_starts, gti_stops):
        gtis.append(Interval(t1, t2, True))

    ints = []

    for t1, t2 in zip(tstarts, tstops):
        ints.append(Interval(t1, t2, True))

    cleaned_ints = []

    for interval in ints:

        for gti in gtis:

            if interval.overlaps_with(gti):
                intersect = interval.intersect(gti)

                cleaned_ints.append(intersect)

    new_starts = [x.start for x in cleaned_ints]
    new_stops = [x.stop for x in cleaned_ints]

    return new_starts, new_stops


class IntervalsDoNotOverlap(RuntimeError):
    pass


class IntervalsNotContiguous(RuntimeError):
    pass


class Interval(object):
    def __init__(self, start, stop, swap_if_inverted=False):

        self._start = float(start)
        self._stop = float(stop)

        # Note that this allows to have intervals of zero duration

        if self._stop < self._start:

            if swap_if_inverted:

                self._start = stop
                self._stop = start

            else:

                raise RuntimeError("Invalid time interval! TSTART must be before TSTOP and TSTOP-TSTART >0. "
                                   "Got tstart = %s and tstop = %s" % (start, stop))

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @classmethod
    def new(cls, *args, **kwargs):

        return cls(*args, **kwargs)

    def _get_width(self):

        return self._stop - self._start

    @property
    def mid_point(self):

        return (self._start + self._stop) / 2.0

    def __repr__(self):

        return " interval %s - %s (width: %s)" % (self.start, self.stop, self._get_width())

    def intersect(self, interval):
        """
        Returns a new time interval corresponding to the intersection between this interval and the provided one.

        :param interval: a TimeInterval instance
        :type interval: Interval
        :return: new interval covering the intersection
        :raise IntervalsDoNotOverlap : if the intervals do not overlap
        """

        if not self.overlaps_with(interval):
            raise IntervalsDoNotOverlap("Current interval does not overlap with provided interval")

        new_start = max(self._start, interval.start)
        new_stop = min(self._stop, interval.stop)

        return self.new(new_start, new_stop)

    def merge(self, interval):
        """
        Returns a new interval corresponding to the merge of the current and the provided time interval. The intervals
        must overlap.

        :param interval: a TimeInterval instance
         :type interval : Interval
        :return: a new TimeInterval instance
        """

        if self.overlaps_with(interval):

            new_start = min(self._start, interval.start)
            new_stop = max(self._stop, interval.stop)

            return self.new(new_start, new_stop)

        else:

            raise IntervalsDoNotOverlap("Could not merge non-overlapping intervals!")

    def overlaps_with(self, interval):
        """
        Returns whether the current time interval and the provided one overlap or not

        :param interval: a TimeInterval instance
        :type interval: Interval
        :return: True or False
        """

        if interval.start == self._start or interval.stop == self._stop:

            return True

        elif interval.start > self._start and interval.start < self._stop:

            return True

        elif interval.stop > self._start and interval.stop < self._stop:

            return True

        elif interval.start < self._start and interval.stop > self._stop:

            return True

        else:

            return False

    def to_string(self):
        """
        returns a string representation of the time interval that is like the
        argument of many interval reading funcitons

        :return:
        """

        return "%f-%f" % (self.start, self.stop)

    def __eq__(self, other):

        if not isinstance(other, Interval):

            # This is needed for things like comparisons to None or other objects.
            # Of course if the other object is not even a TimeInterval, the two things
            # cannot be equal

            return False

        else:

            return self.start == other.start and self.stop == other.stop


def printCommand(cmdname, targs):
    commandLine = cmdname

    for k, v in targs.items():
        commandLine += " %s='%s'" % (k, v)
    pass

    print(("-> %s" % commandLine))


pass


def writeSourceListToFile(sources, outfile):
    colNames = ['name', 'ra', 'dec', 'tstart', 'tstop', 'TS', 'photonIndex', 'photonIndexError',
                'flux', 'fluxError', 'photonFlux', 'photonFluxError', 'roi', 'irf', 'zmax', 'thetamax', 'strategy']

    with open(outfile, 'w+') as f:
        f.write("#%s" % (" ".join(colNames)))
        f.write("\n")
        for src in sources:
            vals = []
            for name in colNames:
                vals.append(str(src.__getattribute__(name)).replace(' ', ''))
            pass
            f.write(" ".join(vals))
            f.write("\n")
        pass
    pass


pass


# def get_likelihood_profile(like, norm):
#
#     # Get the current settings of the parameter so we can restore it later
#     old_scale = float(norm.attrib['scale'])
#     old_value = float(norm.attrib['value'])
#     old_min = norm.attrib['min']
#     old_max = norm.attrib['max']
#
#     # Temporarily set the scale to 1.0 so we can modify boundaries and so on without incurring in errors
#     like['GRB'].src.spectrum().parameter(norm.attrib['name']).setScale(1.0)
#
#     best_fit_value = float(old_value) * old_scale
#
#     min_value = best_fit_value / 100.0
#     max_value = best_fit_value * 100.0
#
#     min_log_value = np.log10(min_value)
#     max_log_value = np.log10(max_value)
#
#     # Make sure we are at the best fit
#
#     l0 = like.fit(verbosity=0, covar=False)
#
#     # Store best fit parameters
#     best_fit_values = []
#
#     for parameter in like.model.params:
#
#         best_fit_values.append(parameter.getValue())
#
#     def lnprob(log_norm):
#
#         if log_norm < min_log_value or log_norm > max_log_value:
#
#             # Outside of boundaries
#
#             return -np.inf
#
#         else:
#
#             # Restore best fit parameters first
#             for i, parameter in enumerate(like.model.params):
#
#                 if parameter.getName() == norm.attrib['name']:
#
#                     # This will be set below
#                     continue
#
#                 parameter.setValue(best_fit_values[i])
#
#             like['GRB'].src.spectrum().parameter(norm.attrib['name']).setFree(True)
#
#             like['GRB'].src.spectrum().parameter(norm.attrib['name']).setValue(float(10**log_norm))
#
#             like['GRB'].src.spectrum().parameter(norm.attrib['name']).setFree(False)
#
#             like.syncSrcParams()
#
#             log_like = -like.fit(verbosity=0, covar=False)
#
#             like['GRB'].src.spectrum().parameter(norm.attrib['name']).setFree(True)
#
#             return log_like
#
#     # Make sure we start from the best fit, and with the appropriate boundaries
#
#     # Temporarily make the boundaries very large, so the next assignment will not fail
#     like['GRB'].src.spectrum().parameter(norm.attrib['name']).setBounds(0.0,
#                                                                         1e9)
#
#     like['GRB'].src.spectrum().parameter(norm.attrib['name']).setValue(best_fit_value)
#
#     # Restore the boundaries
#     like['GRB'].src.spectrum().parameter(norm.attrib['name']).setBounds(min_value,
#                                                                         max_value)
#
#     n_walkers = 30
#
#     # Start all the walkers close to the best fit value
#     p0 = [[best_fit_value + np.random.uniform(-best_fit_value / 50.0, best_fit_value / 50.0)] for i in range(n_walkers)]
#
#     sampler = emcee.EnsembleSampler(n_walkers, 1, lnprob)
#
#     _ = sampler.run_mcmc(np.log10(p0), 100)
#
#     log_norm_values = sampler.flatchain
#     log_like_values = sampler.flatlnprobability
#
#     # Remove infinite values
#     idx = np.isfinite(log_like_values)
#     log_like_values = log_like_values[idx]
#     norm_values = 10**log_norm_values[idx].flatten()
#
#     # Now do a sampling from the minimum to the maximum with 1000 steps (to make sure that we cover also far away
#     # from the minimum)
#     grid = np.logspace(np.log10(min_value), np.log10(max_value), 1000)[1:-1]
#
#     grid_log_like_values = map(lnprob, np.log10(grid))
#
#     log_like_values = np.append(log_like_values, grid_log_like_values)
#     norm_values = np.append(norm_values, grid)
#
#     print("Minimum of -log(like): %s" % (log_like_values.max() * -1))
#     print("Maximum of -log(like): %s" % (log_like_values.min() * -1))
#
#     # Now sort by value
#     idx = np.argsort(norm_values)
#     log_like_values = log_like_values[idx]
#     norm_values = norm_values[idx]
#
#     # Make unique, in the very unlikely scenario where the emcee run and the regular grid provided
#     # some duplicated values
#     norm_values_u, idx = np.unique(norm_values, return_index=True)
#     log_like_values_u = log_like_values[idx]
#
#     # Restore parameter to its original state
#
#     like['GRB'].src.spectrum().parameter(norm.attrib['name']).setFree(True)
#     like['GRB'].src.spectrum().parameter(norm.attrib['name']).setBounds(0, 1e9)
#     like['GRB'].src.spectrum().parameter(norm.attrib['name']).setScale(old_scale)
#     like['GRB'].src.spectrum().parameter(norm.attrib['name']).setValue(old_value)
#     like['GRB'].src.spectrum().parameter(norm.attrib['name']).setBounds(float(old_min), float(old_max))
#
#
#     return -log_like_values_u, norm_values_u

# Main code
if __name__ == "__main__":

    configuration = Configuration()
    irfs = list(IRFS.IRFS.keys())
    irfs.append('auto')
    parser = argparse.ArgumentParser()

    parser.add_argument("triggername", help="Trigger name (in Fermi format: YYMMDDXXX)", type=str)
    parser.add_argument("--outfile", help="File for the results (will be overwritten)", type=str, required=True)
    parser.add_argument("--ra", help="R.A. of the object (J2000)", type=float)
    parser.add_argument("--dec", help="Dec. of the object (J2000)", type=float)
    parser.add_argument("--roi", help="Radius of the Region Of Interest (ROI)", type=float, required=True)
    parser.add_argument("--tstarts", help="Comma-separated list of start times (with respect to trigger)", type=str,
                        required=False, default=None)
    parser.add_argument("--tstops", help="Comma-separated list of stop times (with respect to trigger)", type=str,
                        required=False, default=None)
    parser.add_argument("--log_bins", help="Use logarithmically-spaced bins. Specify [tmin] [tmax] [n]. "
                                           "For example, '--log_bins 1.0 10000.0 30' uses 30 logarithmically spaced "
                                           "bins between 1.0 and 10k seconds. You should either use this, or --tstarts "
                                           "and --tstops", type=float, nargs=3,
                        required=False, default=None)
    parser.add_argument("--bin_file", help="A text file readable by numpy and the columns to read. For example, "
                                           "'--bin_file res.txt start end' will get the start and stop times from the "
                                           "columns 'start' and 'end' in the file res.txt.", type=str, nargs=3,
                        required=False, default=None)
    parser.add_argument("--zmax", help="Zenith cut", type=float, default=100.0)
    parser.add_argument("--emin", help="Minimum energy for the analysis", type=float, default=100.0)
    parser.add_argument("--emax", help="Maximum energy for the analysis", type=float, default=100000.0)
    parser.add_argument("--irf", help="Instrument Function to be used (IRF)", type=str, choices=irfs, required=True)
    parser.add_argument("--galactic_model", help="Galactic model for the likelihood", type=str, required=True,
                        choices=['template (fixed norm.)', 'template', 'none'])
    parser.add_argument("--particle_model", help="Particle model", type=str, required=True,
                        choices=['auto', 'isotr with pow spectrum', 'isotr template', 'none', 'bkge'])
    parser.add_argument("--tsmin", help="Minimum TS to consider a detection", type=float, default=20)
    parser.add_argument("--strategy", help="Strategy for Zenith cut: events or time", type=str,
                        choices=['time', 'events'],
                        default='time')
    parser.add_argument("--thetamax", help="Theta cut", type=float, default=180)
    parser.add_argument("--spectralfiles", help="Produce spectral files to be used in XSPEC?", type=str,
                        choices=['yes', 'no'], default='no')
    parser.add_argument("--liketype", help="Likelihood type (binned or unbinned)", type=str, default="unbinned",
                        choices=['binned', 'unbinned'])
    parser.add_argument("--optimizeposition", help="Optimize position with gtfindsrc?", type=str, default="no",
                        choices=['yes', 'no'])
    parser.add_argument("--datarepository", help="Directory where data are stored",
                        default=configuration.get('dataRepository'))
    parser.add_argument("--ltcube", help="Pre-computed livetime cube", default='', type=str)
    parser.add_argument("--expomap", help="pre-computed exposure map", default='', type=str)
    parser.add_argument('--ulphindex', help="Photon index for upper limits", default=-2, type=float)
    parser.add_argument('--flemin', help="Lower bound energy for flux/upper limit computation", default=None)
    parser.add_argument('--flemax', help="Upper bound energy for flux/upper limit computation", default=None)
    parser.add_argument('--fgl_mode',
                        help="Set 'complete' to use all FGL sources, set 'fast' to use only bright sources",
                        default='fast')
    parser.add_argument("--tsmap_spec",
                        help="A TS map specification of the type half_size,n_side. For example: "
                             "'--tsmap_spec 0.5,8' makes a TS map 1 deg x 1 deg with 64 points",
                        default=None)
    parser.add_argument("--filter_GTI", help="Automatically divide time intervals crossing GTIs", default=False,
                        action='store_true')
    parser.add_argument("--likelihood_profile",
                        help="Produce a text file containing the profile of the likelihood for a "
                             "changing normalization", default=False,
                        action='store_true')
    parser.add_argument("--remove_fits_files",
                        help="Whether to remove the FITS files of every interval in order to save disk space",
                        default=False, action='store_true')

    args = parser.parse_args()

    if args.ltcube != '':
        args.ltcube = os.path.abspath(os.path.expanduser(os.path.expandvars(args.ltcube)))

        assert os.path.exists(args.ltcube), "Livetime cube %s does not exist" % args.ltcube

    if args.expomap != '':
        args.expomap = os.path.abspath(os.path.expanduser(os.path.expandvars(args.expomap)))

        assert os.path.exists(args.expomap), "Exposure map %s does not exist" % args.expomap

    # Parse time interval specification

    if args.log_bins is None and args.bin_file is None:

        assert args.tstarts is not None and args.tstops is not None, "You have to use either --tstarts and --tstops or " \
                                                                     "--log_bins"

    else:

        assert args.tstarts is None and args.tstops is None, "You have to use either --tstarts and --tstops or " \
                                                             "--log_bins"

        if args.log_bins is not None:

            assert len(args.log_bins) == 3

            t1, t2, n_bins = args.log_bins

            edges = np.logspace(np.log10(float(t1)), np.log10(float(t2)), int(n_bins) + 1)

            args.tstarts = ",".join(map(str, edges[:-1]))
            args.tstops = ",".join(map(str, edges[1:]))

        else:

            assert args.bin_file is not None

            assert len(args.bin_file) == 3, "wrong syntax for --bin_file. Use --help for help."

            file_name, col1, col2 = args.bin_file

            d = np.recfromtxt(os.path.expandvars(os.path.expanduser(file_name)), names=True)

            args.tstarts = ",".join(map(str, d[col1]))
            args.tstops = ",".join(map(str, d[col2]))


    if args.likelihood_profile:

        # Check that we can import emcee
        try:

            import emcee

        except ImportError:

            raise ImportError("In order to make the likelihood profile you have to install emcee.")

    # Determine time intervals
    tstarts = numpy.array([float(x.replace('\\', "")) for x in args.tstarts.split(",")])
    tstops = numpy.array([float(x.replace('\\', "")) for x in args.tstops.split(",")])
    print("\nRequested intervals:")
    print("------------------------------------------------------")
    for t1, t2 in zip(tstarts, tstops):
        print(("%-20s - %s" % (t1, t2)))

    # Check if data exists, otherwise download them
    try:
        dataset = dataHandling.getLATdataFromDirectory(
            os.path.join(args.datarepository, 'bn%s' % args.triggername.replace('bn', '')))
    except:
        raise

    if (dataset is None):
        # Download data
        print(("\nData for trigger %s are not available. Let's download them!" % (args.triggername)))
        cmdLine = "gtdownloadLATdata.py triggername=%s timebefore=%s timeafter=%s datarepository=%s" % (
        args.triggername,
        min(tstarts.min(), -5000),
        max(tstarts.max(), 10000),
        args.datarepository)
        subprocess.call(cmdLine, shell=True)
        dataset = dataHandling.getLATdataFromDirectory(
            os.path.join(args.datarepository, 'bn%s' % args.triggername.replace('bn', '')))
    pass

    print("\nData files:")
    print("-----------")
    for k, v in dataset.items():
        print(("%-20s %s" % (k, v)))
    pass

    # Now get R.A. and Dec. if not specified
    if (args.ra is None or args.dec is None):
        header = pyfits.getheader(dataset['eventfile'], 'EVENTS')
        ra, dec = (header['RA_OBJ'], header['DEC_OBJ'])
        args.ra = ra
        args.dec = dec
    else:
        ra = args.ra
        dec = args.dec
    pass

    print("\nROI:")
    print("-----")
    print(("%-20s %s" % ('R.A.', ra)))
    print(("%-20s %s" % ('Dec.', dec)))
    print(("%-20s %s" % ('Radius', args.roi)))

    if args.filter_GTI:

        # Fix the requested time intervals according to the GTIs

        # First select between the first time and the last time, using the Zenith cut and strategy
        # I am going to use
        # Select data
        t1 = min(tstarts)
        t2 = max(tstops)

        targs = {}
        targs['rad'] = args.roi
        targs['eventfile'] = dataset['eventfile']
        targs['zmax'] = args.zmax
        targs['thetamax'] = args.thetamax
        targs['emin'] = args.emin
        targs['emax'] = args.emax
        targs['skymap'] = '%s_LAT_skymap_%s-%s.fit' % (args.triggername, t1, t2)
        targs['rspfile'] = dataset['rspfile']
        targs['strategy'] = args.strategy
        targs['ft2file'] = dataset['ft2file']
        targs['tstart'] = t1
        targs['tstop'] = t2
        targs['ra'] = args.ra
        targs['dec'] = args.dec
        # For this initial selection use the less restrictive IRF
        targs['irf'] = 'p8_transient020e'
        targs['allowEmpty'] = 'no'

        printCommand("gtdocountsmap.py", targs)
        try:
            _, skymap, _, filteredeventfile, _, _, _, _ = gtdocountsmap.run(**targs)
        except:
            raise

        gtis = pyfits.getdata(filteredeventfile, 'GTI')
        trigger_time = dataHandling.getTriggerTime(filteredeventfile)

        tstarts, tstops = clean_intervals(tstarts, tstops,
                                          gtis.field("START") - trigger_time,
                                          gtis.field("STOP") - trigger_time)

        print("\nAfter intersecting with GTI these are the intervals:")
        print("------------------------------------------------------")
        for t1, t2 in zip(tstarts, tstops):
            print(("%-20s - %s" % (t1, t2)))

    results = []
    initialWorkdir = os.getcwd()

    for i, t1, t2 in zip(list(range(1, len(tstarts) + 1)), tstarts, tstops):

        print(("\nInterval # %s (%s-%s):" % (i, t1, t2)))
        print("-----------------------\n")

        # Create a work dir and move there
        dirname = os.path.abspath("interval%s-%s" % (t1, t2))

        if os.path.exists(dirname):
            print(("%s already exists, skipping" % dirname))

            continue

        particle_model = args.particle_model

        if args.irf == 'auto':

            if t2 - t1 <= 100.0:

                args.irf = 'p8_transient020e'

            elif 100.0 < t2 - t1 < 1000.0:

                args.irf = 'p8_transient010e'

            else:

                args.irf = 'p8_source'

            print(("\nAutoselected %s class\n" % args.irf))

        irf = args.irf

        try:
            os.makedirs(dirname)
        except:
            pass
        pass

        try:
            os.chdir(dirname)
        except:
            raise RuntimeError("Could not create/access directory %s" % (dirname))
        pass

        # Select data
        targs = {}
        targs['rad'] = args.roi
        targs['eventfile'] = dataset['eventfile']
        targs['zmax'] = args.zmax
        targs['thetamax'] = args.thetamax
        targs['emin'] = args.emin
        targs['emax'] = args.emax
        targs['skymap'] = '%s_LAT_skymap_%s-%s.fit' % (args.triggername, t1, t2)
        targs['rspfile'] = dataset['rspfile']
        targs['strategy'] = args.strategy
        targs['ft2file'] = dataset['ft2file']
        targs['tstart'] = t1
        targs['tstop'] = t2
        targs['ra'] = args.ra
        targs['dec'] = args.dec
        targs['irf'] = irf
        targs['allowEmpty'] = 'no'

        printCommand("gtdocountsmap.py", targs)
        try:
            _, skymap, _, filteredeventfile, _, _, _, _ = gtdocountsmap.run(**targs)
        except:
            print("\nERROR: could not complete selection of data for this interval.")
            continue

        # Build XML file
        targs = {}
        targs['xmlmodel'] = '%s_LAT_xmlmodel_%s-%s.xml' % (args.triggername, t1, t2)
        targs['filteredeventfile'] = filteredeventfile
        targs['galactic_model'] = args.galactic_model
        targs['particle_model'] = particle_model
        targs['ra'] = args.ra
        targs['dec'] = args.dec
        targs['fgl_mode'] = args.fgl_mode
        targs['ft2file'] = dataset['ft2file']
        targs['source_model'] = 'powerlaw2'
        printCommand("gtbuildxmlmodel", targs)
        _, xmlmodel = gtbuildxmlmodel.run(**targs)

        # Now if the user has specified a specific photon index for upper limits,
        # change the photon index in the XML file

        # Save parameters in comments (ET will strip them out)

        pars_in_comments = {}

        for key in ['OBJECT', 'RA', 'DEC', 'IRF']:
            pars_in_comments[key] = dataHandling._getParamFromXML(xmlmodel, key)

        # Now change the photon index in the XML file

        tree = ET.parse(xmlmodel)
        root = tree.getroot()
        index = root.findall("./source[@name='%s']/spectrum/parameter[@name='Index']" % 'GRB')[0]

        if args.ulphindex == -1.0:
            args.ulphindex += 0.01

        index.set('value', str(args.ulphindex))

        tree.write(xmlmodel)

        # Add the parameters in comments back

        dataHandling._writeParamIntoXML(xmlmodel, **pars_in_comments)

        targs = {}
        targs['spectralfiles'] = args.spectralfiles
        targs['xmlmodel'] = xmlmodel
        targs['liketype'] = args.liketype
        targs['filteredeventfile'] = filteredeventfile
        targs['rspfile'] = dataset['rspfile']
        targs['showmodelimage'] = 'no'
        targs['tsmin'] = args.tsmin
        targs['optimizeposition'] = 'no'
        targs['ft2file'] = dataset['ft2file']
        targs['skymap'] = skymap
        targs['flemin'] = args.flemin
        targs['flemax'] = args.flemax

        if args.ltcube != '':

            if not os.path.exists(args.ltcube):
                raise IOError("Livetime cube %s does not exists!" % (args.ltcube))

            targs['ltcube'] = args.ltcube

        if args.expomap != '':

            if not os.path.exists(args.expomap):
                raise IOError("Exposure map %s does not exists!" % (args.expomap))

            targs['expomap'] = args.expomap

        printCommand("gtdolike.py", targs)
        (_, outfilelike, _, grb_TS,
         _, bestra, _, bestdec,
         _, poserr, _, distance,
         _, sources) = gtdolike.run(**targs)

        # Get root of the name
        root_name = os.path.splitext(os.path.basename(filteredeventfile))[0]

        # Find ltcube

        ltcubes = glob.glob("%s_ltcube.fit*" % root_name)

        assert len(ltcubes) == 1, "Couldn't find ltcube"

        ltcube = ltcubes[0]

        # Find expomap

        expmaps = glob.glob("%s_expomap.fit*" % root_name)

        assert len(expmaps) == 1, "Couldn't find exopmap"

        expmap = expmaps[0]

        # Find XML model output of gtdolike
        xmls = glob.glob("%s_likeRes.xml" % root_name)

        assert len(xmls) == 1, "Couldn't find XML"

        xml_res = xmls[0]

        # If the TS map is required, let's do it

        if args.tsmap_spec is not None:
            half_size, n_side = args.tsmap_spec.replace(" ", "").split(",")

            obs = UnbinnedAnalysis.UnbinnedObs(filteredeventfile, dataset['ft2file'], expMap=expmap, expCube=ltcube)
            like = UnbinnedAnalysis.UnbinnedAnalysis(obs, xml_res, 'MINUIT')

            ftm = FastTSMap(like)
            (bestra, bestdec), maxTS = ftm.search_for_maximum(args.ra, args.dec, float(half_size), int(n_side),
                                                              verbose=False)

        # Now append the results for this interval
        grb = [x for x in sources if x.name.find("GRB") >= 0][0]

        if args.tsmap_spec is not None:

            if maxTS > grb.TS:
                print("\n\n=========================================")
                print(" Fast TS Map has found a better position")
                print("=========================================\n\n")

                # grb.ra                   = float(bestra)
                # grb.dec                  = float(bestdec)
                grb.TS = float(maxTS)

                print(("(R.A., Dec.) = (%.3f, %3f) with TS = %.2f\n" % (bestra, bestdec, grb.TS)))

        else:

            # Do nothing, so that grb.ra and grb.dec will stay what they are already
            pass

        if args.likelihood_profile:
            # Make a profile of the likelihood with a variable normalization
            obs = UnbinnedAnalysis.UnbinnedObs(filteredeventfile, dataset['ft2file'], expMap=expmap, expCube=ltcube)
            like = UnbinnedAnalysis.UnbinnedAnalysis(obs, xml_res, 'MINUIT')

            # Setup likelihood object like is done in gtburst
            like, phIndex_beforeFit, grb_name = dataHandling.LATData.setup_likelihood_object(like, xmlmodel)

            lpw = LikelihoodProfiler(like, xml_res)

            if grb.TS < args.tsmin:

                print("\nForcing index to -2.0 and sampling\n")

                lpw.get_likelihood_profile(forced_photon_index=-2.0)

            else:

                print("\nSampling\n")
                lpw.get_likelihood_profile()

            lpw.save("log_like_profile")
            fig = lpw.plot()

            fig.savefig("log_like_profile.png")

            plt.close(fig)

        if args.remove_fits_files:

            fits_files = glob.glob("*.fit*")

            for ff in fits_files:

                os.remove(ff)

        grb.name = args.triggername
        grb.tstart = t1
        grb.tstop = t2
        grb.roi = args.roi
        grb.irf = irf
        grb.zmax = args.zmax
        grb.thetamax = args.thetamax
        grb.strategy = args.strategy
        results.append(grb)

        os.chdir(initialWorkdir)
    pass

    try:

        writeSourceListToFile(results, args.outfile)

    except IOError:

        print("Looks like the analysis has failed. No output file produced!")
pass
