#!/usr/bin/env python

import argparse

try:

    from threeML import *

except ImportError:

    raise ImportError("3ML is required to run this script")

from threeML.plugins.experimental.CastroLike import CastroLike, IntervalContainer
import glob
import os
import sys


def find_profiles(intervals_dir):

    intervals_dir = os.path.abspath(os.path.expandvars(os.path.expanduser(intervals_dir)))

    intervals = glob.glob(os.path.join(intervals_dir, "interval*"))

    print("Found %i intervals in directory %s" % (len(intervals), intervals_dir))

    if len(intervals) == 0:
        print("No interval found, nothing to do.")

        return {}

    # Find files
    profiles_dict = {}

    for interval in intervals:

        start, stop = map(float, os.path.basename(interval).replace("interval", "").split("-"))

        log_like_profile_file = os.path.join(interval, "log_like_profile.npz")

        if not os.path.exists(log_like_profile_file):
            raise IOError("Interval %s does not contain the log_like_profile.npz file" % interval)

        profiles_dict[(start, stop)] = log_like_profile_file

    return profiles_dict


def go(profiles_dict, tstart, tstop, outroot, local_minimizer="ROOT"):

    containers = []

    for interval in profiles_dict:

        start, stop = interval

        log_like_profile_file = profiles_dict[interval]

        d = np.load(log_like_profile_file)

        if stop - start < 1.0:

            n_integration_points = 3

        elif stop - start < 100.0:

            # One point every two seconds
            n_integration_points = (stop - start) / 2.0

        elif stop - start < 1000.0:

            n_integration_points = (stop - start) / 10.0

        else:

            n_integration_points = (stop - start) / 50.0

        n_integration_points = max(3, int(n_integration_points / 2.0) * 2 + 1)

        containers.append(IntervalContainer(start, stop, d['norm_values'], -d['mlog_like_values'],
                                            n_integration_points))

    # Instance container

    cl = CastroLike("castro", containers)
    n_active_intervals = cl.set_active_measurements(tstart, tstop)

    assert n_active_intervals > 2

    # Get first point
    first_start = cl.active_containers[0].start
    first_stop = cl.active_containers[0].stop
    tc = (first_stop + first_start) / 2.0
    first_flux = cl.active_containers[0].get_measurement()[1]

    # Set up power law model

    po = Powerlaw(K=first_flux, piv=tc)
    po.K.bounds = (first_flux / 1000.0, None)

    model0 = Model(PointSource("GRB", 0.0, 0.0, po))

    # Fit power law model
    jl0 = JointLikelihood(model0, DataList(cl))

    jl0.set_minimizer(local_minimizer)

    po_best_fit, po_log_like = jl0.fit()

    powerlaw_csv = "%s_po.csv" % outroot

    po_best_fit.to_csv(powerlaw_csv)

    fig = cl.plot()

    _ = fig.axes[0].set_ylabel("Flux (erg/cm2/s)")
    _ = fig.axes[0].set_xlabel("Time")

    fig.savefig("%s_pow.png" % outroot)

    # Set up broken power law model

    bknpo = Broken_powerlaw(K=first_flux, piv=tc)

    bknpo.K.bounds = (first_flux / 1000.0, None)
    bknpo.xb.bounds = (cl.start, cl.stop)
    bknpo.alpha.bounds = (-10, 10.0)
    bknpo.beta.bounds = (-10, 10.0)

    bknpo.xb = (cl.stop + cl.start) / 2.0
    bknpo.alpha = -2.0
    bknpo.beta = -1.0

    model1 = Model(PointSource("GRB", 0.0, 0.0, bknpo))

    # Fit
    jl1 = JointLikelihood(model1, DataList(cl))

    # Create an instance of the GRID minimizer
    grid_minimizer = GlobalMinimization("grid")

    # Create an instance of a local minimizer, which will be used by GRID
    local_minimizer = LocalMinimization(local_minimizer)

    # Define a grid for mu as 10 steps between 1 and 80
    my_grid = {model1.GRB.spectrum.main.shape.xb: np.logspace(np.log10(cl.start + 0.1), np.log10(cl.stop - 0.1), 5),
               model1.GRB.spectrum.main.shape.alpha: [-2.5, -2.0, -1.5],
               model1.GRB.spectrum.main.shape.beta: [-2.0, -1.0, -0.5],
               }

    # Setup the global minimization
    # NOTE: the "callbacks" option is useless in a normal 3ML analysis, it is
    # here only to keep track of the evolution for the plot
    grid_minimizer.setup(second_minimization=local_minimizer, grid=my_grid)

    # Set the minimizer for the JointLikelihood object
    jl1.set_minimizer(grid_minimizer)

    bknpo_best_fit, bknpo_log_like = jl1.fit()

    bknpo_csv = "%s_bknpo.csv" % outroot

    bknpo_best_fit.to_csv(bknpo_csv)

    fig = cl.plot()

    _ = fig.axes[0].set_ylabel("Flux (erg/cm2/s)")
    _ = fig.axes[0].set_xlabel("Time")

    fig.savefig("%s_bknpo.png" % outroot)

    TS = 2 * (po_log_like.loc['castro'] - bknpo_log_like.loc['castro']).values[0]

    print("TS is %s" % TS)

    return po_best_fit, bknpo_best_fit, TS


def main(arguments):

    profiles_dict = find_profiles(arguments.intervals_dir)

    go(profiles_dict, arguments.tstart, arguments.tstop, arguments.outroot, arguments.local_minimizer)


# Main code
if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--tstart", help="Use only time intervals with tstart after this time", type=float,
                        required=True, default=None)
    parser.add_argument("--tstop", help="Use only time intervals with tstop before this time", type=float,
                        required=True, default=None)
    parser.add_argument("--local_minimizer", help="Minimizer to use as local minimizer (default: ROOT)", type=str,
                        required=False, default="ROOT")
    parser.add_argument("--intervals_dir", help="Directory containing the intervals produced by doTimeResolvedLike. "
                                                "Default: current directory",
                        type=str, default=os.getcwd(), required=False)

    parser.add_argument("--outroot", help="Root name to use for the output files",
                        type=str, required=True)

    args = parser.parse_args()

    main(args)


