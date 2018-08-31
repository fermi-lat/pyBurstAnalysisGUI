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


# Main code
if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--tstart", help="Use only time intervals with tstart after this time", type=float,
                        required=True, default=None)
    parser.add_argument("--tstop", help="Use only time intervals with tstop before this time", type=float,
                        required=True, default=None)
    parser.add_argument("--intervals_dir", help="Directory containing the intervals produced by doTimeResolvedLike. "
                                                "Default: current directory",
                        type=str, default=os.getcwd(), required=False)

    parser.add_argument("--outroot", help="Root name to use for the output files",
                        type=str, required=True)

    args = parser.parse_args()

    intervals_dir = os.path.abspath(os.path.expandvars(os.path.expanduser(args.intervals_dir)))

    intervals = glob.glob(os.path.join(intervals_dir, "interval*"))

    print("Found %i intervals in directory %s" % (len(intervals), intervals_dir))

    if len(intervals) == 0:

        print("No interval found, nothing to do.")

        sys.exit(0)

    # Load data
    containers = []

    for interval in intervals:

        start, stop = map(float, os.path.basename(interval).replace("interval", "").split("-"))

        log_like_profile_file = os.path.join(interval, "log_like_profile.npz")

        if not os.path.exists(log_like_profile_file):

            raise IOError("Interval %s does not contain the log_like_profile.npz file" % interval)

        d = np.load(log_like_profile_file)

        containers.append(IntervalContainer(start, stop, d['norm_values'], -d['mlog_like_values'], int(stop - start)))

    # Instance container

    cl = CastroLike("castro", containers)
    n_active_intervals = cl.set_active_measurements(args.tstart, args.tstop)

    assert n_active_intervals > 2

    # Set up power law model

    po = Powerlaw(K=1e-4, piv=1.0)
    po.K.bounds = (1e-12, None)

    model0 = Model(PointSource("GRB", 0.0, 0.0, po))

    # Fit power law model
    jl0 = JointLikelihood(model0, DataList(cl))

    jl0.set_minimizer("ROOT")

    po_best_fit, po_log_like = jl0.fit()

    po_best_fit.to_csv("%s_po.csv" % args.outroot)

    fig = cl.plot()

    _ = fig.axes[0].set_ylabel("Flux (erg/cm2/s)")
    _ = fig.axes[0].set_xlabel("Time")

    fig.savefig("%s_pow.png" % args.outroot)

    # Set up broken power law model

    bknpo = Broken_powerlaw(piv=1.0)

    bknpo.K.bounds = (1e-12, None)
    bknpo.xb.bounds = (cl.start, cl.stop)
    bknpo.alpha.bounds = (-10, 0.0)
    bknpo.beta.bounds = (-10, 0.0)

    bknpo.K = 5e-2
    bknpo.xb = (cl.stop + cl.start)/2.0
    bknpo.alpha = -2.0
    bknpo.beta = -1.0

    model1 = Model(PointSource("GRB", 0.0, 0.0, bknpo))

    # Fit
    jl1 = JointLikelihood(model1, DataList(cl))

    # Create an instance of the GRID minimizer
    grid_minimizer = GlobalMinimization("grid")

    # Create an instance of a local minimizer, which will be used by GRID
    local_minimizer = LocalMinimization("ROOT")

    # Define a grid for mu as 10 steps between 1 and 80
    my_grid = {model1.GRB.spectrum.main.shape.xb: np.logspace(np.log10(cl.start), np.log10(cl.stop), 20)}

    # Setup the global minimization
    # NOTE: the "callbacks" option is useless in a normal 3ML analysis, it is
    # here only to keep track of the evolution for the plot
    grid_minimizer.setup(second_minimization=local_minimizer, grid=my_grid)

    # Set the minimizer for the JointLikelihood object
    jl1.set_minimizer(grid_minimizer)

    jl1.set_minimizer("ROOT")

    bknpo_best_fit, bknpo_log_like = jl1.fit()

    bknpo_best_fit.to_csv("%s_bknpo.csv" % args.outroot)

    fig = cl.plot()

    _ = fig.axes[0].set_ylabel("Flux (erg/cm2/s)")
    _ = fig.axes[0].set_xlabel("Time")

    fig.savefig("%s_bknpo.png" % args.outroot)

    TS = 2 * (po_log_like.loc['castro'] - bknpo_log_like.loc['castro']).values[0]

    print("TS is %s" % TS)
