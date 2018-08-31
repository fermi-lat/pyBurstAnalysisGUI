import numpy as np
import emcee
import os
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt


class LikelihoodProfiler(object):

    def __init__(self, likelihood_object, best_fit_xml):

        self._like = likelihood_object
        self._best_fit_xml = os.path.expandvars(os.path.expanduser(best_fit_xml))

        assert os.path.exists(self._best_fit_xml), "XML file %s does not exist" % self._best_fit_xml

        # This will be filled by the get_likelihood_profile method
        self._mlog_like_values = None
        self._norm_values = None

    @staticmethod
    def _photon_flux_to_energy_flux_conv_factor(emin_MeV, emax_MeV, photon_index):

        MeVtoErg = 1.60217646E-6

        emin_erg = emin_MeV * MeVtoErg
        emax_erg = emax_MeV * MeVtoErg

        assert photon_index != -1.0

        if photon_index != -2.0:

            conv = (1. + photon_index) / (2.0 + photon_index) * (
                    pow(emax_erg, photon_index + 2) - pow(emin_erg, photon_index + 2)) / (
                           pow(emax_erg, photon_index + 1) - pow(emin_erg, photon_index + 1))
        else:

            conv = (emin_erg) * (emax_erg) / (emax_erg - emin_erg) * np.log(emax_erg / emin_erg)

        return conv

    def get_likelihood_profile(self, source_name='GRB', energy_flux=False, forced_photon_index=None):

        # Find name of normalization parameter
        tree = ET.parse(self._best_fit_xml)
        root = tree.getroot()
        norm = root.findall("./source[@name='%s']/spectrum/parameter[@name='Integral']" % source_name)[0]

        # Get the current settings of the parameter so we can restore it later
        old_scale = float(norm.attrib['scale'])
        old_value = self._like[source_name].src.spectrum().parameter(norm.attrib['name']).getValue()
        old_min = float(norm.attrib['min'])
        old_max = float(norm.attrib['max'])

        # Make sure we are at the best fit

        _ = self._like.fit(verbosity=0, covar=False)

        # Find best fit value
        best_fit_norm = self._like[source_name].src.spectrum().parameter(norm.attrib['name']).getValue() * old_scale

        if forced_photon_index is None:

            best_fit_index = float(self._like[source_name].src.spectrum().parameter('Index').getValue())

        else:

            best_fit_index = forced_photon_index

        flux_emin_MeV = float(self._like[source_name].src.spectrum().parameter('LowerLimit').getValue())
        flux_emax_MeV = float(self._like[source_name].src.spectrum().parameter('UpperLimit').getValue())

        # Temporarily set the scale to 1.0 so we can modify boundaries and so on without incurring in errors
        self._like[source_name].src.spectrum().parameter(norm.attrib['name']).setScale(1.0)

        # Use the minimum and maximum in the file
        min_value = old_min * old_scale
        max_value = old_max * old_scale

        # Get the logarithm because we will sample in the log space
        min_log_value = np.log10(min_value)
        max_log_value = np.log10(max_value)

        # Store best fit parameters
        best_fit_values = []

        for parameter in self._like.model.params:

            best_fit_values.append(parameter.getValue())

        #####################################
        # Strategy:
        # We apply an adaptive sampling strategy for the profiled likelihood so that we will have many points close to
        # the maximum, where we need very high resolution, and less points farther away.
        # his will allow the interpolation of the profile to have very high fidelity close to the maximum.
        # In order to achieve this, I exploit the capability
        # of Markov Chain Monte Carlo methods, in particular emcee, that will return samples for the normalization in
        # proportion to the likelihood value: many samples close to the maximum, fewer sample away from it.
        # Then, since we need to explore the entire range between minimum and maximum, I will also make a regular grid
        # going from minimum to maximum and append it to the grid returned by emcee
        ####################################

        # Function that will be used to "sample" the likelihood

        like = self._like

        def lnprob(log_norm):

            if log_norm < min_log_value or log_norm > max_log_value:

                # Outside of boundaries

                return -np.inf

            else:

                # Restore best fit parameters first
                for i, parameter in enumerate(like.model.params):

                    if parameter.getName() == norm.attrib['name']:

                        # This will be set below
                        continue

                    parameter.setValue(best_fit_values[i])

                like[source_name].src.spectrum().parameter(norm.attrib['name']).setFree(True)

                like[source_name].src.spectrum().parameter(norm.attrib['name']).setValue(float(10**log_norm))

                like[source_name].src.spectrum().parameter(norm.attrib['name']).setFree(False)

                like[source_name].src.spectrum().parameter('Index').setValue(best_fit_index)
                like[source_name].src.spectrum().parameter('Index').setFree(True)

                like.syncSrcParams()

                log_like = -like.fit(verbosity=0, covar=False)

                like[source_name].src.spectrum().parameter(norm.attrib['name']).setFree(True)
                like[source_name].src.spectrum().parameter('Index').setFree(True)

                return log_like

        # Make sure we start from the best fit, and with the appropriate boundaries

        # Temporarily make the boundaries very large, so the next assignment will not fail
        self._like[source_name].src.spectrum().parameter(norm.attrib['name']).setBounds(0.0,
                                                                            1e9)

        self._like[source_name].src.spectrum().parameter(norm.attrib['name']).setValue(best_fit_norm)

        # Restore the boundaries
        self._like[source_name].src.spectrum().parameter(norm.attrib['name']).setBounds(min_value,
                                                                            max_value)

        n_walkers = 30

        # Start all the walkers close to the best fit value
        p0 = [[best_fit_norm + np.random.uniform(-best_fit_norm / 50.0, best_fit_norm / 50.0)]
              for i in range(n_walkers)]

        # Sample
        sampler = emcee.EnsembleSampler(n_walkers, 1, lnprob)
        _ = sampler.run_mcmc(np.log10(p0), 100)

        log_norm_values = sampler.flatchain
        log_like_values = sampler.flatlnprobability

        # Remove infinite values
        idx = np.isfinite(log_like_values)
        log_like_values = log_like_values[idx]
        photon_flux = 10**log_norm_values[idx].flatten()

        # Now do a sampling from the minimum to the maximum with 1000 steps (to make sure that we cover also far away
        # from the minimum)
        grid = np.logspace(np.log10(min_value), np.log10(max_value), 1000)[1:-1]

        grid_log_like_values = map(lnprob, np.log10(grid))

        log_like_values = np.append(log_like_values, grid_log_like_values)
        photon_flux = np.append(photon_flux, grid)

        # Now sort by value
        idx = np.argsort(photon_flux)
        log_like_values = log_like_values[idx]
        photon_flux = photon_flux[idx]

        # Make unique, in the very unlikely scenario where the emcee run and the regular grid provided
        # some duplicated values
        photon_flux_u, idx = np.unique(photon_flux, return_index=True)
        log_like_values_u = log_like_values[idx]

        # Restore parameter to its original state

        self._like[source_name].src.spectrum().parameter(norm.attrib['name']).setFree(True)
        self._like[source_name].src.spectrum().parameter(norm.attrib['name']).setBounds(0, 1e9)
        self._like[source_name].src.spectrum().parameter(norm.attrib['name']).setScale(old_scale)
        self._like[source_name].src.spectrum().parameter(norm.attrib['name']).setValue(old_value)
        self._like[source_name].src.spectrum().parameter(norm.attrib['name']).setBounds(float(old_min), float(old_max))

        conv = self._photon_flux_to_energy_flux_conv_factor(
            flux_emin_MeV, flux_emax_MeV,
            best_fit_index)

        if energy_flux:

            # Convert the photon flux to the energy flux
            norm_values = photon_flux_u * conv

        else:

            norm_values = photon_flux_u

        self._mlog_like_values = -log_like_values_u
        self._norm_values = norm_values
        self._conv_factor = conv
        self._photon_index = best_fit_index

        return self._mlog_like_values, self._norm_values

    def save(self, outfile):

        assert self._mlog_like_values is not None and self._norm_values is not None, \
            "You need to run get_likelihood_profile first"

        # Save to file
        np.savez(outfile,
                 mlog_like_values=self._mlog_like_values,
                 norm_values=self._norm_values,
                 conversion_factor=self._conv_factor,
                 photon_index=self._photon_index)

    def plot(self):

        fig = plt.figure()

        _ = plt.loglog(self._norm_values, self._mlog_like_values - self._mlog_like_values.min() + 1.0)

        plt.ylabel("-log(likelihood) - %.3g" % (self._mlog_like_values.min() + 1.0))
        plt.xlabel("Flux")

        return fig
