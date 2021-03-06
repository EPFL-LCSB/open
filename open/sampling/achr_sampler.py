# -*- coding: utf-8 -*-

"""Provide ACHR sampler."""

from __future__ import absolute_import, division

import numpy as np
import pandas

#from cobra.sampling.hr_sampler import HRSampler, step
from open.sampling.samplers import HRSampler, step

class ACHRSampler(HRSampler):
    """Artificial Centering Hit-and-Run sampler.

    A sampler with low memory footprint and good convergence.

    Parameters
    ----------
    model : cobra.Model
        The cobra model from which to generate samples.
    thinning : int, optional
        The thinning factor of the generated sampling chain. A thinning of 10
        means samples are returned every 10 steps.
    nproj : int > 0, optional
        How often to reproject the sampling point into the feasibility space.
        Avoids numerical issues at the cost of lower sampling. If you observe
        many equality constraint violations with `sampler.validate` you should
        lower this number.
    seed : int > 0, optional
        Sets the random number seed. Initialized to the current time stamp if
        None.

    Attributes
    ----------
    model : cobra.Model
        The cobra model from which the samples get generated.
    thinning : int
        The currently used thinning factor.
    n_samples : int
        The total number of samples that have been generated by this
        sampler instance.
    problem : collections.namedtuple
        A python object whose attributes define the entire sampling problem in
        matrix form. See docstring of `Problem`.
    warmup : numpy.matrix
        A matrix of with as many columns as reactions in the model and more
        than 3 rows containing a warmup sample in each row. None if no warmup
        points have been generated yet.
    retries : int
        The overall of sampling retries the sampler has observed. Larger
        values indicate numerical instabilities.
    seed : int > 0, optional
        Sets the random number seed. Initialized to the current time stamp if
        None.
    nproj : int
        How often to reproject the sampling point into the feasibility space.
    fwd_idx : numpy.array
        Has one entry for each reaction in the model containing the index of
        the respective forward variable.
    rev_idx : numpy.array
        Has one entry for each reaction in the model containing the index of
        the respective reverse variable.
    prev : numpy.array
        The current/last flux sample generated.
    center : numpy.array
        The center of the sampling space as estimated by the mean of all
        previously generated samples.

    Notes
    -----
    ACHR generates samples by choosing new directions from the sampling space's
    center and the warmup points. The implementation used here is the same
    as in the Matlab Cobra Toolbox [2]_ and uses only the initial warmup points
    to generate new directions and not any other previous iterates. This
    usually gives better mixing since the startup points are chosen to span
    the space in a wide manner. This also makes the generated sampling chain
    quasi-markovian since the center converges rapidly.

    Memory usage is roughly in the order of (2 * number reactions)^2
    due to the required nullspace matrices and warmup points. So large
    models easily take up a few GB of RAM.

    References
    ----------
    .. [1] Direction Choice for Accelerated Convergence in Hit-and-Run Sampling
       David E. Kaufman Robert L. Smith
       Operations Research 199846:1 , 84-95
       https://doi.org/10.1287/opre.46.1.84
    .. [2] https://github.com/opencobra/cobratoolbox

    """

    def __init__(self, model, thinning=100, nproj=None, seed=None):
        """Initialize a new ACHRSampler."""

        super(ACHRSampler, self).__init__(model, thinning, nproj=nproj,
                                          seed=seed)
        self.generate_fva_warmup()
        self.prev = self.center = self.warmup.mean(axis=0)
        np.random.seed(self._seed)

    def __single_iteration(self):
        pi = np.random.randint(self.n_warmup)

        # mix in the original warmup points to not get stuck
        delta = self.warmup[pi, ] - self.center
        self.prev = step(self, self.prev, delta)

        if self.problem.homogeneous and (self.n_samples *
                                         self.thinning % self.nproj == 0):
            self.prev = self._reproject(self.prev)
            self.center = self._reproject(self.center)

        self.center = ((self.n_samples * self.center) / (self.n_samples + 1) +
                       self.prev / (self.n_samples + 1))
        self.n_samples += 1

    def sample(self, n, fluxes=True):
        """Generate a set of samples.

        This is the basic sampling function for all hit-and-run samplers.

        Parameters
        ----------
        n : int
            The number of samples that are generated at once.
        fluxes : boolean
            Whether to return fluxes or the internal solver variables. If set
            to False will return a variable for each forward and backward flux
            as well as all additional variables you might have defined in the
            model.

        Returns
        -------
        numpy.matrix
            Returns a matrix with `n` rows, each containing a flux sample.

        Notes
        -----
        Performance of this function linearly depends on the number
        of reactions in your model and the thinning factor.

        """

        samples = np.zeros((n, self.warmup.shape[1]))

        for i in range(1, self.thinning * n + 1):
            self.__single_iteration()

            if i % self.thinning == 0:
                samples[i//self.thinning - 1, ] = self.prev

        if fluxes:
            names = [r.id for r in self.model.reactions]

            return pandas.DataFrame(
                samples[:, self.fwd_idx] - samples[:, self.rev_idx],
                columns=names)
        else:
            names = [v.name for v in self.model.variables]

            return pandas.DataFrame(samples, columns=names)
