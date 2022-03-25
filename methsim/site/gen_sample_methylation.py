from typing import List, Union
import numpy as np
from methsim.site import simulate_site_matrix, generate_sample_methylation, generate_site_coefs
from methsim.sample import gen_pheno_matrix
from methsim.utilities import Scaler


class GenerateSampleMethylation:

    def __init__(self):
        self.methylation_sites = None
        self.sim_site_order = {}

    def generate_methylation_sites(self, deviation_low=0.0, deviation_high=0.1,
                                   delta_low=0.0, delta_high=0.5, number_of_sites=5000):
        # reset site order if generating methylation sites
        self.sim_site_order = {}
        self.methylation_sites = simulate_site_matrix(deviation_low, deviation_high,
                                                      delta_low, delta_high,
                                                      number_of_sites)

    def generate_sample_methylation(self, samples, sim_conditions: List[List[Union[List[str], int, np.ndarray]]],
                                    scale_min=None, scale_max=None):
        site_values, site_error = [], []
        site = 0
        for phenotypes, sim_sites, pheno_weights in sim_conditions:
            pheno_values = self.get_phenotype_values(samples, phenotypes, scale_min, scale_max)
            for count in range(sim_sites):
                if site in self.sim_site_order:
                    meth_site, meth_phenotypes, coefs = self.sim_site_order[site]
                    assert phenotypes == meth_phenotypes
                else:
                    meth_site = np.random.randint(0, self.methylation_sites.shape[0])
                    coefs = generate_site_coefs(self.methylation_sites[meth_site][0],
                                                self.methylation_sites[meth_site][1],
                                                state_weights=pheno_weights)
                    self.sim_site_order[site] = (meth_site, phenotypes, coefs)
                meth, error = generate_sample_methylation(pheno_values, self.methylation_sites[meth_site][0],
                                                          self.methylation_sites[meth_site][2], coefs)
                site_values.append(meth)
                site_error.append(error)
                site += 1
        return np.array(site_values), np.array(site_values)

    def get_phenotype_values(self, samples, phenotypes, scale_min=None, scale_max=None, verbose=False):
        if not isinstance(phenotypes, list):
            values = gen_pheno_matrix(samples, [phenotypes], key='value')
        else:
            values = gen_pheno_matrix(samples, phenotypes, key='value')
        unit_vectors = []
        for pheno in range(values.shape[1]):
            if np.isnan(np.sum(values[:, pheno])):
                if verbose:
                    print(f'{phenotypes} not set in all samples, returning unit vector')
                values[:, pheno] = np.ones((len(samples)))
                # save unit vector pos to avoid scaling
                unit_vectors.append(pheno)
        return self.scale_outputs(values, unit_vectors, scale_min, scale_max)

    @staticmethod
    def scale_outputs(values, unit_vectors=None, scale_min=None, scale_max=None):
        s_values = values
        scaling_phenos = [x for x in range(values.shape[1]) if x not in unit_vectors]
        if not scaling_phenos:
            return values
        scaler = Scaler(X_min=np.zeros(len(scaling_phenos)), X_max=np.ones(len(scaling_phenos)))
        if scale_min is not None:
            scaler.X_min = scale_min
        if scale_max is not None:
            scaler.X_max = scale_max
        s_values[:, scaling_phenos] = scaler.fit_transform(s_values[:, scaling_phenos])
        return s_values
