from typing import Dict
import numpy as np
from methsim.site import simulate_site_matrix, generate_sample_methylation
from methsim.sample import gen_pheno_matrix
from methsim.utilities import Scaler


class GenerateSampleMethylation:

    def __init__(self):
        self.methylation_sites = None
        self.sim_site_order = {}
        self.scaler = Scaler(X_min=np.array([0.0]), X_max=np.array([1.0]))

    def generate_methylation_sites(self, deviation_low=0.0, deviation_high=0.1,
                                   delta_low=0.0, delta_high=0.5, number_of_sites=5000):
        self.methylation_sites = simulate_site_matrix(deviation_low, deviation_high,
                                                      delta_low, delta_high,
                                                      number_of_sites)

    def generate_sample_methylation(self, samples, sim_phenotypes=Dict[str, int]):
        site_values, site_error = [], []
        site = 0
        for phenotype, sim_sites in sim_phenotypes.items():
            pheno_values = self.get_phenotype_values(samples, phenotype)
            for count in range(sim_sites):
                if site in self.sim_site_order:
                    meth_site = self.sim_site_order[site][0]
                else:
                    meth_site = np.random.randint(0, self.methylation_sites.shape[0])
                    self.sim_site_order[site] = (site, phenotype)
                site += 1
                meth, error = generate_sample_methylation(pheno_values, * self.methylation_sites[meth_site])
                site_values.append(meth)
                site_error.append(error)
        return np.array(site_values), np.array(site_values)

    def get_phenotype_values(self, samples, phenotype):
        values = gen_pheno_matrix(samples, [phenotype])
        if np.isnan(np.sum(values)):
            print(f'{phenotype} not set in all samples, returning unit vector')
            return np.ones((len(samples), 1))
        return self.scaler.transform(values)
