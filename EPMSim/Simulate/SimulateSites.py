import random
from typing import Dict, List, Tuple
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from EPMSim.MethylationSite.SiteGenerator import SiteGenerator
from EPMSim.MethylationSite.ValueGenerator import generate_sample_methylation
from EPMSim.Simulate.SimulateSamples import retrieve_phenotypes


class GenerateSites:

    def __init__(self, phenotype_weights: Dict[str, float],
                 site_conditions: Dict[str, Tuple[List[str], int]],
                 noise_sites: int = 5000):
        self.phenotype_weights = phenotype_weights
        self.site_conditions = site_conditions
        self.noise_sites = noise_sites
        self.pheno_site_generators = {}
        self.noise_site_generators = {}
        self.get_site_generators()

    def get_site_generators(self):
        for deviation in [0.1, 0.15, 0.2]:
            self.pheno_site_generators[f'phenotype_{deviation}'] = SiteGenerator(site_deviation=deviation)
            self.noise_site_generators[f'noise_{deviation}'] = SiteGenerator(max_delta=0.0,
                                                                             site_deviation=deviation + 0.05)

    def generate_matrix(self, samples):
        meth_matrix = []
        sample_phenotypes, phenotype_key = retrieve_phenotypes(samples)
        scaler = MinMaxScaler()
        for site_cond_label, site_cond_info in self.site_conditions.items():
            site_phenos, site_count = site_cond_info
            sample_values = sum(np.array([sample_phenotypes[:, phenotype_key[pheno]] * self.phenotype_weights[pheno] for
                                          pheno in site_phenos]))
            scaled_values = scaler.fit_transform(sample_values.reshape(-1, 1))
            for _ in range(site_count):
                generator_label = random.choice(list(self.pheno_site_generators.keys()))
                site_gen = self.pheno_site_generators[generator_label]
                m_not_i, rate, std = site_gen.simulate_site()
                meth_values, meth_error = [], []
                for value in scaled_values:
                    meth, error = generate_sample_methylation(value[0], m_not_i, rate, std)
                    meth_values.append(meth)
                    meth_error.append(error)
                meth_matrix.append(((site_cond_info, generator_label, m_not_i, rate, std),
                                    np.array(meth_values), np.array(meth_error)))
        for _ in range(self.noise_sites):
            generator_label = random.choice(list(self.noise_site_generators.keys()))
            site_gen = self.noise_site_generators[generator_label]
            m_not_i, rate, std = site_gen.simulate_site()
            meth_values, meth_error = [], []
            for _ in range(len(samples)):
                meth, error = generate_sample_methylation(0.0, m_not_i, rate, std)
                meth_values.append(meth)
                meth_error.append(error)
            meth_matrix.append((('noise', generator_label, m_not_i, rate, std),
                                np.array(meth_values), np.array(meth_error)))
        random.shuffle(meth_matrix)
        site_info = {count: site[0] for count, site in enumerate(meth_matrix)}
        meth_values = np.array([site[1] for site in meth_matrix])
        meth_errors = np.array([site[2] for site in meth_matrix])
        return site_info, meth_values, meth_errors




