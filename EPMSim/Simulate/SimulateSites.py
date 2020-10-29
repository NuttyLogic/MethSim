import random
from typing import Dict, List, Tuple
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from EPMSim.MethylationSite.SiteGenerator import SiteGenerator
from EPMSim.MethylationSite.ValueGenerator import generate_sample_methylation
from EPMSim.Simulate.SimulateSamples import retrieve_phenotypes


class GenerateSites:

    def __init__(self,
                 site_conditions: Dict[str, Tuple[List[str], int]],
                 phenotype_weights: Dict[str, float] = None,
                 noise_sites: int = 5000,
                 phenotype_site_generators: Dict[str, SiteGenerator] = None,
                 noise_site_generators: Dict[str, SiteGenerator] = None):
        self.site_conditions = site_conditions
        self.phenotype_weights = phenotype_weights if phenotype_weights else self.gen_phenotype_weights
        self.noise_sites = noise_sites
        self.pheno_site_generators = phenotype_site_generators if phenotype_site_generators else {}
        self.noise_site_generators = noise_site_generators if noise_site_generators else {}
        self.get_site_generators()

    @property
    def gen_phenotype_weights(self):
        weights = {}
        for condition, count in self.site_conditions.values():
            for pheno in condition:
                if pheno not in weights:
                    weights[pheno] = 1.0
        return weights

    def get_site_generators(self):
        deviations = [0.01, 0.025, 0.05]
        if not self.pheno_site_generators:
            for deviation in deviations:
                self.pheno_site_generators[f'phenotype_{deviation}'] = SiteGenerator(site_deviation=deviation)
        if not self.noise_site_generators:
            for deviation in deviations:
                self.noise_site_generators[f'noise_{deviation}'] = SiteGenerator(max_delta=0.0,
                                                                                 site_deviation=deviation + 0.05)

    def generate_matrix(self, samples) -> Tuple[dict, np.ndarray, np.ndarray]:
        meth_matrix = []
        sample_phenotypes, phenotype_key = retrieve_phenotypes(samples)
        scaler = MinMaxScaler()
        for site_cond_label, site_cond_info in self.site_conditions.items():
            site_phenos, site_count = site_cond_info
            sample_values = sum(np.array([sample_phenotypes[:, phenotype_key[pheno]] * self.phenotype_weights[pheno] for
                                          pheno in site_phenos]))
            scaled_values = [x[0] for x in scaler.fit_transform(sample_values.reshape(-1, 1))]
            meth_matrix.extend(self.generate_sites(scaled_values, self.pheno_site_generators,
                                                   site_count, site_phenos))
        meth_matrix.extend(self.generate_sites([0.0 for _ in range(len(samples))],
                                               self.noise_site_generators,
                                               self.noise_sites, ['noise']))
        random.shuffle(meth_matrix)
        site_info = {count: site[0] for count, site in enumerate(meth_matrix)}
        meth_values = np.array([site[1] for site in meth_matrix])
        meth_errors = np.array([site[2] for site in meth_matrix])
        return site_info, meth_values, meth_errors

    @staticmethod
    def generate_sites(phenotype_values: np.array,
                       site_generators: Dict[str, SiteGenerator],
                       site_count: int, site_info: List[str]) -> List[Tuple[dict, np.ndarray, np.ndarray]]:
        meth_sites = []
        for _ in range(site_count):
            generator_label = random.choice(list(site_generators.keys()))
            site_gen = site_generators[generator_label]
            m_not_i, rate, std = site_gen.simulate_site()
            meth_values, meth_error = [], []
            for value in phenotype_values:
                meth, error = generate_sample_methylation(value, m_not_i, rate, std)
                meth_values.append(meth)
                meth_error.append(error)
            site_meta = dict(site_phenotypes=sorted(site_info), generator=generator_label,
                             m_not=m_not_i, rate=rate, std=std)
            meth_sites.append((site_meta,
                              np.array(meth_values), np.array(meth_error)))
        return meth_sites
