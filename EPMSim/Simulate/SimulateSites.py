import random
from typing import Dict, List, Tuple, Union
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
        self.sites = None
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

    def generate_sites(self):
        self.sites = self.get_condition_sites(self.site_conditions, self.pheno_site_generators)
        self.sites.extend(self.get_condition_sites({'noise': (['noise'], self.noise_sites)},
                                                   self.noise_site_generators))

    def generate_sample_methylation(self, samples) -> Union[None, Tuple[np.ndarray, np.ndarray]]:
        meth_matrix, meth_error = [], []
        sample_phenotypes, phenotype_key = retrieve_phenotypes(samples)
        scaler = MinMaxScaler()
        if not self.sites:
            print('must generate sites')
            return None
        for site in self.sites:
            sample_values = sum(np.array([sample_phenotypes[:, phenotype_key[pheno]] * self.phenotype_weights[pheno] for
                                          pheno in site['site_phenotypes']]))
            scaled_values = [x[0] for x in scaler.fit_transform(sample_values.reshape(-1, 1))]
            values, error = self.generate_site_values(scaled_values, site)
            meth_matrix.append(values)
            meth_error.append(error)
        return np.array(meth_matrix), np.array(meth_error)

    @staticmethod
    def get_condition_sites(site_conditions: Dict[str, Tuple[List[str], int]],
                            site_generators: Dict[str, SiteGenerator]) -> List[Dict[str, Union[List[str], str, float]]]:
        meth_sites = []
        for cond_label, site_info in site_conditions.items():
            for _ in range(site_info[1]):
                generator_label = random.choice(list(site_generators.keys()))
                site_gen = site_generators[generator_label]
                m_not_i, rate, std = site_gen.simulate_site()
                meth_sites.append(dict(site_phenotypes=sorted(site_info[0]), generator=generator_label,
                                       m_not=m_not_i, rate=rate, std=std))
        return meth_sites

    @staticmethod
    def generate_site_values(phenotype_values: np.array,
                             site: Dict) -> Tuple[np.ndarray, np.ndarray]:
        meth_values, meth_error = [], []
        for value in phenotype_values:
            meth, error = generate_sample_methylation(value, site['m_not'], site['rate'], site['std'])
            meth_values.append(meth)
            meth_error.append(error)
        return np.array(meth_values), np.array(meth_error)
