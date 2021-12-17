from typing import Dict
import numpy as np


class Sample:

    def __init__(self, age: float = 0.0, health: float = 0.0,
                 sample_id: int = 0):
        self.age = float(age)
        self.health = float(health)
        self.sample_id = sample_id
        self.pheno_detail = {'age': {'has_trait': 1, 'exposure': 1.0, 'value': age},
                             'health': {'has_trait': 1, 'exposure': 1.0, 'value': health}}

    def add_phenotype(self, phenotype_label: str, phenotype):
        has_trait, exposure, value = phenotype(self.age, self.health)
        if hasattr(self, phenotype_label):
            raise Exception(f'Sample {phenotype_label} already set')
        setattr(self, phenotype_label, value)
        self.pheno_detail[phenotype_label] = {'has_trait': has_trait, 'exposure': exposure, 'value': value}

    def add_phenotypes(self, phenotypes):
        for phenotype_label, phenotype in phenotypes.items():
            self.add_phenotype(phenotype_label, phenotype)


def simulate_samples(min_age: float = 0.0, max_age: float = 100.0,
                     health_variation: float = 0.05, sample_count: int = 1000,
                     phenotypes=None):
    sim_samples = {}
    for sample_id in range(sample_count):
        age = np.random.uniform(min_age, max_age)
        health = np.random.normal(0, health_variation)
        sample = Sample(age=age, health=health, sample_id=sample_id)
        if phenotypes:
            sample.add_phenotypes(phenotypes)
        sim_samples[sample_id] = sample
    return sim_samples


def gen_pheno_matrix(samples, phenotypes: list = None):
    pheno_values = np.zeros((len(samples), len(phenotypes)))
    for count, sample in enumerate(samples.values()):
        pheno_values[count, :] = np.array([sample.pheno_detail.get(pheno, {}).get('value', np.nan)
                                           for pheno in phenotypes])
    return pheno_values
