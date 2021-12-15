from typing import Dict, List, Tuple
import numpy as np
from EPMSim.phenotype.PhenotypeBase import PhenotypeBase
from EPMSim.Sample import Sample, SamplePhenotype


def simulate_samples(min_age: float = 0.0, max_age: float = 100.0,
                     health_variation: float = 0.05, sample_count: int = 1000,
                     phenotypes: List[PhenotypeBase] = None):
    sim_samples = {}
    for sample in range(sample_count):
        age = np.random.uniform(min_age, max_age)
        health = np.random.normal(0, health_variation)
        sim_samples[sample] = Sample(age=age, health=health, epigenetic_state=epigenetic_state)
    for sample in range(sample_count):
        current_sample = sim_samples[sample]
        if phenotypes:
            set_sample_phenotypes(current_sample, phenotypes=phenotypes)
    return sim_samples


def set_sample_phenotypes(sample: Sample, phenotypes: List[PhenotypeBase] = None):
    for count, phenotype in enumerate(phenotypes):
        pheno_repr = str(phenotype) if str(phenotype) != 'NA' else str(count)
        has_trait, trait_value, expected_trait_value = phenotype.get_phenotype(sample.age,
                                                                               sample.health)
        pheno_info = SamplePhenotype(has_trait=has_trait, trait_value=trait_value,
                                     expected_trait_value=expected_trait_value,
                                     binary=phenotype.binary, health_effect=phenotype.health_effect,
                                     age_repr=phenotype.age_repr, pheno_repr=pheno_repr)
        sample.add_phenotype(pheno_info)


def retrieve_phenotypes(samples: Dict[str, Sample]) -> Tuple[np.array, Dict[str, int]]:
    phenotypes = list(next(iter(samples.values())).phenotypes.keys())
    phenotype_values = np.zeros((len(samples), len(phenotypes)))
    for sample_index, sample_info in enumerate(samples.values()):
        phenotype_values[sample_index, :] = np.array([sample_info.phenotypes[ph].trait_value for ph in phenotypes])
    return phenotype_values, {phenotype: index for index, phenotype in enumerate(phenotypes)}
