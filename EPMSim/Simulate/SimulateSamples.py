from typing import Dict, List, Tuple
import numpy as np
from EPMSim.Phenotype.PhenotypeBase import PhenotypeBase
from EPMSim.Sample import Sample


def simulate_samples(min_age: float = 0.0, max_age: float = 100.0,
                     health_variation: float = 0.05, sample_count: int = 1000,
                     epigenetic_state=lambda x: 10.0 * x**0.5, phenotypes: List[PhenotypeBase] = None):
    sim_samples = {}
    for sample in range(sample_count):
        age = np.random.uniform(min_age, max_age)
        health = np.random.normal(0, health_variation)
        sim_samples[sample] = Sample(age=age, health=health, epigenetic_state=epigenetic_state)
    for sample in range(sample_count):
        current_sample = sim_samples[sample]
        if phenotypes:
            for count, phenotype in enumerate(phenotypes):
                pheno_repr = str(phenotype) if str(phenotype) != 'NA' else str(count)
                has_trait, trait_value = phenotype.get_phenotype(current_sample.age,
                                                                 current_sample.health)
                pheno_info = dict(has_trait=has_trait, trait_value=trait_value,
                                  binary=phenotype.binary, health_effect=phenotype.health_effect,
                                  age_repr=phenotype.age_repr, pheno_repr=str(phenotype))
                current_sample.add_phenotype(pheno_repr, pheno_info)
    return sim_samples


def retrieve_phenotypes(samples: Dict[str, Sample]) -> Tuple[np.array, Dict[str, int]]:
    phenotypes = list(next(iter(samples.values())).phenotypes.keys())
    phenotype_values = np.zeros((len(samples), len(phenotypes)))
    for sample_index, sample_info in enumerate(samples.values()):
        for phenotype_index, phenotype in enumerate(phenotypes):
            phenotype_values[sample_index, phenotype_index] = sample_info.phenotypes[phenotype]['trait_value']
    return phenotype_values, {phenotype: index for index, phenotype in enumerate(phenotypes)}
