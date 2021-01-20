import numpy as np

from EPMSim.Phenotype.AgeAssociations import construct_age_association
from EPMSim.Phenotype.BinaryPhenotype import BinaryPhenotype
from EPMSim.Phenotype.ContinuousPhenotype import ContinuousPhenotype
from EPMSim.Simulate.SimulateSamples import simulate_samples

sqrt_age_assoc, sqrt_age_assoc_repr = construct_age_association(age_exponent=0.5, return_repr=True)

phenotypes = {}

for count in range(3):
    mean = np.random.uniform(0.0, 2.0)
    bp = BinaryPhenotype(mean=mean, std=0.0, binary_prob=0.5,
                         age_association=sqrt_age_assoc, age_repr=sqrt_age_assoc_repr,
                         health_dist=True, representation=f'{count}_bp')
    phenotypes[bp.representation] = bp
    cp = ContinuousPhenotype(mean=mean, std=0.0, age_association=sqrt_age_assoc,
                             age_repr=sqrt_age_assoc_repr, health_dist=True,
                             representation=f'{count}_cp')
    phenotypes[cp.representation] = cp


samples = simulate_samples(phenotypes=list(phenotypes.values()), health_variation=.1,
                           sample_count=100, epigenetic_state=lambda x: 10.0 * x**0.5)


def test_sample_expected_epigenetic_state():
    for sample_info in samples.values():
        assert sample_info.expected_epigenetic_state == 10.0 * sample_info.age**0.5


def test_sample_phenotypes():
    for sample_info in samples.values():
        for phenotype_label, phenotype_info in sample_info.phenotypes.items():
            pheno_mean = phenotypes[phenotype_label].mean
            trait_value = (pheno_mean + sample_info.health) * sample_info.age ** 0.5
            if phenotype_info.binary:
                if phenotype_info.has_trait:
                    assert phenotype_info.trait_value == trait_value
                else:
                    # if binary trait site only dependent on age
                    assert phenotype_info.trait_value == sample_info.age ** 0.5
            else:
                assert phenotype_info.trait_value == trait_value


def test_epigenetic_age():
    for sample_info in samples.values():
        all_trait_values = 0.0
        for phenotype_label, phenotype_info in sample_info.phenotypes.items():
            pheno_mean = phenotypes[phenotype_label].mean
            if phenotype_info.has_trait:
                trait_value = (pheno_mean + sample_info.health) * sample_info.age ** 0.5
                expected_value = pheno_mean * sample_info.age ** 0.5
            else:
                trait_value = sample_info.age ** 0.5
                expected_value = trait_value
            all_trait_values += trait_value - expected_value
        assert np.isclose(sample_info.epigenetic_state, sample_info.expected_epigenetic_state + all_trait_values)
