import numpy as np

from EPMSim.Phenotype.AgeAssociations import construct_age_association
from EPMSim.Phenotype.BinaryPhenotype import BinaryPhenotype
from EPMSim.Phenotype.ContinuousPhenotype import ContinuousPhenotype

sqrt_age_assoc, sqrt_age_assoc_repr = construct_age_association(age_exponent=0.5, return_repr=True)

phenotypes = []

for _ in range(3):
    mean = np.random.uniform(0.0, 2.0)
    bp = BinaryPhenotype(mean=mean, std=0.0, binary_prob=0.5,
                         age_association=sqrt_age_assoc, age_repr=sqrt_age_assoc_repr,
                         health_dist=True)
    cp = ContinuousPhenotype(mean=mean, std=0.0, age_association=sqrt_age_assoc,
                             age_repr=sqrt_age_assoc_repr, health_dist=True)
    phenotypes.extend([bp, cp])


def test_phenotypes():
    for phenotype in phenotypes:
        for _ in range(100):
            age = np.random.uniform(1, 100)
            health = np.random.uniform(-1.0, 1.0)
            has_pheno, pheno_value = phenotype.get_phenotype(age, health)
            if has_pheno:
                assert pheno_value == (phenotype.mean + health) * np.sqrt(age)
            else:
                assert isinstance(phenotype, BinaryPhenotype)
                assert pheno_value == 0.0

