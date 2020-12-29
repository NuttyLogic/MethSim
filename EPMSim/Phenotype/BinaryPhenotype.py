from typing import Tuple

import numpy as np

from EPMSim.Phenotype.PhenotypeBase import PhenotypeBase
from EPMSim.Phenotype.AgeAssociations import construct_age_association

default_age_association, default_age_repr = construct_age_association(return_repr=True)


class BinaryPhenotype(PhenotypeBase):
    """
    Params:
        * *representation (str)*:
        * *param mean (float)*:
        * *param std (float)*:
        * *param binary_prob (float)*:
        * *age_limit (float)*:
        * *age_association (EPMSim.Phenotypes.PhenotypeAgeAssociation)*:
        * *age_repr (str)*:
        * *health_effect (bool)*:
        * *health_dist (bool)*:
    """

    def __init__(self, representation='NA', mean: float = 1.0,
                 std: float = 0.05,
                 binary_prob: float = 0.5,
                 age_limit: float = 0.0,
                 age_association=default_age_association,
                 age_repr: str = default_age_repr,
                 health_effect=True,
                 health_dist=True):
        PhenotypeBase.__init__(self)
        self.representation = representation
        self.mean = mean
        self.std = std
        self.binary = 1
        self.binary_prob = binary_prob
        self.age_limit = age_limit
        self.age_association = age_association
        self.age_repr = age_repr
        self.health_dist = health_dist
        self.health_effect = health_effect

    def __str__(self):
        """Return user defined phenotype representation"""
        return self.representation

    def get_phenotype(self, age: float, health: float) -> Tuple[int, float]:
        """Get sample phenotype information given age and health
        Params:
            * *age (float)*:
            * *health (float)*:
        Returns:
            * *Phenotype Info (Tuple[int, float])*: 1 if sample has trait 1 else 0,
                                                    phenotype value if sample has phenotype
            """
        has_trait = 0
        if age > self.age_limit:
            has_trait = 1 if np.random.uniform(0.0, 1.0) > self.binary_prob else 0
        if self.health_dist:
            phenotype = np.random.normal(loc=self.mean + health, scale=self.std) if has_trait else 0.0
        else:
            phenotype = np.random.normal(loc=self.mean, scale=self.std) if has_trait else 0.0
        return has_trait, self.age_association(age, phenotype)
