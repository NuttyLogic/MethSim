from typing import Tuple

import numpy as np

from EPMSim.phenotype.AgeAssociations import construct_age_association
from EPMSim.phenotype.PhenotypeBase import PhenotypeBase

default_age_association, default_age_repr = construct_age_association(return_repr=True)


class ContinuousPhenotype(PhenotypeBase):
    """
       :param mean:
       :param std:
       :param age_association:
    """

    def __init__(self, representation='NA', mean: float = 1.0,
                 std: float = 0.05,
                 age_association=default_age_association,
                 age_repr: str = default_age_repr,
                 health_effect=True,
                 health_dist=True):
        PhenotypeBase.__init__(self)
        self.representation = representation
        self.mean = mean
        self.std = std
        self.binary = 0
        self.age_association = age_association
        self.age_repr = age_repr
        self.health_effect = health_effect
        self.health_dist = health_dist

    def __str__(self):
        return self.representation

    def get_phenotype(self, age: float, health: float) -> Tuple[int, float, float]:
        if self.health_dist:
            phenotype = np.random.normal(loc=self.mean + health, scale=self.std)
        else:
            phenotype = np.random.normal(loc=self.mean, scale=self.std)
        return 1, self.age_association(age, phenotype), self.age_association(age, self.mean)
