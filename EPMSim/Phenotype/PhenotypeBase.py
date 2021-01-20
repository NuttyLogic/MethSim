from typing import Tuple


class PhenotypeBase:
    """Abstract base class for simulation phenotypes

    Attributes:

    * self.representation (str)*: function representation of phenotype behavior
    * self.binary (int)*: sample has phenotype
    * self.mean (float)*: mean phenotype value
    * self.std (std)*: phenotype standard deviation
    * self.health_effect (bool)*: trait affects sample health
    * self.health_dist (bool)*: phenotype sampling distribution effected by sample health
        """

    def __init__(self):
        self.representation = 'NA'
        self.binary = None
        self.mean = None
        self.std = None
        self.health_effect = None
        self.health_dist = None
        self.age_repr = None

    def __str__(self):
        return self.representation

    def get_phenotype(self, age: float, health: float) -> Tuple[int, float, float]:
        """Return phenotype status, phenotype value, expected phenotype value"""
        pass

    def age_association(self, age: float, phenotype: float):
        """Return phenotype age association"""
        pass
