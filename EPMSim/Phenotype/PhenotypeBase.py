from typing import Tuple


class PhenotypeBase:

    def __init__(self):
        self.representation = 'NA'
        self.binary = None
        self.mean = None
        self.std = None
        self.health_effect = None
        self.age_repr = None

    def __repr__(self):
        return self.representation

    def get_phenotype(self, age: float, health: float) -> Tuple[int, float]:
        pass

    def age_association(self, age: float, phenotype: float):
        pass
