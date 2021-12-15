from typing import Union


class Sample:

    def __init__(self, age: float = 0.0, health: float = 0.0):
        self.age = float(age)
        self.health = float(health)

    def add_phenotype(self, phenotype: str, value: Union[int, float]):
        setattr(self, phenotype, value)

