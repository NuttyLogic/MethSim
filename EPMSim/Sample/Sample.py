from typing import Dict


class Sample:

    def __init__(self, age: float = 0.0, health: float = 0.0, epigenetic_state=lambda x: 10.0 * x**0.5):
        self.age = float(age)
        self.health = float(health)
        self.epigenetic_state = epigenetic_state(age)
        self.expected_epigenetic_state = epigenetic_state(age)
        self.phenotypes = {}

    def add_phenotype(self, phenotype_name: str, phenotype_info: Dict):
        self.phenotypes[phenotype_name] = phenotype_info
        if phenotype_info['health_effect']:
            self.epigenetic_state += phenotype_info['trait_value']
