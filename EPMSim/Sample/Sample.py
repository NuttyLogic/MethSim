from EPMSim.Sample import SamplePhenotype


class Sample:

    def __init__(self, age: float = 0.0, health: float = 0.0, epigenetic_state=lambda x: 10.0 * x**0.5):
        self.age = float(age)
        self.health = float(health)
        self.epigenetic_state = epigenetic_state(age)
        self.expected_epigenetic_state = epigenetic_state(age)
        self.phenotypes = {}

    def add_phenotype(self, phenotype: SamplePhenotype):
        self.phenotypes[str(phenotype)] = phenotype
        if phenotype.health_effect:
            self.epigenetic_state += (self.expected_epigenetic_state - phenotype.trait_actual_expected_diff * self.expected_epigenetic_state)
