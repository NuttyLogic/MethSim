

class SamplePhenotype:
    
    def __init__(self, has_trait: int = 1, trait_value: float = 1.0,
                 expected_trait_value: float = 1.0, binary: int = 0,
                 health_effect: bool = False):
        self.has_trait = has_trait
        self.trait_value = trait_value
        self.expected_trait_value = expected_trait_value
        self.binary = binary
        self.health_effect = health_effect

    def __add__(self, other):
        if isinstance(other, SamplePhenotype):
            return self.trait_value + other.trait_value
        return self.trait_value + other

    def __mul__(self, other):
        if isinstance(other, SamplePhenotype):
            return self.trait_value * other.trait_value
        return self.trait_value * other

    def __eq__(self, other):
        if isinstance(other, SamplePhenotype):
            return self.trait_value == other.trait_value
        return self.trait_value == other

    def __ge__(self, other):
        if isinstance(other, SamplePhenotype):
            return self.trait_value >= other.trait_value
        return self.trait_value >= other

    def __gt__(self, other):
        if isinstance(other, SamplePhenotype):
            return self.trait_value > other.trait_value
        return self.trait_value > other

    def __le__(self, other):
        if isinstance(other, SamplePhenotype):
            return self.trait_value <= other.trait_value
        return self.trait_value <= other

    def __lt__(self, other):
        if isinstance(other, SamplePhenotype):
            return self.trait_value < other.trait_value
        return self.trait_value < other
