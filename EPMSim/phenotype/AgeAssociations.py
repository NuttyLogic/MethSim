

class AgeAssociation:

    def __init__(self, age_weight: float = 1.0, age_exponent: float = 0.0):
        """age_phenotype = (age_weight)age**(age_exp) * phenotype"""
        self.age_weight = age_weight
        self.age_exponent = age_exponent

    def age_association(self, age: float, phenotype: float):
        return self.age_weight * age ** self.age_exponent * phenotype


def construct_age_association(age_weight: float = 1.0, age_exponent: float = 0.0):
    age_assoc = AgeAssociation(age_weight=age_weight, age_exponent=age_exponent)
    return age_assoc.age_association
