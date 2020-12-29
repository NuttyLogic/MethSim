

class PhenotypeAgeAssociation:

    def __init__(self, age_weight: float = 1.0, age_exponent: float = 0.0):
        """age_phenotype = (age_weight)age**(age_exp) * phenotype"""
        self.age_weight = age_weight
        self.age_exponent = age_exponent
        self.phenotype_representation = f'$p_a = {age_weight:.2f}\\times a^{{{age_exponent:.2f}}}\\times q$'

    def age_association(self, age: float, phenotype: float):
        return self.age_weight * age ** self.age_exponent * phenotype


def construct_age_association(age_weight: float = 1.0, age_exponent: float = 0.0, return_repr=False):
    age_assoc = PhenotypeAgeAssociation(age_weight=age_weight, age_exponent=age_exponent)
    if return_repr:
        return age_assoc.age_association, age_assoc.phenotype_representation
    return age_assoc.age_association
