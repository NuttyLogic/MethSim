import numpy as np
from EPMSim.phenotype.AgeAssociations import construct_age_association, PhenotypeAgeAssociation

default_phenotype = PhenotypeAgeAssociation(age_weight=1.0, age_exponent=0.0)

default_age_association, default_repr = construct_age_association(return_repr=True)
default_age_association_no_repr = construct_age_association()
sqrt_age_association, sqrt_age_repr = construct_age_association(age_weight=1.0, age_exponent=0.5, return_repr=True)


def test_default_association_init():
    assert default_age_association(1.0, 1.0) == default_age_association_no_repr(1.0, 1.0) == 1.0
    assert default_repr == '$p_a = 1.00\\times a^{0.00}\\times q$' == default_phenotype.phenotype_representation


def test_phenotype_values():
    for age in [x + 1.0 for x in range(100)]:
        assert np.sqrt(age) == sqrt_age_association(age, 1.0)
        # default has no age association so always 1
        assert 1.0 == default_age_association(age, 1.0)
        # if age is the phenotype then age == pheno
        assert age == default_age_association_no_repr(age, age)
