from methsim.phenotype import single_exp_assoc
from methsim.phenotype import binary_normal, continuous_normal

linear_time_assoc = single_exp_assoc(age_weight=1, age_exponent=1.0)

binary_pheno1 = binary_normal(mean=.5, std=0.0, binary_prob=1.0, time_association=linear_time_assoc)

# 1, 0.5, 5.0, test pheno setting
binary_value_1 = binary_pheno1(10, 1)

binary_pheno_ha = binary_normal(mean=.5, std=0.0, binary_prob=1.0, time_association=linear_time_assoc,
                                health_association=True, age_limit=5)
# 1, 2.0, 20, test_health_assoc
binary_value_2 = binary_pheno_ha(10, 1.5)

# 0, 1.0, 1.0, test_age_limit
binary_value_3 = binary_pheno_ha(1, 1)

binary_non_carrier = binary_normal(mean=.5, std=0.0, binary_prob=0.0, time_association=linear_time_assoc,
                                   health_association=True)
# 0, 1, 10, test non carriers
binary_value_4 = binary_non_carrier(10, 11)

#


con1 = continuous_normal(mean=1, std=0.0, time_association=linear_time_assoc)
# 1, 1, 1
con_value_1 = con1(1, 5)
# 1, 1, 10
con_value_2 = con1(10, 5)
# 1, 1, 50
con_value_3 = con1(50, 5)

con2 = continuous_normal(mean=.5, std=0.0, time_association=linear_time_assoc, health_association=True)
# 1, 1, 1
con_value_4 = con1(1, .5)
# 1, 1, 10
con_value_5 = con1(10, .5)
# 1, 1, 50
con_value_6 = con1(50, .5)


def test_binary_phenotypes():
    assert binary_value_1 == (1, 0.5, 5.0)
    assert binary_value_2 == (1, 2.0, 20.0)
    assert binary_value_3 == (0, 1.0, 1.0)
    assert binary_value_4 == (0, 1.0, 10.0)


def test_continuous_phenotypes():
    assert con_value_1 == con_value_4 == (1, 1.0, 1.0)
    assert con_value_2 == con_value_5 == (1, 1.0, 10.0)
    assert con_value_3 == con_value_6 == (1, 1.0, 50.0)
