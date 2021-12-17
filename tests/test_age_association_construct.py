from methsim.phenotype import single_exp_assoc

sqrt_time_assoc = single_exp_assoc(age_weight=1, age_exponent=0.5)
no_time_assoc = single_exp_assoc(age_weight=1, age_exponent=0)
linear_time_assoc = single_exp_assoc(age_weight=1, age_exponent=1.0)


def test_default_association_init():
    assert sqrt_time_assoc(4, 4) == 8
    assert no_time_assoc(4, 4) == 4
    assert linear_time_assoc(4, 4) == 16


