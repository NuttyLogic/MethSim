

def single_exp_assoc(age_weight=1.0, age_exponent=0.5):
    """"""
    def time_association(age=0.0, phenotype=0.0, return_var=False):
        if return_var:
            return {'age_weight': age_weight, 'age_exponent':age_exponent}
        return age_weight * age ** age_exponent * phenotype
    return time_association


def double_exp_assoc(age_weight=1.0, exp1=0.5, exp2=0.5):
    """"""
    def time_association(age=0.0, phenotype=0.0, return_var=False):
        if return_var:
            return {'age_weight': 1.0, 'exp1': exp1, 'exp2': exp2}
        return age_weight * age ** (exp1 - exp2) * phenotype
    return time_association
