import numpy as np
from methsim.phenotype.time import single_exp_assoc

rng = np.random.default_rng()


def continuous_normal(mean, std, health_association=False, time_association=single_exp_assoc()):
    def con_pheno(age=0.0, health=0.0, return_var=False):
        if return_var:
            pheno_info = {'mean': mean, 'std': std, 'health_assoc': health_association}
            pheno_info.update(time_association(return_var=True))
            return pheno_info
        pheno_mean = mean + health if health_association else mean
        pheno = rng.normal(loc=pheno_mean, scale=std)
        return 1, pheno, time_association(age, pheno)
    return con_pheno


def binary_normal(mean, std, health_association=False, binary_prob=0.5,
                  age_limit=0, time_association=single_exp_assoc()):
    def binary_pheno(age=0.0, health=0.0, return_var=False):
        if return_var:
            pheno_info = {'mean': mean, 'std': std, 'health_assoc': health_association,
                          'binary_prob': binary_prob, 'age_limit': age_limit}
            pheno_info.update(time_association(return_var=True))
            return pheno_info
        pheno_mean = mean + health if health_association else mean
        has_trait = 0
        if age > age_limit:
            has_trait = 1 if np.random.uniform(0.0, 1.0) <= binary_prob else 0
        pheno = 1.0 if not has_trait else rng.normal(loc=pheno_mean, scale=std)
        return has_trait, pheno, time_association(age, pheno)
    return binary_pheno


def continuous_uniform(low, high, time_association=single_exp_assoc()):
    def con_uniform_pheno(age=0.0, health=0.0, return_var=False):
        if return_var:
            pheno_info = {'low': low, 'high': high}
            pheno_info.update(time_association(return_var=True))
            return pheno_info
        pheno = rng.uniform(low, high)
        return 1.0, pheno, time_association(age, pheno)
    return con_uniform_pheno


def continuous_exponential(scale=1.0, time_association=single_exp_assoc()):
    def con_exponential_pheno(age=0.0, health=0.0, return_var=False):
        if return_var:
            pheno_info = {'scale': scale}
            pheno_info.update(time_association(return_var=True))
            return pheno_info
        pheno = rng.exponential(scale)
        return 1.0, pheno, time_association(age, pheno)
    return con_exponential_pheno
