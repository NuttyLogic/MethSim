import numpy as np

from methsim.phenotype import single_exp_assoc, continuous_normal
from methsim.site import GenerateSampleMethylation
from methsim.sample import simulate_samples

sqrt_time_association = single_exp_assoc(age_weight=1.0, age_exponent=0.5)
continuous_pheno = continuous_normal(mean=1.0, std=0.0, time_association=sqrt_time_association)

phenotypes = {'continuous': continuous_pheno}

sim_samples = simulate_samples(min_age=0, max_age=100,
                               phenotypes=phenotypes,
                               sample_count=500)

meth_phenotype_counts = {'continuous': 100,
                         'noise': 100}

site_gen = GenerateSampleMethylation()
site_gen.generate_methylation_sites(deviation_low=0.00, deviation_high=0.00,
                                    delta_low=0.1, delta_high=0.4, number_of_sites=500)

meth_matrix, meth_error = site_gen.generate_sample_methylation(sim_samples, meth_phenotype_counts)



def test_sample_setting():
    for sample in sim_samples.values():
        assert sample.continuous == sample.age ** 0.5



