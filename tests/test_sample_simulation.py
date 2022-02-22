import numpy as np
from methsim.phenotype import single_exp_assoc, continuous_normal
from methsim.site import GenerateSampleMethylation
from methsim.sample import gen_pheno_matrix, simulate_samples
from methsim.utilities import pearson_correlation

sqrt_time_association = single_exp_assoc(age_weight=1.0, age_exponent=0.5)
continuous_pheno = continuous_normal(mean=1.0, std=0.0, time_association=sqrt_time_association)

phenotypes = {'continuous': continuous_pheno}

sim_samples = simulate_samples(min_age=0, max_age=100,
                               phenotypes=phenotypes,
                               sample_count=500)
sim_samples_2 = simulate_samples(min_age=0, max_age=100,
                                 phenotypes=phenotypes,
                                 sample_count=250)

meth_phenotype_counts = [[['continuous'], 100, np.ones(1)],
                         ['noise', 100, np.ones(1)]]

sim_continuous = gen_pheno_matrix(sim_samples, ['continuous'])

site_gen = GenerateSampleMethylation()
site_gen.generate_methylation_sites(deviation_low=0.001, deviation_high=0.00,
                                    delta_low=0.1, delta_high=0.2, number_of_sites=50)

meth_matrix, meth_error = site_gen.generate_sample_methylation(sim_samples, meth_phenotype_counts)

# generate second set of methylation matrices
meth_matrix_2, meth_error_2 = site_gen.generate_sample_methylation(sim_samples_2, meth_phenotype_counts)


site_corr = pearson_correlation(sim_continuous, meth_matrix)


def test_sample_setting():
    for sample in sim_samples.values():
        assert sample.continuous == sample.age ** 0.5


def test_continuous_site_corr():
    for corr in site_corr[0, 0:100]:
        assert abs(corr) > .98


def test_noise_site_corr():
    for corr in site_corr[0, 101:]:
        assert abs(corr) < .5

