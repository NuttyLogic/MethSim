from methsim.celltype import get_ct_meth_matrix
from methsim.celltype import set_blood_ct_comp
from methsim.phenotype import single_exp_assoc, continuous_normal
from methsim.site import GenerateSampleMethylation
from methsim.sample import gen_pheno_matrix, simulate_samples

import numpy as np

sqrt_time_association = single_exp_assoc(age_weight=1.0, age_exponent=0.5)
continuous_pheno = continuous_normal(mean=1.0, std=0.0, time_association=sqrt_time_association)

phenotypes = {'continuous': continuous_pheno}

sim_samples = simulate_samples(min_age=0, max_age=100,
                               phenotypes=phenotypes,
                               sample_count=50)

meth_phenotype_counts = {'continuous': 100,
                         'noise': 100}

sim_continuous = gen_pheno_matrix(sim_samples, ['continuous'])

site_gen = GenerateSampleMethylation()
site_gen.generate_methylation_sites(deviation_low=0.001, deviation_high=0.00,
                                    delta_low=0.1, delta_high=0.2, number_of_sites=500)

meth_matrix, meth_error = site_gen.generate_sample_methylation(sim_samples, meth_phenotype_counts)

set_blood_ct_comp(sim_samples)

test_matrix, test_ct_ref = get_ct_meth_matrix(meth_matrix, sim_samples, ct_pheno_prop=.1)


def test_ct_sim():
    for sample in range(len(sim_samples)):
        cts = np.array([sim_samples[sample].pheno_detail[cat]['value'] for cat in
                        list(dir(sim_samples[sample])) if '_ct' in cat])
        for site_val, ref, val in zip(meth_matrix[:, sample], test_ct_ref, test_matrix[:, sample]):
            site_ct_contrib = np.array([ref[x] if not np.isnan(ref[x]) else site_val for x in range(len(cts))])
            assert np.isclose(sum(site_ct_contrib * cts), val)

