import numpy as np
from methsim.sample import gen_pheno_matrix


def get_blood_ct_estimates_by_age(ages: np.ndarray = None, return_refs=False) -> np.ndarray:
    """Returns expected cell type composition of blood cell types by input age"""
    ref_cts = ('mono_ct', 'gran_ct', 'cd4t_ct', 'cd8t_ct', 'nk_ct', 'bcell_ct')
    age_matrix = np.ones((len(ages), 2))
    age_matrix[:, 0] = np.array(ages)
    ref_matrix = np.array([[-1.46906370e-04,  3.77031614e-04, -5.12501622e-04,
                            -3.69454688e-04,  9.20123006e-04, -6.56921275e-04],
                           [9.17914076e-02,  5.74901202e-01,  1.43470171e-01,
                            1.19758713e-01,  3.57934305e-02,  8.87064038e-02]])
    ct_estimates = np.dot(age_matrix, ref_matrix)
    if return_refs:
        return ct_estimates, ref_cts
    return ct_estimates


def set_blood_ct_comp(samples, alpha_scale=500, ct_estimator=get_blood_ct_estimates_by_age):
    sample_ages = np.array([sample.age for sample in samples.values()])
    ct_comps, cts = ct_estimator(sample_ages, return_refs=True)
    for sample, row in zip(samples.values(), ct_comps):
        for ct, ct_val in zip(cts, np.random.dirichlet(alpha_scale*row)):
            sample.set_phenotype(ct, 1.0, ct_val, ct_val)


def get_ct_site_contribution(site_values, ct_count=0, ct_pheno_prop=0.5):
    ct_contribution = np.random.uniform(size=ct_count) < ct_pheno_prop
    site_vals = site_values.reshape(-1,1) * np.ones((len(site_values), ct_count))
    ct_ref = np.zeros(ct_count)
    ct_ref[:] = np.nan
    for count, contrib in enumerate(ct_contribution):
        if not contrib:
            ct_val = np.random.beta(.1,.1)
            site_vals[:, count] = ct_val
            ct_ref[count] = ct_val
    return site_vals, ct_ref


def get_ct_meth_matrix(meth_matrix, samples, ct_pheno_prop=.8):
    cell_types = [x for x in dir(samples[list(samples.keys())[0]]) if '_ct' in x]
    cell_type_dist = gen_pheno_matrix(samples, cell_types)
    ct_adj_values, ct_contrib = np.zeros(meth_matrix.shape), np.zeros((meth_matrix.shape[0], len(cell_types)))
    for count, site in enumerate(meth_matrix):
        site_contribution, site_ref = get_ct_site_contribution(site, len(cell_types), ct_pheno_prop=ct_pheno_prop)
        values = np.sum(cell_type_dist * site_contribution, axis=1)
        ct_adj_values[count] = values
        ct_contrib[count] = site_ref
    return ct_adj_values, ct_contrib

