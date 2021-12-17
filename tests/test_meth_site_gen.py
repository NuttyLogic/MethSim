import numpy as np
from methsim.site import generate_sample_methylation, simulate_site, simulate_site_matrix

phenotypes = np.linspace(0, 1, 100).reshape(-1, 1)

site = simulate_site(0.00, .2)
meth_values, error = generate_sample_methylation(phenotypes, *site)

print(meth_values.shape)
