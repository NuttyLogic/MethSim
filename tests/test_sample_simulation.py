import numpy as np

from EPMSim.Phenotype.AgeAssociations import construct_age_association
from EPMSim.Phenotype.BinaryPhenotype import BinaryPhenotype
from EPMSim.Phenotype.ContinuousPhenotype import ContinuousPhenotype
from EPMSim.Sample.Sample import Sample
from EPMSim.Simulate.SimulateSamples import retrieve_phenotypes, simulate_samples

sqrt_age_assoc, sqrt_age_assoc_repr = construct_age_association(age_exponent=0.5, return_repr=True)

phenotypes = {}

for count in range(3):
    mean = np.random.uniform(0.0, 2.0)
    bp = BinaryPhenotype(mean=mean, std=0.0, binary_prob=0.5,
                         age_association=sqrt_age_assoc, age_repr=sqrt_age_assoc_repr,
                         health_dist=True)
    cp = ContinuousPhenotype(mean=mean, std=0.0, age_association=sqrt_age_assoc,
                             age_repr=sqrt_age_assoc_repr, health_dist=True)
    phenotypes[f'{count}_bp'] = bp
    phenotypes[f'{count}_cp'] = cp

sample_one = Sample(age=10, health=0.0)


