import numpy as np


def gen_depth_adj_estimates(meth_matrix, read_mean, read_scale):
    adjusted_meth_values = np.array(meth_matrix.shape)
    for count, row in enumerate(meth_matrix):
        read_numbers = round(np.random.normal(read_mean, read_scale, len(row)))
        adjusted_meth_values[count, :] = np.random.binomial(read_numbers, row) / read_numbers
    return adjusted_meth_values