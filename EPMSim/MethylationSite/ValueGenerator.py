import numpy as np


def generate_sample_methylation(state: float, m_i: float,
                                r_i: float, std_error: float = 0.0):
    error = 0.0 if not std_error else np.random.normal(0, std_error)
    m_hat = m_i + r_i * state
    if m_hat > 1.0:
        m_hat = 1.0
    elif m_hat < 0.0:
        m_hat = 0.0
    if m_hat + error > 1.0:
        error = 0
    elif m_hat + error < 0.0:
        error = 0
    return m_hat + error, error
