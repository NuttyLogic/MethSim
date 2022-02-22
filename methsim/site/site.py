import numpy as np


def simulate_site(site_deviation=0.0, delta=1.0, beta_params=(.3, .3)):
    # pick random initial methylation value
    m_not_i = np.random.beta(*beta_params)
    # sample std deviation for the site
    std = np.random.uniform(0, site_deviation)
    pos = m_not_i < .3
    if .3 <= m_not_i < .7:
        pos = np.random.rand() > .5
    if pos:
        target = m_not_i + delta
        target = 1.0 if target > 1.0 else target
    else:
        target = m_not_i - delta
        target = 0.0 if target < 0.0 else target
    return np.array([m_not_i, target, std])


def simulate_site_matrix(deviation_low=0.001, deviation_high=0.1,
                         delta_low=0.0, delta_high=0.5, number_of_sites=5000,
                         beta_params=(.3, .3)):
    sites = np.zeros((number_of_sites, 3))
    for site in range(number_of_sites):
        site_dev = np.random.uniform(deviation_low, deviation_high)
        site_delta = np.random.uniform(delta_low, delta_high)
        sites[site, :] = simulate_site(site_deviation=site_dev, delta=site_delta,
                                       beta_params=beta_params)
    return sites


def normalize_weights(state_weights=None):
    """Normalize weights for coef fitting"""
    adjusted_weights = np.array(state_weights) / np.sum(abs(np.array(state_weights)))
    if np.sum(adjusted_weights) < 0.5:
        return abs(adjusted_weights)
    return adjusted_weights


def generate_site_coefs(m_i: float, m_target: float, state_weights: np.ndarray = np.ones(1)):
    weights = normalize_weights(state_weights)
    site_delta = m_target - m_i
    return site_delta * weights


def generate_sample_methylation(states: np.ndarray, m_i: float,
                                std_error: float = 0.0,
                                coefs: np.ndarray = np.ones((1, 1))):
    m_hat = m_i + np.dot(states, coefs)
    # add random error
    error = np.random.normal(0, std_error, m_hat.shape[0])
    m_hat[m_hat > 1.0] = 1.0
    m_hat[m_hat < 0.0] = 0.0
    error[m_hat + error > 1.0] = 0.0
    error[m_hat + error < 0.0] = 0.0
    return m_hat + error, error
