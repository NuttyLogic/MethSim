from typing import Tuple
import numpy as np


class SiteGenerator:

    def __init__(self, state_range: float = 1.0,
                 max_delta: float = 0.4,
                 site_deviation: float = 0.1):
        self.state_range = abs(float(state_range))
        self.max_delta = abs(float(max_delta))
        self.site_deviation = abs(float(site_deviation))

    def simulate_site(self) -> Tuple[float, float, float]:
        # pick random initial methylation value
        m_not_i = np.random.beta(.15, .15)
        # sample std deviation for the site
        std = np.random.uniform(0, self.site_deviation)
        pos = m_not_i < .6
        if .6 >= m_not_i > .4:
            pos = np.random.rand() > .5
        if pos:
            m_target = np.random.uniform(m_not_i + self.max_delta, 1)
        else:
            m_target = np.random.uniform(0, m_not_i - self.max_delta)
        # set site rate
        try:
            rate = (m_target - m_not_i) / self.state_range
        except ZeroDivisionError:
            rate = 0
        return m_not_i, rate, std
