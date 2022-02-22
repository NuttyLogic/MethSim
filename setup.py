
try:
    from setuptools import setup
except ImportError as e:
    print('Please install python setuptools, '
          'https://packaging.python.org/tutorials/installing-packages/#use-pip-for-installing ')
    raise e

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='methsim',
      version='0.0.2',
      description='Utility for Simulating Epigenetic Aging Matrices',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/NuttyLogic/EpigeneticPacemakerSimulation',
      author='Colin P. Farrell',
      author_email='colinpfarrell@gmail.com',
      license='MIT',
      packages=['methsim',
                'methsim.sample',
                'methsim.site',
                'methsim.phenotype',
                'methsim.utilities'],
      classifiers=['Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   'Programming Language :: Python :: 3.8'],
      platforms=["Linux", "Mac OS-X", "Unix"],
      requires=['numpy', 'sklearn'],
      install_requires=['numpy>=1.16.3', 'setuptools>=46.0.0'],
      python_requires='>=3.6',
      test_suite='tests',
      # include_package_data=True,
      zip_safe=True
      )
