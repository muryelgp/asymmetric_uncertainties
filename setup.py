from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='asymmetric_uncertainties',
    version='1.0.0',
    description='A package for dealing and propagating asymmetric uncertainties in measurements using Monte Carlo simulations.',
    long_description=readme,
    author='Muryel Guolo Pereira',
    author_email='muryel@astro.ufsc.br',
    url='https://github.com/muryelgp/asymmetric_uncertainties',
    download_url = 'https://github.com/muryelgp/asymmetric_uncertainties/archive/master.zip',
    license=license,
    keywords = ['statistics', 'uncertanties', 'propagation'],
    packages=['asymmetric_uncertainties'],
    install_requires=[
          'numpy',
          'matplotlib',
          'scipy',
      ],
)
