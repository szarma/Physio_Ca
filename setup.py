from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open('README.md') as f:
    long_description = f.read()

with open('version') as f:
    version = f.read().strip()

setup(
    name='islets',
    version=version,
    description='Module for analysing islets.',
    packages=find_packages(),
    long_description=long_description,
    install_requires=requirements,
    author='Srdjan Sarikas & Johannes Pfabe'
)
