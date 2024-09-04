#!/bin/env python

from setuptools import setup, find_packages

setup(
    name='genomic-mlp',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'tensorflow',
        'pandas',
        'pandas-plink',
        'keras-tuner',
        'dask',
        'numpy',
        'scikit-learn',
        'matplotlib'
    ],
    entry_points={
        'console_scripts': [
            'genomic-mlp=main:main',
        ],
    },
)
