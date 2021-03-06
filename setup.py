#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import glob

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'numpy',
    'pandas',
    'dask',
    'toolz',
    'cloudpickle',
    'xmltodict'
]

test_requirements = [
    # TODO: put package test requirements here
]

scripts = ['scripts/TECAN', 'scripts/LLP-DB', 'scripts/LLP']

long_desc = """
General bioinformatic pipelines associated with the Ley Lab at the MPI in Tuebingen.
For more information, see the README.
"""


setup(
    name='leylab_pipelines',
    version='0.2.1',
    description="General bioinformatic pipelines associated with the Ley Lab at the MPI in Tuebingen",
    long_description=long_desc + '\n\n' + history,
    author="Nick Youngblut",
    author_email='leylabmpi@gmail.com',
    url='https://github.com/leylabmpi/leylab_pipelines',
    packages=find_packages(),
    package_dir={'leylab_pipelines':
                 'leylab_pipelines'},
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='leylab_pipelines',
    classifiers=[
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    scripts=scripts
)
