#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
import glob

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'pandas'
]

test_requirements = [
    # TODO: put package test requirements here
]

#scripts = glob.glob('scripts/TECAN')
scripts = ['scripts/TECAN', 'scripts/DB']

long_desc = """
General bioinformatic pipelines associated with the Ley Lab at the MPI in Tuebingen.
For more information, see the README.
"""


setup(
    name='leylab_pipelines',
    version='0.1.0',
    description="General bioinformatic pipelines associated with the Ley Lab at the MPI in Tuebingen",
    long_description=long_desc + '\n\n' + history,
    author="Nick Youngblut",
    author_email='nyoungb2@gmail.com',
    url='https://github.com/leylabmpi/leylab_pipelines',
    packages=[
        'leylab_pipelines',
    ],
    package_dir={'leylab_pipelines':
                 'leylab_pipelines'},
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='leylab_pipelines',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
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
