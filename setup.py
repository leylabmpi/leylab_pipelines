#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
import glob

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'pandas'
]

test_requirements = [
    # TODO: put package test requirements here
]

scripts = glob.glob('scripts/*.py')

setup(
    name='leylab_pipelines',
    version='0.1.0',
    description=" General bioinformatic pipelines associated with the Ley Lab at the MPI in Tuebingen",
    long_description=readme + '\n\n' + history,
    author="Nick Youngblut",
    author_email='nyoungb2@gmail.com',
    url='https://github.com/nick-youngblut/leylab_pipelines',
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
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
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
