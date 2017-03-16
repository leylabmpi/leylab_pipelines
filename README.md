[![Build Status](https://travis-ci.org/leylabmpi/leylab_pipelines.svg?branch=master)](https://travis-ci.org/leylabmpi/leylab_pipelines)
[![PyPI version](https://badge.fury.io/py/leylab_pipelines.svg)](http://badge.fury.io/py/leylab_pipelines)

leylab_pipelines
=================

 General bioinformatic pipelines associated with the Ley Lab at the MPI in Tuebingen


#### Sections

- [Contents](#contents)
- [Examples](#examples)
- [Installation](#installation)
- [Usage](#usage)
- [Changelog](#changelog)
- [License](#license)


## Contents

[[top](#sections)]

### TECAN robot helper scripts

#### `map2robot.py`

Convert QIIME formatted mapping file to NGS library prep commands for the TECAN robot. 

#### `dilute.py`

Based on a table of sample concentrations, create TECAN robot commands for sample dilution.

#### `qPCR_setup.py`

Based on the plate layout of a qPCR experiment, create TECAN robot commands to set up the 
PCR reactions by aliquoting samples and reagents into the final PCR plate. 

## Examples

[[top](#sections)]

### TECAN robot helper scripts

* See the bash scripts in:
  * `./examples/map2robot/`
  * `./examples/dilute/`
  * `./examples/qPCR_setup/`


## Installation

[[top](#sections)]

### Testing

`python setpy.py test`

### Build

`python setup.py install`


## Changelog

[[top](#sections)]


# License

[[top](#sections)]

* Free software: MIT license
* Documentation: https://leylab-pipelines.readthedocs.io.
