# coinsparse: A simple command line tool to parse phenotypic data from the Healthy Brain Network

This repository contains simple command line tools to summarize phenotypic data
from the Healthy Brain Network (HBN) study. To gain access to the HBN phenotypic data visit [the HBN Pheno Access page](http://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network/Pheno_Access.html). This repository is a work in progress.

## Building this software

### Clone and pip install

To use this option, clone this repository and run
```bash
python -m pip install .
```
from the root directory.

### Docker installation

To use this option, you must have [Docker
installed](https://docs.docker.com/get-docker/). After that, you can install
the software by cloning this repository to your local computer and using
```bash
docker compose build
```
from the root directory of the repository.

## Usage

For usage information on the command line tool type
```bash
docker compose run coinsparse --help
```

All of the subcommands and parameters discussed therein can be appended to the command
```bash
docker compose run coinsparse ...
```
without the ellipses.

For example, to summarize only the "FSQ" and "WIAT" assessments, with data dictionaries located in the `data/data_dictionaries_coins` directory and phenotypic data located in the `data/coins_assessment_data_20210329` directory, and computing the intersection of subjects with tractometry data data located in `data/combined_tract_profiles.csv`, you would use the command
```bash
docker compose run coinsparse summarize -a FSQ -a WIAT --deriv-file data/combined_tract_profiles.csv --deriv-name tractometry data/data_dictionaries_coins data/coins_assessment_data_20210329
```

Or to list all of the subjects with the aforesaid tractometry data and responses to the "FSQ_08" and "WIAT_RC_Stnd" variables, you would use the command
```bash
docker compose run coinsparse list-subjects -v FSQ_08 -v WIAT_RC_Stnd -d data/combined_tract_profiles.csv data/data_dictionaries_coins data/coins_assessment_data_20210329
```

Or perhaps you want to subset the phenotypic and derivatives data with only the subjects output by the previous command. Then you'd use
```bash
docker compose run coinsparse merge
```
with the same arguments as the last command.