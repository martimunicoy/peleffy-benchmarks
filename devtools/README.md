# Development tools
This directory contains a collection of tools for developers.

## Conda environments
A standard conda environment to run the benchmarks is found in the `conda-envs`folder.
This environment can be set up by executing the commands from below:
```
conda config --add channels conda-forge --add channels omnia --add channels martimunicoy
conda update --all
conda env create -f devtools/conda-envs/environment.yaml
```
