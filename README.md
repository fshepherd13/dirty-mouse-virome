# Natural rodent model for the study of virus transmission reveals diverse factors affect virus populations
Elizabeth J. Fay, Keir M. Balla, Shanley Roach, Frances Shepherd, Dira Putri, Talia Wiggen, Stephen A. Goldstein, Mark Pierson, Andrew Tucker, Sergei V. Kotenko, Ryan C. Hunter, David Masopust, Nels C. Elde, and Ryan A. Langlois

### Getting started

This analysis requires several dependencies which can be installed as a conda environment using the environment.yaml file.
```
$ conda install --file environment.yaml
$ conda activate dirty_mouse
```

The pipelines are based off of code originally written by Karthik Gangavarapu in 2018 (https://github.com/andersen-lab/ivar). See the Andersen Lab's github page for further details on ivar, which is code written to call variants in amplicon sequencing reads. Snakemake is used to run the pipelines. For details on how to use snakemake, see https://snakemake.readthedocs.io/en/stable/.

### Further details
The `analysis` folder contains two types of pipelines, one for calling variants as viruses transmit from a reservoir hosts(i.e. pet store mice) to recipient hosts (i.e. SPF lab mice), called `transmission_analysis`. The second directory, named `dissemination_analysis` is used for calling variants as viruses disseminate from small intestine to the liver. They work slightly differently, and there are readme files within each directory to explain how to run the pipeline. 

### Raw data availability
Raw sequencing data is deposited in NCBI under BioProject ID: 
