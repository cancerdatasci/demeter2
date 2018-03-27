# DEMETER2

DEMETER2 version 1.0

## Installation

Pull the docker image:

```bash
docker pull jmmcfarl/demeter2_pub
```

## Setting up input files

To run DEMETER2 you need to provide the following data as input

* **data files**: matrix of log-fold-change values for each cell line and shRNA
* **data file paths**: list of paths to data files
* **sh_targets file**: file that maps shRNAs to genes 
* **positive control gene list**: list of genes to be used as positive controls
* **negative control gene list**: list of genes to be used as negative controls

Each element is described further below:

### Data files
Each of the data files must be csv file containing a matrix of log-fold-change values, where each row represents a separate shRNA, and each column a separate cell lines. Column headers should contain unique cell line identifiers, and the first column must contain the (21mer) targeting sequences for each shRNA. If you are providing multiple 'batches' of data, cell lines screened in multiple datasets will be matched by the cell line names.

### Data file paths
Create a text file which has the path to each data file you want to include on a separate line, with a single header line "data_file_paths". It should look like this (providing two matrices of LFC values as inputs):

data_file_paths
/data/achilles-55k-batch1-repcollapsed-lfc.csv
/data/achilles-55k-batch2-repcollapsed-lfc.csv
/data/achilles-98k-repcollapsed-lfc.csv

Note: these paths will all be relative paths inside of the docker container, so it's easiest to just put the data files in a folder, and then map that directory to `/data/` inside the container.

### sh_targets
This should be a csv file with the following fields:

* **Barcode Sequence**: This should contain the 21mer targeting sequence of each shRNA in the dataset
* **Gene ID**: This should be the Entrez ID of a gene targeted by the shRNA.
shRNAs that target multiple genes can be represented by multiple rows with the same Barcode Sequence. 
shRNAs that have no gene targets should be represented by "NO_CURRENT" as the gene ID.
Other kinds of gene identifiers could be used, but the values need to be consistent everywhere.

### Positive control genes
This should be a csv file with a field "Gene_ID" that contains the Gene identity of each gene in the positive control set.

### Negative control genes
Same as above for negative control genes.


## Running the model

### Launch the Docker container 
From inside the demeter2 directory run:

```bash
docker run -it -v ./:/demeter2 -v /path/to/data/:/data jmmcfarl/demeter2_pub /bin/bash
```

This will launch the docker container, map the current directory to `/demeter2` inside the container, and your data directory to `/data`, then launch a bash prompt inside the container.

Then, from within the container, you can run the model by e.g.:

```bash
python3 /demeter2/modules/fit_demeter2.py \
	--model_dir output_path \ 
	--data_file /data/data_files_Achilles_paths \
	--sh_targets /data/shRNA-mapping.csv \
	--pos_cons /data/Hart-pos-controls.csv \
	--neg_cons /data/Hart-neg-controls.csv 
```

## Model outputs
* **gene_data**: Matrix of ML estimates of gene effects for each cell line and gene (based on initial alternating ML optimization stage)
* **gene_means**: Matrix of posterior mean estimates of gene effects (from variational Bayesian inference stage)
* **gene_SDs**: Matrix of posterior standard deviation estimates of gene effects (from variational Bayesian inference stage)
* **seed_data**: Matrix of ML seed effect estimates (analogous to gene_data)
* **seed_means**: Matrix of posterior mean seed effects (analogous to gene_means)
* **seed_SDs**: Matrix of posterior std dev seed effects (analogous to gene_SDs)
* **CL_batch_data**: Model parameters inferred per screen (cell line and batch). 
	- CL_offset: additive offset (ML estimate)
	- CL_slope: Overal multiplicative scaling term 
	- noise_vars: Estimated noise variance 
	- offset_mean: posterior mean of CL_offset
	- offset_sd: posterior std dev of CL_offset
* **CL_data**: Data per cell line 
	- 'gene_slope' aka the screen signal term. Tells how the gene KD effects are scaled in each cell line relative to other cell lines
* **hp_batch_data**: Model parameters inferred for each shRNA in each screen
	- hairpin_offset: additive offset term (ML estimate)
	- hairpin_offset_mean: additive offset (posterior mean)
	- hairpin_offset_sd: additive offset (posterior SD)
* **hp_data**: Model parameters inferred globally for each shRNA
	- Geff: gene KD efficacy
	- Seff: seed efficacy
	- unpred_offset: additive offset (shared across screens); ML estimate
	- unpred_offset_mean: posterior mean of above
	- unpred_offset_sd: posterior SD of above
* **per_gene_data**: stores across-cell line average and SD of gene effects for each gene
* **other_info.json**: file that keeps track of some other info from the model fit.
