# AllenSDK

[![DOI](https://zenodo.org/badge/104984017.svg)](https://zenodo.org/badge/latestdoi/104984017)
[![Twitter](https://img.shields.io/twitter/url/https/twitter.com/bendfulcher.svg?style=social&label=Follow%20%40bendfulcher)](https://twitter.com/bendfulcher)

This repository contains code for:
1. Retrieving gene-expression data from the AllenSDK; and
2. processing it into nice structures for further analysis in Matlab.

Requires Matlab and python.
The [AllenSDK package](http://alleninstitute.github.io/AllenSDK/install.html) for python must be installed.

Commands to install AllenSDK (As incompatible with newer versions)
```
conda create -n allensdk python=3.7
conda activate allensdk
pip install allensdk
```

If anything is unclear or needs improvement, please send questions by [raising an Issue](https://docs.github.com/en/github/managing-your-work-on-github/creating-an-issue) or [sending me an email](mailto:ben.d.fulcher@gmail.com).

This pipeline is based on code developed for [Fulcher and Fornito, _PNAS_ (2016)](https://doi.org/10.1073/pnas.1513302113), and used for [Fulcher et al., _PNAS_ (2019)](https://doi.org/10.1073/pnas.1814144116).
If you find this code useful, consider citing these papers if relevant to your work, or you can cite this code directly using its [DOI](https://doi.org/10.5281/zenodo.3951756).

## Constructing a brain region x gene matrix

### Retrieve full gene information
You first need to get a full list of genes, by running `AllGenes.py`.

This outputs you generic information about the genes:
* `sectionDatasetInfo.csv` (all section data)
* `geneInfo.csv` (gene information: acronym, entrez_id, gene_id, name)
* `geneEntrezID.csv` (just the list of EntrezIDs)

### Preparing inputs for a specific region x gene matrix

#### 1. Retrieve IDs for all brain regions, `structIDs` and `structInfo`

Retrieve all structure IDs of interest directly by adapting `WriteStructureInfo.py` to retrieve a custom set of structures.

If you already have structure IDs in Matlab, you can alternatively to this step using `WriteStructureIDs.m` -> `structIDs_Oh.csv` and `structInfo_Oh.csv`.

#### 2. Retrieve gene entrez IDs

Save a list of gene entrez IDs for the genes you're interested in.
For all genes, you can use the `geneEntrezID.csv` file produced from `AllGenes.py` above.
For a subset of genes, you can adapt something like `subsetGenes.py`.

#### 3. Run retrieve the expression data from the Allen API

Now you've defined the structures and genes you're interested in, you can run the queries to get all combinations of expression data (of brain regions and genes).
This is done using `RetrieveGene.py`.

Note that in `RetrieveGene.py`, variables need to be set.

First the input files need to match the IDs saved in Steps 1 and 2 above.

___Input files___
* `structIDSource`: name of the `.csv` file of Allen structure IDs
* `entrezSource`: name of the `.csv` file of gene entrez IDs to retrieve

___Output filenames___

__To set:__
* `structInfoFilename`: saves retrieved information for the structure IDs specified.
* `allDataFilename`: saves detailed expression information out to this file.

__Generated:__
* `expression_energy_AxB`: expression energy values for the A structures and B section datasets
* `expression_density_AxB`: expression density values for the A structure and B section datasets
* `dataSetIDs_Columns.csv`: dataset IDs representing each column in the above matrices

## Importing data into Matlab

Then you can import the resulting data into Matlab as:
```matlab
[GeneExpData,sectionDatasetInfo,geneInfo,structInfo] = ImportAllenToMatlab();
```

In this function, you must specify the filenames to read in:
* `fileNames.struct`: the structure info file specified above (`structInfoFilename`)
* `fileNames.sectionDatasets`: full information about all datasets retrieved (`allDataFilename`)
* `fileNames.geneInfo`:
* `fileNames.energy`:
* `fileNames.density`:
* `fileNames.columns`:

Outputs a processed .mat file: `AllenGeneDataset_X.mat` containing information about X unique genes.

## Computing a structure mask

Example pipeline:
First generate `.csv` files for structure IDs and matching to structure info (for interpretation)
E.g., for the Oh et al. 213-region parcellation:
```matlab
WriteStructureIDs
```
This generates `structIDs_Oh.csv` and `structInfo_Oh.csv`.
In the python file `MakeCCFMasks`, these files are listed as inputs, such that
```python
MakeCCFMasks
```
generates a mask for these, saving as `mask_Oh.h5`.
