# AllenSDK

[![DOI](https://zenodo.org/badge/104984017.svg)](https://zenodo.org/badge/latestdoi/104984017)
[![Twitter](https://img.shields.io/twitter/url/https/twitter.com/bendfulcher.svg?style=social&label=Follow%20%40bendfulcher)](https://twitter.com/bendfulcher)

This repository contains code for:
1. Retrieving gene-expression data from the AllenSDK; and
2. processing it into nice structures for further analysis in Matlab.

Requires Matlab and python.
The [AllenSDK package](http://alleninstitute.github.io/AllenSDK/install.html) for python must be installed.

If anything is unclear or needs improvement, please send questions by [raising an Issue](https://docs.github.com/en/github/managing-your-work-on-github/creating-an-issue) or [sending me an email](mailto:ben.d.fulcher@gmail.com).

This pipeline is based on code developed for [Fulcher and Fornito, _PNAS_ (2016)](https://doi.org/10.1073/pnas.1513302113), and used for [Fulcher et al., _PNAS_ (2019)](https://doi.org/10.1073/pnas.1814144116).
If you find this code useful, consider citing these papers if relevant to your work, or you can cite this code directly using its [DOI](https://doi.org/10.5281/zenodo.3951756).

## Getting a region-by-gene matrix

### Getting full gene information first
This is kind of badly worked up, where you first need to get a full list of genes, by running `AllGenes.py`.

This gives you generic information:
* `geneInfo.csv`
* `geneEntrezID.csv`

### Preparing inputs

The basic workflow is:
1. Get all structure IDs using `WriteStructureInfo.py` (kind of cheat by doing this from Matlab, `WriteStructureIDs.m` -> `structIDs_Oh.csv` and `structInfo_Oh.csv`)
2. Get all gene entrez IDs, e.g., `AllGenes.py` -> `allGenes.csv`
3. Run `RetrieveGene.py` and retrieve the combinations from Allen API

Note that in `RetrieveGene.py`, three variables need to be set.

Input files:
* `structIDSource`: name of the .csv file of Allen structure IDs
* `entrezSource`: name of the .csv file of gene entrez IDs to retrieve

### Preparing outputs

Output files (specified filenames):
* `structInfoFilename`: saves retrieved information for the structure IDs specified
* `allDataFilename`: saves detailed expression information out to this file

Output files (generated filenames):
* `expression_energy_AxB`: expression energy values for the A structures and B section datasets
* `expression_density_AxB`: expression density values for the A structure and B section datasets
* `dataSetIDs_Columns.csv`: dataset IDs representing each column in the above matrices

### Importing python outputs into Matlab

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
