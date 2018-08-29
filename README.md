# AllenSDK

Requires Matlab, python, and the AllenSDK, [install here](http://alleninstitute.github.io/AllenSDK/install.html).
Questions: ben.d.fulcher@gmail.com;
Twitter: @bendfulcher

## Getting a region-by-gene matrix

The basic workflow is:
1. Get all structure IDs (kind of cheat by doing this from Matlab, `WriteStructureIDs.m`) -> `structIDs_Oh.csv` and `structInfo_Oh.csv`
2. Get all gene entrez IDs, e.g., `AllGenes.py` -> `allGenes.csv`
3. Run `RetrieveGene.py` and retrieve the combinations from Allen API

Then you can import the resulting data into Matlab as:
```matlab
[GeneExpData,sectionDatasetInfo,geneInfo,structInfo] = ImportAllenToMatlab();
```

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
