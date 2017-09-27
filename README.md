# AllenSDK

Requires Matlab, python, and the AllenSDK, [install here](http://alleninstitute.github.io/AllenSDK/install.html).
Questions: ben.d.fulcher@gmail.com; Twitter: @bendfulcher

The basic workflow is:
1. Get all structure IDs (kind of cheat by doing this from Matlab, `WriteStructureIDs.m`) -> `structureIDs.csv`
2. Get all gene entrez IDs, e.g., `AllGenes.py` -> `allGenes.csv`
3. Run `RetrieveGene.py` and retrieve the combinations from Allen API

Then you can import the resulting data into Matlab as:
```matlab
[GeneExpData,sectionDatasetInfo,geneInfo,structInfo] = ImportAllenToMatlab();
```
