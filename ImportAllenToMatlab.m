function [GeneExpData,sectionDatasetInfo,geneInfo,structInfo] = ImportAllenToMatlab()
% ------------------------------------------------------------------------------
% Import the .csv files obtained from python scripts querying the Allen API
% --RetrieveGene.py
%-------------------------------------------------------------------------------

% Plot results to figures:
doPlot = false;

% Layer-specific expression:
fileNames = struct();
fileNames.struct = 'structureInfoLayers.csv';
fileNames.sectionDatasets = 'sectionDatasetInfo.csv';
fileNames.geneInfo = 'geneInfo.csv';
fileNames.energy = 'expression_energy_204x280.csv';
fileNames.density = 'expression_density_204x280.csv';
fileNames.columns = 'dataSetIDs_Columns.csv';
% fileNames.struct = 'structureInfo.csv';
% fileNames.sectionDatasets = 'sectionDatasetInfo.csv';
% fileNames.geneInfo = 'geneInfo.csv';
% fileNames.energy = 'expression_energy_213x25469.csv';
% fileNames.density = 'expression_density_213x25469.csv';

% ------------------------------------------------------------------------------
%% Assemble region information as a table -> structInfo
% ------------------------------------------------------------------------------
fprintf(1,'Working with structures (from %s)...\n',fileNames.struct);
% Import the data:
structInfo = ImportStructures(fileNames.struct);
numStructures = size(structInfo,1);
numInfo = size(structInfo,2); % variables for each structure
% Import major region labels assigned by Oh et al.:
load('Mouse_Connectivity_Data.mat','regionAcronyms','MajorRegionLabels')
% Match and add to the table:
[~,~,match_ix] = intersect(structInfo.acronym,regionAcronyms,'stable');
if length(match_ix)==numStructures
    % All matched:
    structInfo.divisionLabel = MajorRegionLabels(match_ix);
else
    warning('Could not match regions to Oh et al. major region labels')
end

% ------------------------------------------------------------------------------
%% Section datasets and gene information
% Retrieve section dataset info -> sectionDatasetInfo
% Assemble gene info -> geneInfo
% Reads in outputs from AllGenes.py
% ------------------------------------------------------------------------------

% Import mapping from section datasets (including plane of section info) to gene expression:
sectionDatasetInfo = ImportSectionDatasetGeneMapping(fileNames.sectionDatasets);

fprintf(1,'Imported metadata for %u gene expression section datasets from %s\n',...
                        size(sectionDatasetInfo,1),fileNames.sectionDatasets);

% Import gene information:
% (info about all genes from any of the section datasets)
geneInfo = ImportGeneInformation(fileNames.geneInfo);

% Add Cahoy gene types:
minFoldEnrichment = 10;
[GeneCellType,GeneCellTypeName] = CahoyEnrichedGenes(geneInfo.acronym,minFoldEnrichment);
geneInfo.CahoyCellTypeName = GeneCellTypeName;
geneInfo.CahoyCellTypeLabel = GeneCellType(:,1);
geneInfo.CahoyCellTypezscore = GeneCellType(:,2);
fprintf(1,'Added Cahoy et al. cell types for each gene\n');

%-------------------------------------------------------------------------------
% Import gene expression data
%-------------------------------------------------------------------------------
fprintf(1,'Reading from %s and %s\n',fileNames.energy,fileNames.density);
theFields = {'energy','density'};
GeneExpData = struct();
structureIDs = struct();
dataSetIDs = struct();
for k = 1:2
    GeneExpData.(theFields{k}) = dlmread(fileNames.(theFields{k}),',',1,1);
    [numStructuresRead,numDatasetsRead] = size(GeneExpData.(theFields{k}));
    structureIDs.(theFields{k}) = dlmread(fileNames.(theFields{k}),',',[1,0,numStructuresRead,0]);
    dataSetIDs.(theFields{k}) = dlmread(fileNames.(theFields{k}),',',[0,1,0,numDatasetsRead]);
end

if ~all(structureIDs.energy==structureIDs.density) || ~all(dataSetIDs.energy==dataSetIDs.density)
    error('Data for expression energy and density don''t match...');
end
if ~all(ismember(dataSetIDs.energy,sectionDatasetInfo.section_id))
    error('More sections in expression matrix (%u) then sectiondataset info (%u) :-[',...
                    length(datasetIDs.energy),numDatasets);
end
if ~all(ismember(structureIDs.energy,structInfo.id))
    error('Structures in expression matrix (%u) don''t match imported structure info (%u)',...
                    length(structureIDs.energy),numStructures);
end

% Gene data (columns) sorted by dataset IDs, and structure data (rows) sorted
% by structure ID
fprintf(1,'Data Imported!\n');

%-------------------------------------------------------------------------------
% Match so we have tables with info for rows (structures) and columns (sections)
%-------------------------------------------------------------------------------
[~,ia,ib] = intersect(dataSetIDs.energy,sectionDatasetInfo.section_id,'stable');
if ~all(ia==(1:length(ia))')
    error('Mismatch');
end
sectionDatasetInfo = sectionDatasetInfo(ib,:);
numDatasets = height(sectionDatasetInfo);

[~,ia,ib] = intersect(structureIDs.energy,structInfo.id,'stable');
if ~all(ia==(1:length(ia))')
    error('Mismatch');
end
structInfo = structInfo(ib,:);
numStructures = height(structInfo);

%-------------------------------------------------------------------------------
% So we're there in terms of having the expression data (energy and density)
% as well as tables summarizing details of the datasets and the structures
%-------------------------------------------------------------------------------
% Next is to process the sections that repeat genes
%-------------------------------------------------------------------------------
allEntrezIDs = unique(sectionDatasetInfo.entrez_id);
numGenes = length(allEntrezIDs);
fprintf(1,'We have %u genes across %u section datasets\n',...
                length(allEntrezIDs),height(sectionDatasetInfo));
% Let's take the mean across repeated measures
GeneExpData.gene_energy = zeros(numStructures,numGenes);
GeneExpData.gene_density = zeros(numStructures,numGenes);
didMean = false(numGenes,1);
for i = 1:numGenes
    isG = (sectionDatasetInfo.entrez_id == allEntrezIDs(i));
    if sum(isG)==1
        GeneExpData.gene_energy(:,i) = GeneExpData.energy(:,isG);
        GeneExpData.gene_density(:,i) = GeneExpData.density(:,isG);
    else
        GeneExpData.gene_energy(:,i) = nanmean(GeneExpData.energy(:,isG),2);
        GeneExpData.gene_density(:,i) = nanmean(GeneExpData.density(:,isG),2);
        didMean(i) = true;
    end
end
fprintf(1,'Took mean for %u genes with multiple section datasets\n',sum(didMean));

% Make a gene info table that matches these
[~,ia,ib] = intersect(allEntrezIDs,geneInfo.entrez_id);
if ~all(ia==(1:length(ia))')
    error('Mismatch');
end
geneInfo = geneInfo(ib,:);

%-------------------------------------------------------------------------------
% If the 213 regions, match to the Oh et al. ordering
%-------------------------------------------------------------------------------
if numStructures==213
    fprintf(1,'Matching to connectivity data for 213 regions\n');
    C = load('Mouse_Connectivity_Data.mat','RegionStruct');
    regionIDsConnectivity = [C.RegionStruct.id]';
    % Match to the connectivity region ordering
    [~,ia,ib] = intersect(regionIDsConnectivity,structInfo.id,'stable');
    if ~all(ia==(1:length(ia))') || ~(length(ia)==length(ib))
        error('Mismatch');
    end
    % Now reorder to match the connectivity ordering:
    structInfo = structInfo(ib,:);
    GeneExpData.gene_energy = GeneExpData.gene_energy(ib,:);
    GeneExpData.gene_density = GeneExpData.gene_density(ib,:);
    GeneExpData.energy = GeneExpData.energy(ib,:);
    GeneExpData.density = GeneExpData.density(ib,:);
end

%-------------------------------------------------------------------------------
% Save
%-------------------------------------------------------------------------------
fileNames.output = sprintf('AllenGeneDataset_%u.mat',numGenes);
save(fileNames.output,'GeneExpData','sectionDatasetInfo','geneInfo','structInfo');
fprintf(1,'Saved %s\n',fileNames.output);
