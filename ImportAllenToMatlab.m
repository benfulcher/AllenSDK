function [geneExpData,sectionDatasetInfo,geneInfo,structInfo] = ImportAllenToMatlab(whatExpressionData)
% ------------------------------------------------------------------------------
% Import the .csv files obtained from python scripts querying the Allen API
% --RetrieveGene.py
%-------------------------------------------------------------------------------
if nargin < 1
    fprintf(1,'Processing layer-specific cortical expression data by default\n');
    whatExpressionData = 'layersCortex'; % 'full', 'layersCortex'
end

% Plot results to figures:
doPlot = false;

%-------------------------------------------------------------------------------
% Define the filenames to import data from:
fileNames = struct();
fileNames.sectionDatasets = 'sectionDatasetInfo.csv';
fileNames.geneInfo = 'geneInfo.csv';

switch whatExpressionData
case 'palomero'
    fileNames.struct = 'structureInfo_Palomero.csv';
    fileNames.energy = 'expression_energy_10x73.csv';
    fileNames.density = 'expression_density_10x73.csv';
    doMatch = false;
case 'layersCortex'
    %---Layer-specific expression:
    fileNames.struct = 'structureInfoLayers.csv';
    fileNames.energy = 'expression_energy_214x295.csv'; % +grik
    fileNames.density = 'expression_density_214x295.csv'; % +grik
    doMatch = true;
case 'full'
    %---Full genome,brain expression:
    fileNames.struct = 'structureInfo.csv';
    fileNames.energy = 'expression_energy_213x25469.csv';
    fileNames.density = 'expression_density_213x25469.csv';
    doMatch = true;
case 'cortexSubset'
    fileNames.struct = 'structInfo_cortex.csv';
    fileNames.energy = 'expression_energy_51x295.csv';
    fileNames.density = 'expression_density_51x295.csv';
    doMatch = true;
case 'cortexAll'
    fileNames.struct = 'structInfo_cortex.csv';
    fileNames.energy = 'expression_energy_51x25460.csv';
    fileNames.density = 'expression_density_51x25460.csv';
    doMatch = true;
otherwise
    error('Unknown expression dataset for %s',whatExpressionData);
end

% ------------------------------------------------------------------------------
%% Assemble region information as a table -> structInfo
% ------------------------------------------------------------------------------
fprintf(1,'Working with structures (from %s)...\n',fileNames.struct);
% Import the data:
structInfo = ImportStructures(fileNames.struct);
numStructures = size(structInfo,1);
numInfo = size(structInfo,2); % variables for each structure
if doMatch
    % Import major region labels assigned by Oh et al.:
    load('Mouse_Connectivity_Data.mat','regionAcronyms','MajorRegionLabels')
    % Match and add to the table:
    [~,~,match_ix] = intersect(structInfo.acronym,regionAcronyms,'stable');
    if length(match_ix)==numStructures
        % All matched:
        structInfo.divisionLabel = MajorRegionLabels(match_ix);
    else
        warning('Could not match regions to Oh et al. major region labels: making all ISOCORTEX!')
        allCortex = cell(height(structInfo),1);
        for i = 1:height(structInfo)
             allCortex{i} = 'Isocortex';
        end
        structInfo.divisionLabel = allCortex;
    end
end

%-------------------------------------------------------------------------------
% Custom (additional) structure filtering:
if ismember(whatExpressionData,{'cortexSubset','cortexAll'})
    [structInfoFilt,ix_ABA40] = StructureFilter(structInfo,'ABAcortex40');
    warning('Filtered from %u to %u cortical structures',height(structInfo),height(structInfoFilt))
    structInfo = structInfoFilt;
    numStructures = size(structInfo,1);
end

%===============================================================================
%% Section datasets and gene information
% Retrieve section dataset info -> sectionDatasetInfo
% Assemble gene info -> geneInfo
% Reads in outputs from AllGenes.py
%===============================================================================

%-------------------------------------------------------------------------------
% Import mapping from section datasets (including plane of section info) to gene expression:
%-------------------------------------------------------------------------------
sectionDatasetInfo = ImportSectionDatasetGeneMapping(fileNames.sectionDatasets);
fprintf(1,'Imported metadata for %u gene expression section datasets from %s\n',...
                        size(sectionDatasetInfo,1),fileNames.sectionDatasets);
fprintf(1,'(%u coronal sections, %u sagittal sections; %u unique genes)\n',...
                sum(sectionDatasetInfo.plane_of_section_id==1),...
                sum(sectionDatasetInfo.plane_of_section_id==2),...
                length(unique(sectionDatasetInfo.entrez_id)));

%-------------------------------------------------------------------------------
% Import gene information:
%-------------------------------------------------------------------------------
% (info about all genes from any of the section datasets)
geneInfo = ImportGeneInformation(fileNames.geneInfo);
fprintf(1,'Imported gene information from Allen for %u genes\n',height(geneInfo));

% Add Cahoy gene types:
minFoldEnrichment = 10;
[GeneCellType,GeneCellTypeName] = CahoyEnrichedGenes(geneInfo.acronym,minFoldEnrichment);
geneInfo.CahoyCellTypeName = GeneCellTypeName;
geneInfo.CahoyCellTypeLabel = GeneCellType(:,1);
geneInfo.CahoyCellTypezscore = GeneCellType(:,2);
fprintf(1,'Added Cahoy et al. cell types for each gene (at %u-fold enrichment threshold)\n',...
                minFoldEnrichment);

%-------------------------------------------------------------------------------
% Import gene expression data
%-------------------------------------------------------------------------------
fprintf(1,'Reading from %s and %s\n',fileNames.energy,fileNames.density);
theFields = {'energy','density'};
geneExpData = struct();
structureIDs = struct();
dataSetIDs = struct();
for k = 1:2
    geneExpData.section.(theFields{k}) = dlmread(fileNames.(theFields{k}),',',1,1);
    [numStructuresRead,numDatasetsRead] = size(geneExpData.section.(theFields{k}));
    structureIDs.(theFields{k}) = dlmread(fileNames.(theFields{k}),',',[1,0,numStructuresRead,0]);
    dataSetIDs.(theFields{k}) = dlmread(fileNames.(theFields{k}),',',[0,1,0,numDatasetsRead]);
    if numStructuresRead~=numStructures
        warning('Taking subsets of gene expression data to match reduced Allen-cortex-40')
        [~,ia,ib] = intersect(structInfo.id,structureIDs.(theFields{k}),'stable');
        geneExpData.section.(theFields{k}) = geneExpData.section.(theFields{k})(ib,:);
        structureIDs.(theFields{k}) = structureIDs.(theFields{k})(ib);
    end
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
allEntrezIDs = unique(sectionDatasetInfo.entrez_id);
numGenes = length(allEntrezIDs);

[~,ia,ib] = intersect(structureIDs.energy,structInfo.id,'stable');
if ~all(ia==(1:length(ia))')
    error('Mismatch');
end
structInfo = structInfo(ib,:);
numStructures = height(structInfo);

fprintf(1,'We have %u genes across %u section datasets, %u structures\n',...
            numGenes,numDatasets,numStructures);

%-------------------------------------------------------------------------------
% So we're there in terms of having the expression data (energy and density)
% as well as tables summarizing details of the datasets and the structures
%-------------------------------------------------------------------------------
% Next is to process the sections that repeat genes
%-------------------------------------------------------------------------------

theFields = {'comb',...% Mean across repeated measures (combining both types of section datasets)
            'combZ',... % Mean across z-score of repeated datasets
            'coronal',... % Take the mean across repeated measures (using only coronal section datasets):
            'sagittal',... % Take the mean across repeated measures (using only sagittal section datasets):
            'replicated',... % Takes mean of (z-scored) section datasets that have been replicated to high correlation
            'benCombo'}; % Takes replicated result if exists, otherwise takes coronal result
% In order to be designated as 'replicated', multiple measurements need to be correlated at least at this level
thresholdReplicated = 0.5;
numFields = length(theFields);
% Initialize geneExpData (energy and density) as NaNs:
for i = 1:numFields
    geneExpData.(theFields{i}).energy = nan(numStructures,numGenes);
    geneExpData.(theFields{i}).density = nan(numStructures,numGenes);
end

numCoronal = zeros(numGenes,1);
numSagittal = zeros(numGenes,1);
isCoronal = (sectionDatasetInfo.plane_of_section_id==1);
isSagittal = (sectionDatasetInfo.plane_of_section_id==2);
for i = 1:numGenes
    isG = (sectionDatasetInfo.entrez_id == allEntrezIDs(i));
    % Count coronal versus sagittal planes:
    numCoronal(i) = sum(isG & isCoronal);
    numSagittal(i) = sum(isG & isSagittal);

    % Combined coronal/sagittal section data ('comb')
    if sum(isG)==1
        % Single section dataset for this gene: easy (no combination required)
        geneExpData.comb.energy(:,i) = geneExpData.section.energy(:,isG);
        geneExpData.comb.density(:,i) = geneExpData.section.density(:,isG);
        geneExpData.combZ.energy(:,i) = BF_NormalizeMatrix(geneExpData.section.energy(:,isG),'zscore');
        geneExpData.combZ.density(:,i) = BF_NormalizeMatrix(geneExpData.section.density(:,isG),'zscore');
        if isCoronal(isG)
            % The only data is from a coronal section:
            geneExpData.coronal.energy(:,i) = geneExpData.section.energy(:,isG);
            geneExpData.coronal.density(:,i) = geneExpData.section.density(:,isG);
            % We can also use this for benCombo:
            geneExpData.benCombo.energy(:,i) = geneExpData.section.energy(:,isG);
            geneExpData.benCombo.density(:,i) = geneExpData.section.density(:,isG);
        else
            if ~isSagittal(isG)
                error('Error labeling coronal/sagittal sections');
            end
            % The only data is from a sagittal section:
            geneExpData.sagittal.energy(:,i) = geneExpData.section.energy(:,isG);
            geneExpData.sagittal.density(:,i) = geneExpData.section.density(:,isG);
        end
    else
        % Multiple section datasets for this gene: need to combine

        % (raw):
        geneExpData.comb.energy(:,i) = nanmean(geneExpData.section.energy(:,isG),2);
        geneExpData.comb.density(:,i) = nanmean(geneExpData.section.density(:,isG),2);

        % (zscored):
        geneExpData.combZ.energy(:,i) = nanmean(BF_NormalizeMatrix(geneExpData.section.energy(:,isG),'zscore'),2);
        geneExpData.combZ.density(:,i) = nanmean(BF_NormalizeMatrix(geneExpData.section.density(:,isG),'zscore'),2);

        % zscore by default:
        if numCoronal(i) > 0
            % Take mean from coronal sections
            geneExpData.coronal.energy(:,i) = nanmean(BF_NormalizeMatrix(geneExpData.section.energy(:,isG & isCoronal),'zscore'),2);
            geneExpData.coronal.density(:,i) = nanmean(BF_NormalizeMatrix(geneExpData.section.density(:,isG & isCoronal),'zscore'),2);
        end
        if numSagittal(i) > 0
            % Take mean from sagittal sections
            geneExpData.sagittal.energy(:,i) = nanmean(BF_NormalizeMatrix(geneExpData.section.energy(:,isG & isSagittal),'zscore'),2);
            geneExpData.sagittal.density(:,i) = nanmean(BF_NormalizeMatrix(geneExpData.section.density(:,isG & isSagittal),'zscore'),2);
        end

        % Now repeated (i.e., there exist multiple section datasets with >0.5 correlation in expression energy):
        r = corr(geneExpData.section.energy(:,isG),'type','Spearman','rows','pairwise');
        r(tril(true(size(r)))) = 0; % exclude lower triangle (+diagonal)
        [ind1,ind2] = find(r > thresholdReplicated);
        if ~isempty(ind1)
            isGindx = find(isG);
            doesReplicate = isGindx(unique([ind1;ind2])); % a bit shit in case of many repeats (e.g., two clusters), but should be ok
            geneExpData.replicated.energy(:,i) = nanmean(BF_NormalizeMatrix(geneExpData.section.energy(:,doesReplicate),'zscore'),2);
            geneExpData.replicated.density(:,i) = nanmean(BF_NormalizeMatrix(geneExpData.section.density(:,doesReplicate),'zscore'),2);
        end

        if ~isempty(ind1)
            geneExpData.benCombo.energy(:,i) = geneExpData.replicated.energy(:,i);
            geneExpData.benCombo.density(:,i) = geneExpData.replicated.density(:,i);
        elseif numCoronal(i) > 0
            geneExpData.benCombo.energy(:,i) = geneExpData.coronal.energy(:,i);
            geneExpData.benCombo.density(:,i) = geneExpData.coronal.density(:,i);
        end
    end
end
fprintf(1,'Took mean for %u genes with multiple section datasets\n',...
                    sum(numCoronal + numSagittal > 1));

% Make a gene info table that matches these
[~,ia,ib] = intersect(allEntrezIDs,geneInfo.entrez_id);
if ~all(ia==(1:length(ia))')
    error('Mismatch');
end
geneInfo = geneInfo(ib,:);
geneInfo.numCoronal = numCoronal;
geneInfo.numSagittal = numSagittal;

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
    theFields = fieldnames(geneExpData);
    numFields = length(theFields);
    for i = 1:numFields
        fprintf(1,'Matches structure ordering for %s\n',theFields{i});
        geneExpData.(theFields{i}).energy = geneExpData.(theFields{i}).energy(ib,:);
        geneExpData.(theFields{i}).density = geneExpData.(theFields{i}).density(ib,:);
    end
end

%-------------------------------------------------------------------------------
% Save
%-------------------------------------------------------------------------------
fileNames.output = sprintf('AllenGeneDataset_%u_%u.mat',numStructures,numGenes);
fileNames.output = fullfile('Data',fileNames.output);
save(fileNames.output,'geneExpData','sectionDatasetInfo','geneInfo','structInfo');
fprintf(1,'Saved processed Allen data to ''%s''!\n',fileNames.output);

end
