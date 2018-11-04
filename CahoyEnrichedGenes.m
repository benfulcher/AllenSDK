function [geneCellType,geneCellTypeName] = CahoyEnrichedGenes(ABA_GeneAcronyms,minEnrichment)
% Imports excel files from Cahoy et al. (J. Neurosci., 2008) for genes enriched in
% astrocytes, neurons, and ogligodendrocytes.
%
% Excel files are Supplements S4--S6 from the paper:
% http://www.jneurosci.org/content/28/1/264/suppl/DC1
%
% cf. Analysis in French et al. (front. Neurosci., 2011)
%-------------------------------------------------------------------------------

if nargin < 2
    % Minimum enrichment level to label a gene according to a given cell type:
    minEnrichment = 10; % (10 by default, as used in Tan et al., Front. Neurosci., 2013)
end
% ------------------------------------------------------------------------------


% Transpose ABA_GeneAcronyms if a column vector:
if size(ABA_GeneAcronyms,1) > size(ABA_GeneAcronyms,2)
    ABA_GeneAcronyms = ABA_GeneAcronyms';
end

% Set file names for excel files to import data from:
fileNames = {'340_Cahoy_S_Table_S4_AvX_2007-09-11.xls', ...
            '350_Cahoy_S_Table_S5_OvX_2007-09-11.xls', ...
            '360_Cahoy_S_Table_S6_NvX_2007-09-11.xls'};

% Genes statistically enriched in the following cell types:
CellType = {'astrocyte','ogligodendrocyte','neuron'};

% Genes from ABA
% GeneStruct = G.GeneStruct;
numGenes = length(ABA_GeneAcronyms);
% ABA_GeneAcronyms = {GeneStruct.gene_acronym};
geneCellType = zeros(numGenes,2); % first column: type; second column: fold enriched

% In format RowNum,Probe_Set_ID,Gene_Name,Fold_Enriched
% We want to extract the list of gene names, then to match them to genes from the Allen Brain Atlas

for k = 1:3
    fprintf(1,'Loading from %s...',fileNames{k});
    [AllNum,AllText] = xlsread(fileNames{k});
    fprintf(1,' Done.\n');
    GeneNames_k = AllText(3:end,3);
    FoldEnriched_k = AllNum(:,4);

    % See if any ABA acronyms match names in this list:
    isCellType_k = cellfun(@(x)ismember(x,GeneNames_k),ABA_GeneAcronyms);

    % Generate matrix the same size ABA, that puts fold enriched for each gene implicated:
    theFold = zeros(size(isCellType_k));
    for i = find(isCellType_k)
        here = find(strcmp(ABA_GeneAcronyms{i},GeneNames_k));
        theFold(i) = FoldEnriched_k(here);
    end

    % Threshold the fold level according to minEnrichment:
    excludeThese = (theFold < minEnrichment);
    theFold(excludeThese) = 0;
    isCellType_k(excludeThese) = 0;

    % Check not already assigned:
    if any(geneCellType(isCellType_k,1))~=0

        % Distinguish overlapping, and non-overlapping genes:
        f_overlap = find(geneCellType(:,1) ~= 0 & isCellType_k');
        f_nooverlap = find(geneCellType(:,1) == 0 & isCellType_k');

        % First assign non-overlapping genes
        geneCellType(f_nooverlap,1) = k;
        geneCellType(f_nooverlap,2) = theFold(f_nooverlap);

        for i = 1:length(f_overlap) % each gene with an overlapping classification:
            fprintf(1,'Multiple cell type assignment: %s, %s, fold = %.2f; %s, fold = %.2f\n', ...
                                                ABA_GeneAcronyms{f_overlap(i)},...
                                                CellType{geneCellType(f_overlap(i),1)}, ...
                                                geneCellType(f_overlap(i),2), ...
                                                CellType{k}, ...
                                                theFold(f_overlap(i)));
            % Assign to the higher fold:
            if geneCellType(f_overlap(i),2) < theFold(f_overlap(i))
                geneCellType(f_overlap(i),1) = k;
                geneCellType(f_overlap(i),2) = theFold(f_overlap(i));
            end
        end

    else
        % No overlap, just assign all to the current type
        geneCellType(isCellType_k,1) = k;
        geneCellType(isCellType_k,2) = theFold(isCellType_k);
    end

    fprintf(1,'Found %u genes in the ABA associated with %s\n\n',sum(isCellType_k),CellType{k});
end


% Define the CahoyCellTypeLabel:
CahoyCellTypeLabel = geneCellType(:,1);

fprintf(1,'We labeled: %u as astrocyte, %u as ogligodendrocyte, %u as neuron, %u as other\n', ...
                            sum(CahoyCellTypeLabel==1), ...
                            sum(CahoyCellTypeLabel==2), ...
                            sum(CahoyCellTypeLabel==3), ...
                            sum(CahoyCellTypeLabel==0));

% Convert labels to strings
geneCellTypeName = cell(numGenes,1);
for k = 0:3
    IsTypek = find(geneCellType(:,1)==k);
    for l = 1:length(IsTypek)
        if k==0
            geneCellTypeName{IsTypek(l)} = 'other';
        else
            geneCellTypeName{IsTypek(l)} = CellType{k};
        end
    end
end

end
