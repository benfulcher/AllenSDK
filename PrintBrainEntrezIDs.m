function PrintBrainEntrezIDs()
% Idea is to read in and print out entrez_ids from brain
%-------------------------------------------------------------------------------

% First get the relevant brain genes:
[~,~,geneInfo] = GiveMeGeneData(G,'brainGenes','energy');
entrezIDs = geneInfo.entrez_id;
numGenes = length(entrezIDs);
fprintf(1,'Writing %u entrez IDs of brain-related genes to file\n',numGenes);
fid = fopen(fullfile('AllenSDK','brainGeneEntrezID.csv'),'w');
for i = 1:numGenes
    fprintf(fid,'%u\n',entrezIDs(i));
end
fclose(fid)

end
