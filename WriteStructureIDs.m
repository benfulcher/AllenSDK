%-------------------------------------------------------------------------------
% Simple script to output a set of Allen structure IDs:

fileName = struct();
fileName.ID = 'structIDs_Oh.csv';
fileName.info = 'structInfo_Oh.csv';

% If you already have a Matlab file containing structure IDs, you can load it from here:
% (e.g., 'AllenGeneDataset_19419.mat')
fprintf(1,'Loading full gene data (NEW: FROM ALLEN SDK)...');
dataFile = 'AllenGeneDataset_19419.mat';
load(dataFile,'structInfo');
fprintf(1,' Done.\n');

% Load data and filter by Isocortex:
numStructs = height(structInfo);

% Write IDs to csv:
dlmwrite(fileName.ID,structInfo.id,'precision','%u');
fprintf(1,'Wrote IDs for %u structures to %s\n',numStructs,fileName.ID);

% Now the full structure info:
index = (1:numStructs)';
acronym = structInfo.acronym;
name = structInfo.name;
id = structInfo.id;
outputStructInfo = table(index,id,acronym,name);
writetable(outputStructInfo,fileName.info);
fprintf(1,'Wrote info for %u structures to %s\n',numStructs,fileName.info);
