function expressionEnergy = getExpressionEnergy(geneEntrezID,structureIDs)
% Retrieves expression energy from the Allen SDK
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Input checking:
%-------------------------------------------------------------------------------
if nargin < 1
    geneEntrezID = 20604;
end
if nargin < 2 || isempty(structureIDs)
    C = load('Mouse_Connectivity_Data.mat','RegionStruct');
    isCortex = ([C.RegionStruct.OhRegionIndex]==1);
    % IDs of all cortical brain regions:
    structureIDs = [C.RegionStruct(isCortex).id];
end
%-------------------------------------------------------------------------------

% Write the gene entrez ID to geneEntrezID.csv
csvwrite(fullfile('AllenSDK','geneEntrezID.csv'),geneEntrezID);

% Write the structure list to structureIDs.csv
csvwrite(fullfile('AllenSDK','structureIDs.csv'),structureIDs);

% Run the python script
% (seems like it's possible to do this through matlab without file i/o, but
% modules need to sit in the python search path?
% cf. https://au.mathworks.com/help/matlab/matlab_external/call-user-defined-custom-module.html)
% 
% P = py.sys.path;
% if count(P,'AllenSDK') == 0
%     insert(P,int32(0),'AllenSDK');
% end

cd('AllenSDK');
system('python RetrieveGene.py');
cd('../')

% Read in the csv output
csvFileName = sprintf('expressionEnergy_gene%u.csv',geneEntrezID);
expressionEnergy = csvread(csvFileName);

% Match output structures to input structures, just to be extra safe:
[~,ia,ib] = intersect(expressionEnergy(:,1),structureIDs);

if length(ia) < length(structureIDs)
    error('Could not retrieve expression data for all structures??');
end

% Output the expression energy vector
expressionEnergy = expressionEnergy(ia,2);

end
