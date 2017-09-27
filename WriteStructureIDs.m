%-------------------------------------------------------------------------------
% Simple script to output a set of Allen structure IDs:
%-------------------------------------------------------------------------------
C = load('Mouse_Connectivity_Data.mat','RegionStruct');
isCortex = ([C.RegionStruct.OhRegionIndex]==1);
corticalStructureIDs = [C.RegionStruct(isCortex).id];
dlmwrite('structureIDs.csv',corticalStructureIDs);
% fid = fopen('structureIDs.csv','w');
