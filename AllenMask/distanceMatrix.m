coOrds = csvread('centre.csv');
fid = fopen('acronym.csv','r');
regAcronyms = textscan(fid,'%s');
fclose(fid);
regAcronyms = regAcronyms{1};

%-------------------------------------------------------------------------------
dataFile = '/Users/benfulcher/GoogleDrive/Work/CurrentProjects/CellTypesMouse/Code/AllenGeneDataset_19419.mat';
fprintf(1,'New Allen SDK-data from %s\n',dataFile);
load(dataFile,'structInfo');

% Match:
[~,ia,ib] = intersect(structInfo.acronym,regAcronyms,'stable');
structInfo = structInfo(ia,:);
coOrds = coOrds(ib,:);
numRegions = height(structInfo);

% coOrds(120,:) = [];
% structInfo(120,:) = [];
% numRegions = numRegions-1;

%-------------------------------------------------------------------------------
% Get Euclidean distances and rescale to 2d
%-------------------------------------------------------------------------------
distMat = squareform(pdist(coOrds,'euclidean'));
score = mdscale(distMat,2);
xData = score(:,1);
yData = score(:,2);

%-------------------------------------------------------------------------------
% Plot
%-------------------------------------------------------------------------------
f = figure('color','w');
dotColors = arrayfun(@(x)rgbconv(structInfo.color_hex_triplet{x})',...
                                        1:numRegions,'UniformOutput',0);
dotColors = [dotColors{:}]';

nodeSize = 50;
% scatter(xData,yData,zData,nodeSize,dotColors,'fill','MarkerEdgeColor','k')
scatter3(coOrds(:,1),coOrds(:,2),coOrds(:,3),nodeSize,dotColors,'fill','MarkerEdgeColor','k')

% Add labels:
xDataRange = range(xData);
for i = 1:numRegions
    text(xData(i)+0.04*xDataRange,yData(i),structInfo.acronym{i},...
                        'color',brighten(dotColors(i,:),-0.3))
end

% And for major divisions:
divisionLabels = categorical(structInfo.divisionLabel);
theDivisions = unique(divisionLabels);
numDivisions = length(theDivisions);
for i = 1:numDivisions
    % Put each major region in the center of those points
    centrePoint = [mean(xData(divisionLabels==theDivisions(i))),mean(yData(divisionLabels==theDivisions(i)))];
    find_1 = find(divisionLabels==theDivisions(i),1);
    text(centrePoint(1),centrePoint(2),char(theDivisions(i)), ...
                'color','k','FontSize',14,'BackgroundColor',dotColors(find_1,:))
end
