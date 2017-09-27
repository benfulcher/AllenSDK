function structureInfo = ImportStructures(fileName)
%% Import data from structure csv file output from python

%% Initialize variables.
if nargin < 1
    fileName = 'structureInfo.csv';
end

delimiter = ',';
startRow = 2;

%% Format for each line of text:
%   column2: text (%q)
%	column3: text (%q)
%   column4: double (%f)
%	column5: text (%q)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*q%q%q%f%q%[^\n\r]';

%% Open the text file.
fileID = fopen(fileName,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
structureInfo = table(dataArray{1:end-1}, 'VariableNames', {'acronym','color_hex_triplet','id','name'});
