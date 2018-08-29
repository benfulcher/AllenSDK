function geneInfo = ImportGeneInformation(fileName, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   GENEINFO = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
%   GENEINFO = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   geneInfo = importfile('geneInfo.csv', 2, 25536);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/11/16 19:29:47

%% Initialize variables.
delimiter = ',';
if nargin < 1
    fileName = 'geneInfo.csv';
end
if nargin <= 2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column2: text (%q)
%	column3: double (%f)
%   column4: double (%f)
%	column5: text (%q)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*q%q%f%f%q%[^\n\r]';

%% Open the text file.
fileID = fopen(fileName,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
geneInfo = table(dataArray{1:end-1}, 'VariableNames', {'acronym','entrez_id','gene_id','name'});