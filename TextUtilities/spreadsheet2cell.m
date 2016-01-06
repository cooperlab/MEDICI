function [row column data] = spreadsheet2cell(filename, delimiter)
%same as text2cell, except assumes one header row and one header column,
%rest are numbers.  Avoids expensive 'num2str' conversions of massive cell
%arrays.

fid = fopen(filename, 'r');

if(fid >= 3) %file opened successfully
    %get first line, determine number of columns
    row = textscan(fgetl(fid), '%s', 'Delimiter', delimiter); row = {row{1}{2:end}}; 
    columns = length(row);
      
    %create format string
    format = ['%s '];
    for i = 1:columns
        format = [format '%f '];
    end
    
    %scan remainder of file
    text = textscan(fid, format, 'delimiter', delimiter);
        
    %extract header column
    column = {text{1}{:}};
    
    %extract data
    data = [text{2:end}];

    %close file
    fclose(fid);
    
else %file error
    A = {};
end