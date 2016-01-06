function CNV = CCLEBuildCNV(CNVFile, Output)
%Extracts integral copy number values from GISTIC analysis and labels genes
%in significant focal and broad events identified by GISTIC.
%inputs:
%CNVFile - filename and path to gene-centric CNV file for cell lines. 
%           (see Data folder for example file: CCLE_copynumber_byGene_2013-12-03.txt)
%Output - filename and path to store output variables.
%outputs:
%A structure 'CN' with the following fields:
%   CN.Lines - N-length cell array of strings describing cell line identifiers.
%   CN.Symbols - an M-length cell array of strings describing gene symbols.
%   CN.Chromosome - an M-length cell array of strings describing the chromosome
%               where the gene is located.
%   CN.Start - an M-length vector containing gene start base location.
%   CN.Stop - an M-length vector containing gene stop base location.
%   CN.Values - an M x N array of copy number values.

%import copy number data
Contents = text2cell(CNVFile, '\t');
Symbols = Contents(2:end,2);
Chromosome = Contents(2:end,3);
Start = cellfun(@(x)str2double(x), Contents(2:end,4));
Stop = cellfun(@(x)str2double(x), Contents(2:end,5));
CN = cellfun(@(x)str2double(x), Contents(2:end, 6:end));
Lines = Contents(1,6:end);

%build output structure
CN.Symbols = Symbols;
CN.Chromosome = Chromosome;
CN.Start = Start;
CN.Stop = Stop;
CNV.CNV = Values;
CN.Lines = Lines;

%save data
save(Output, 'CN');