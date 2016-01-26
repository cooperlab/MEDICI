function CN = CCLEBuildCNV(CNVFile, Output)
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
%   CN.Labels - an M x N cell array of strings describing which genes have
%   copynumber less than -1 

%import copy number data
Contents = text2cell(CNVFile, '\t');
Symbols = Contents(2:end,2);
Chromosome = Contents(2:end,3);
Start = cellfun(@(x)str2double(x), Contents(2:end,4));
Stop = cellfun(@(x)str2double(x), Contents(2:end,5));
Values = cellfun(@(x)str2double(x), Contents(2:end, 6:end));
Lines = Contents(1,6:end);


Labels = cell(1,length(Lines));
z = {''};
marker = repmat(z,size(Values,1),1);

%Find genes with copy number less than -1 
for i = 1: length(Lines)
    ind = find(Values(:,i)<-1);
    c = {'CN'};
    temp = marker;
    temp(ind) = c;
    Labels{i} = temp;
end

%build output structure
CN.Symbols = Symbols;
CN.Chromosome = Chromosome;
CN.Start = Start;
CN.Stop = Stop;
CN.CNV = Values;
CN.Lines = Lines;
CN.Labels = Labels;


%save data
save(Output, 'CN','-v7.3');