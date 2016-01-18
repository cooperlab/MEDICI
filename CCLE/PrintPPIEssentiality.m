function Table = PrintPPIEssentiality(Essentialities, Source, Target,...
                                        Novel, Pathways, Lines, Histology)
%Generates a spreadsheet containing PPI essentialities
%inputs:
%Essentialities - M x N matrix containing PPI essentiality values. Cell
%                 lines in columns, PPIs in rows.
%Source - M-length cell array of strings, listing source of PPI
%         interaction.
%Target - M-length cell array of strings, listing target of PPI
%         interaction.
%Novel - M-length logical vector (0,1) indicating which interactions are
%        novel.
%Pathways - M-length cell array of strings containing comma-separated 
%           sub-pathways where each interaction is found.
%Lines - N-length cell array 
%outputs:
%Contents - M+2 x N+4 cell array of strings that can be written to disk
%           using cell2text.m.

%get dimensions of input 'Essentialities'
M = length(Source); %number of PPIs
N = length(Lines); %number of cell lines

%allocate output table
Table = cell(M+2, N+4);

%capture source, target, novelty and pathway strings
Table(2,1:4) = {'Source', 'Target', 'Novel', 'Pathways'};
Table(3:end, 1) = Source;
Table(3:end, 2) = Target;
Table(3:end, 3) = cellfun(@(x)sprintf('%g', x), num2cell(double(Novel)),...
                            'UniformOutput', false);
Table(3:end, 4) = Pathways;

%capture cell line names and histology
Table(1, 5:end) = Lines;
Table(2, 5:end) = Histology;

%convert essentialities to strings
% Table(3:end, 5:end) = cellfun(@(x)sprintf('%g', x),...
%                          num2cell(Essentialities), 'UniformOutput', false);

essen = cell(M,N);
for i=1:M
    
   for j=1:N
     essen{i,j} =  Essentialities{i}{j}{1};
   end
end

Table(3:end, 5:end) = essen;


end

