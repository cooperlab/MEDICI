function GeneCentricEssentiality(GCTFile, CellDescriptionFile, HUGOFile, B, ESNullFile, OutputFile)
%Generate gene-centric essentiality values from shRNA or CRISPR screens with HUGO nomenclature.
%Inputs:
%GCTFile (string) - filename and path to a GCT file containing essentiality values from 
%					shRNA or CRISPR screens. The first column contains the probe identifier.
%					The second column contains the corresponding symbol of the target gene.
%					The third column and beyond contains essentiality values representing the
%					screen readout of log(after/before). Here we assume 0 --> neutral, negative
%					--> lethal, and positive --> proliferation.
%					(see Data folder for example file: Achilles.QCv2.4.3.rnai.gct).
%					The first row (columns 3 and beyond) contain cell line identifiers.
%CellDescriptionFile (string) - filename and path to a tab-delimited text file containing cell-line
%								information in columns. The cancer type column is indicated by the 
%								header 'Type' in the first row, the subtype column is indicated by 
%								'Subtype', and the cell line name is indicated by 'Name'.
%								(see Data folder for example file: Achilles.CellLineData.txt)
%HUGOFile (string) - filename and path to a HUGO gene set file containing official gene symbols
%					 and past aliases. Example is available here: ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
%B (scalar) - number of resamplings used in estimating shRNA significance
%ESNullFile (string) - filename and path to a file containing null-distribution values for 
%			  Kolmogorov-Smirnov enrichment scores. If the file exists from a previous 
%			  analysis, providing the path and filename will reduce computation time.
%			  If it does not already exist, the null enrichment scores will be saved at the
%			  provided location.
%OutputFile (string) - filename and path to store function outputs in .mat format. Variables
%					   describing the probes, gene symbols, unique gen symbols, correspondence
%						between probes and symbols will be stored, along with gene-centric 
%						essentiality scores aggregated by three methods: most lethal probe per gene, 
%						second most lethal probe per gene, and KS enrichment of all probes per gene.

%build HUGO structure
HUGO = ParseHUGO(HUGOFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read in shRNA values from .gct - Achilles or the like

Contents = text2cell(GCTFile, '\t');
Probes = Contents(2:end,1);
Symbols = Contents(2:end,2);
Lines = Contents(1,3:end);
Contents = Contents(2:end,3:end);
Symbols = HUGOLookup(Symbols, HUGO, 'Symbol');
Missing = cellfun(@isempty, Symbols);
Symbols(Missing) = [];
Contents(Missing,:) = [];
[UniqueSymbols, ~, Correspondence] = unique(Symbols);
Scores = cellfun(@(x)str2double(x), Contents);
Scores = knnimpute(Scores, 2);
clear Contents;

%read in site and histology information, map to cell lines
Contents = text2cell(CellDescriptionFile, '\t');
TypeCol = find(strcmp(Contents(1,:), 'Type'));
SubtypeCol = find(strcmp(Contents(1,:), 'Subtype'));
NameCol = find(strcmp(Contents(1,:), 'Name'));
Strings = cellfun(@(x,y) [x '.' y], Contents(2:end, TypeCol),...
                        Contents(2:end, SubtypeCol), 'UniformOutput', false);
Mapping = StringMatch(Lines, Contents(2:end, NameCol));
Histology = cell(size(Lines));
Histology(~cellfun(@isempty, Mapping)) = Strings([Mapping{:}]);

%calculate scores using top ranking shRNA (most lethal) %%%%%%%%%%%%%%%%%%%

%sort for each gene and each cell line
ScoresMax = zeros(length(UniqueSymbols), size(Scores,2));
for i = 1:length(UniqueSymbols)
    List = Correspondence == i;
    ScoresMax(i,:) = min(Scores(List,:), [], 1);
end


%calculate scores using second best rank (second most lethal) %%%%%%%%%%%%%

%generate null distribution of second-rank values
MaxProbes = max(hist(Correspondence, [1:length(UniqueSymbols)])); %find max number of probes
Null = zeros(B, MaxProbes);
for i = 1:MaxProbes
    Indices = ceil(size(Scores,1)*size(Scores,2) * rand(B*i,1));
    Values = reshape(Scores(Indices), [B i]);
    if i > 1
        Values = sort(Values, 2);
        Null(:,i) = Values(:,2);
    else
        Null(:,i) = Values;
    end
end

%scan through each gene and each line, calculate p-value for second most
%lethal shRNA
ScoresSecond = zeros(length(UniqueSymbols), length(Lines));
for i = 1:length(UniqueSymbols)
    List = Correspondence == i; %get indices of shRNAs for gene 'i'
    N = sum(List); %get number of shRNAs for gene 'i'
    Sorted = sort(Scores(List,:), 1);
    Second = Sorted(min(2, N), :);
    for j = 1:length(Lines)
        ScoresSecond(i,j) = sum(Null(:,N) <= Second(j)) / B;
    end
end


%calculate scores using KS statistic of all shRNAs %%%%%%%%%%%%%%%%%%%%%%%%

%generate null distributions of KS statistics for gene sets with varying sizes
N = numel(Scores); %total number of shRNA measurements across all lines
[~, Sorted] = sort(Scores(:)); %rank all shRNA scores
Ranks(Sorted) = 1:N;
Ranks = reshape(Ranks, size(Scores));
if(~exist(ESNullFile, 'file')) %generate null ES scores for gene sets of size 1:MaxProbes
    NullES = zeros(B, MaxProbes);
    for i = 1:MaxProbes
        NG = i; %number of shRNAs in set
        Indices = ceil(size(Scores,1)*size(Scores,2) * rand(B*i,1)); %random sampling of B sets of NG shRNAs
        NullData = reshape(Ranks(Indices), [B i]);
        
        for j = 1:B %calculate ES for each sampled shRNA set
            
            [i j]
            
            CDF = zeros(1,N); %rank CDF for 'in' shRNAs
            CDF(NullData(j,:)) = 1;
            Pg = cumsum(CDF) / NG;
            
            CDF = ones(1,N); %rank CDF for control 'out' shRNAs
            CDF(NullData(j,:)) = 0;
            Pn = cumsum(CDF) / (N - NG);
            
            Difference = Pg - Pn; %calculate enrichment score
            [~, Index] = max(abs(Difference));
            NullES(j,i) = Difference(Index);
            
        end
    end
    save(ESNullFile, 'NullES', 'MaxProbes');
else
    load(ESNullFile);
end

%calculate actual enrichment scores using KS of shRNAs corresponding to each gene
N = size(Scores,1);
[~, Sorted] = sort(Scores,1); %rank all shRNA scores
Ranks = zeros(N, length(Lines));
for i = 1:length(Lines)
    Ranks(Sorted(:,i),i) = 1:N;
end
ScoresES = zeros(length(UniqueSymbols), length(Lines));
for i = 1:length(UniqueSymbols)
    
    i
    
	List = Correspondence == i; %get indices of shRNAs for gene 'i'
    NG = sum(List); %get number of shRNAs for gene 'i'
    
    for j = 1:length(Lines)
        
        CDF = zeros(1,N); %rank CDF for 'in' shRNAs
        CDF(Ranks(List,j)) = 1;
        Pg = cumsum(CDF) / NG;
        
        CDF = ones(1,N); %rank CDF for control 'out' shRNAs
        CDF(Ranks(List,j)) = 0;
        Pn = cumsum(CDF) / (N - NG);
        
        Difference = Pg - Pn; %calculate enrichment score
        [~, Index] = max(abs(Difference));
        ES = Difference(Index);
        
        ScoresES(i,j) = sum(NullES(:,NG) <= ES) / B; %calculate p-value
        
    end
end

%write outputs
save(OutputFile, 'Probes', 'Symbols', 'UniqueSymbols', 'Correspondence',...
        'ScoresMax', 'ScoresSecond', 'ScoresES', 'Lines');
