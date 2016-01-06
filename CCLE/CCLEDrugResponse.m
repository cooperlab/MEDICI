function CCLEDrugResponse(DrugFile, LineFile, HUGOFile)
%generates CCLE treatment table from drug data and responses (EC50, IC05,
%AUC)
%inputs:
%DrugFile (string) - filename and path of the drug treatment data from CCLE
%                    (see Data folder for example: CCLE_NP24.2009_Drug_data_2012.02.20.txt).
%LineFile (string) - filename and path of the cell line annotation file
%                       from CCLE (see data folder for example: CCLE_sample_info_file_2012-10-18.txt).
%HUGOFile (string) - filename and path to a HUGO gene set file containing official gene symbols
%					 and past aliases. Example is available here: ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
%OutputFile (string) - filename and path to store outputs.
%outputs:
%CCLE (structure) - a Matlab structure with the following fields:
%						CCLE.Lines - N-length cell array of strings describing cell line names
%						CCLE.Compounds - M-length cell array of strings describing compound names
%						CCLE.Targets - M-length cell array of strings describing compound target proteins
%						CCLE.EC50 - MxN matrix of EC50 values
%						CCLE.IC50 - MxN matrix of IC50 values
%						CCLE.Amax - MxN matrix of activation function maximum values
%						CCLE.AUC - MxN matrix of activation function area-under-curve values

%load CCLE descriptions
Descriptions = text2cell(LineFile, '\t');

%get cell line names, determine hybrid sequence status
Lines = Descriptions(2:end, find(strcmp(Descriptions(1,:), 'CCLE name')));

%read in drug response file
Responses = text2cell(DrugFile, '\t');
TreatedLines = Responses(2:end, find(strcmp(Responses(1,:), 'CCLE Cell Line Name')));
TreatedCompounds = Responses(2:end, find(strcmp(Responses(1,:), 'Compound')));
TreatedTargets = Responses(2:end, find(strcmp(Responses(1,:), 'Target')));
TreatedEC50 = cellfun(@(x)str2double(x), Responses(2:end, find(strcmp(Responses(1,:), 'EC50 (uM)'))));
TreatedIC50 = cellfun(@(x)str2double(x), Responses(2:end, find(strcmp(Responses(1,:), 'IC50 (uM)'))));
TreatedAmax = cellfun(@(x)str2double(x), Responses(2:end, find(strcmp(Responses(1,:), 'Amax'))));
TreatedArea = cellfun(@(x)str2double(x), Responses(2:end, find(strcmp(Responses(1,:), 'ActArea'))));

%get list of unique compounds
[Compounds, Indices] = unique(TreatedCompounds);
Targets = TreatedTargets(Indices);

%fill in tables
EC50 = nan(length(Compounds), length(Lines));
IC50 = nan(length(Compounds), length(Lines));
Amax = nan(length(Compounds), length(Lines));
ActArea = nan(length(Compounds), length(Lines));
for i = 1:length(Lines)
    
    %identify treatment entries in current line
    Hits = find(strcmp(TreatedLines, Lines{i}));
    
    %identify treatments
    HitCompounds = TreatedCompounds(Hits);
    Mapping = StringMatch(HitCompounds, Compounds);
    EC50([Mapping{:}], i) = TreatedEC50(Hits);
    IC50([Mapping{:}], i) = TreatedIC50(Hits);
    Amax([Mapping{:}], i) = TreatedAmax(Hits);
    ActArea([Mapping{:}], i) = TreatedArea(Hits);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correct drug target list

%table to correct problems with compound target descriptions
CorrectionTable = {'AZD6244', 'MEK', {'MAPK1', 'MAPK2'};
                    'L-685458', 'GS', 'NOTCH1';
                    'PD-0325901', 'MEK', {'MAPK1', 'MAPK2'};
                    'PF2341066', 'c-MET', 'MET';
                    'PHA-665752', 'c-MET', 'MET';
                    'PLX4720', 'RAF', 'BRAF';
                    'RAF265', 'RAF', 'BRAF';
                    'Sorafenib', 'RTK', {'RAF1', 'BRAF', 'KIT', 'FLT3', 'KDR', 'FLT4', 'PDGFRB'};
                    'TKI258', 'FGFR', {'FGFR1', 'FGFR2', 'FGFR3'}};

%generate HUGO structure
HUGO = ParseHUGO(HUGOFile);

%map compound targets to HUGO symbols
TargetsHUGO = HUGOLookup(Targets, HUGO, 'Symbol');

%fill in missing info from correction table
Mapping = StringMatch(Compounds, CorrectionTable(:,1));
TargetsHUGO(~cellfun(@isempty, Mapping)) = CorrectionTable([Mapping{:}], 3);

%build output structure
CCLE.Lines = Lines;
CCLE.Compounds = Compounds;
CCLE.Targets = TargetsHUGO;
CCLE.EC50 = EC50;
CCLE.IC50 = IC50;
CCLE.Amax = Amax;
CCLE.AUC = ActArea;

%save tables
save(OutputFile, 'CCLE');