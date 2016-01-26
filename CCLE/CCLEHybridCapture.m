function CCLEHybridCapture(MAFFile, LineFile, OutputFile)
%Generates CCLE mutation table from hybrid capture sequencing data.
%inputs:
%MAFFile (string) - filename and path of the MAF file containing hybrid capture data
%					(see Data folder for example: CCLE_hybrid_capture1650_hg19_damaging.maf).
%LineFile (string) - filename and path of the cell line annotation file
%                       from CCLE (see Data folder for example: CCLE_sample_info_file_2012-10-18.txt).
%OutputFile (string) - filename and path to store outputs.
%Capture (structure) - a Matlab structure containing the following fields:
%						Capture.Lines - N-length cell array of strings describing cell line names
%						Capture.Lineage - N-length cell array of strings describing cell line lineage
%						Capture.Mutations - MxN binary matrix indicating mutations (1 --> mutated)
%						Capture.Symbols - M-length cell array of strings describing gene symbols 
%											from hybrid capture platform.
%                       Capture.Labels - MxN cell array contaning the type
%                       of mutation for each gene in each cell line

%load CCLE descriptions
Descriptions = text2cell(LineFile, '\t');

%get cell line names, determine hybrid sequence status
Lines = Descriptions(2:end, find(strcmp(Descriptions(1,:), 'CCLE name')));
Captures = Descriptions(2:end, find(strcmp(Descriptions(1,:), 'Hybrid Capture Sequencing')));
Lineage = Descriptions(2:end, find(strcmp(Descriptions(1,:), 'Site Primary')));
Histology = Descriptions(2:end, find(strcmp(Descriptions(1,:), 'Histology')));
HistSubtype = Descriptions(2:end, find(strcmp(Descriptions(1,:), 'Hist Subtype1')));
Hits = find(~cellfun(@isempty, HistSubtype));
Lineage(Hits) = cellfun(@(x,y)[x '.' y], Lineage(Hits), HistSubtype(Hits),...
                    'UniformOutput', false);

%read in .maf file
CaptureMAF = text2cell(MAFFile, '\t');
MutatedLines = CaptureMAF(2:end, find(strcmp(CaptureMAF(1,:), 'Tumor_Sample_Barcode')));
Mutations = CaptureMAF(2:end, find(strcmp(CaptureMAF(1,:), 'Hugo_Symbol')));

%get list of unique gene symbols
Symbols = unique(CaptureMAF(2:end, find(strcmp(CaptureMAF(1,:), 'Hugo_Symbol'))));

%remove silent mutations
VariantClass = CaptureMAF(2:end, find(strcmp(CaptureMAF(1,:), 'Variant_Classification')));
Silent = strcmp(VariantClass, 'Silent');
CaptureMAF([false; Silent],:) = [];

%fill in tables
CaptureMutations = zeros(length(Symbols), length(Lines));
Mutation_Labels = cell(1,length(Lines));
no_capture = {'no_capture'};
nocapture_marker = repmat(no_capture,length(Symbols),1);

for i = 1:length(Lines)
    if(strcmp(Captures{i}, 'yes'))
        muts_type = VariantClass;
        Hits = find(strcmp(MutatedLines, Lines{i}));
        HitSymbols = Mutations(Hits);
        muts_type = muts_type(Hits);
        muts_value = zeros(size(muts_type));
        
        z = {''};
        muts_type_label = repmat(z,size(muts_type));
        muts_marker = repmat(z,length(Symbols),1);
        
        for j= 1:length(muts_value)
            if(strcmp(muts_type{j},'Frame_Shift_Del') || strcmp(muts_type{j}, 'Frame_Shift_Ins') || strcmp(muts_type{j}, 'Nonsense_Mutation'))
                muts_value(j) = 1;
                muts_type_label{j} = muts_type{j};
            end
        end
        
        Mapping = StringMatch(HitSymbols, Symbols);
        indx = cell2mat(Mapping)';
        CaptureMutations(indx,i) = muts_value;
        
        muts_marker(indx) = muts_type_label;
        Mutation_Labels{i} = muts_marker;
    else
        CaptureMutations(:, i) = NaN;
        Mutation_Labels{i} = nocapture_marker;
    end
end

%build output structure
Capture.Lines = Lines;
Capture.Lineage = Lineage;
Capture.Mutations = CaptureMutations;
Capture.Symbols = Symbols;
Capture.Labels = Mutation_Labels;
Capture.Histology = Histology;

%save tables
save(OutputFile, 'Capture');
