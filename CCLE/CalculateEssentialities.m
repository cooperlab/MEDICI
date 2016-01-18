function E = CalculateEssentialities(superpathwayFile,PathwaysTableFile,GeneEssentialitiesFile,TemplateReactionFile,alpha,w, Output)
%Builds context specific pathway for each cell line according to the
%mutaions and copy number values. 
%inputs:
%superpathwayFile - filename and path to superpathway file containing
%                   Sources and Targets of interactions.
%CNFile - filename and path to CN file which is the output file of the 
%         'CCLEBuildCNV.m' function. 
%CaptureFile - filename and path to the CaptureFile which is the output of
%              the 'CCLEHybridCapture.m'             
%Output - filename and path to store output variables.
%outputs
%A structure 'PathwaysTab

load(superpathwayFile)

load(PathwaysTableFile)
% load('SuperPathways_v2.mat') %% load pathway(7906*1016) SuperPathway_Lines(1016)
path_Lines = PathwaysTable.Lines;
path_Labels = PathwaysTable.Labels;
path_Values = PathwaysTable.Values;

load(TemplateReactionFile) % load template reactions source target

load('Histology')
histology = Histology;
hist_lines = Lines;

% load('CCLE.CopyNumber.mat')
% load('CCLE.CopyNumber.mat', 'CN') %% load variable CN (23316(genes)*1043(Lines))


% load('FullSuperpathway.strict.mat') %% load Proteins(1726) Source(7906) Target(7906) Adjacency(1726*2726) Connectivity(1726*7906)

load(GeneEssentialitiesFile)
% load('Achilles.mat') %% load Lines(216) ScoreSecond(12066*216)
ach_Lines = Lines;
scores = ScoresSecond;

%Remove cell lines in GeneEssentialitiesFile which are not in PathwaysTableFile
ach_unmapped_cellLines = find(~ismember(ach_Lines,path_Lines));
scores(:,ach_unmapped_cellLines) = [];
ach_Lines(ach_unmapped_cellLines) = [];

%Remove cell lines in PathwaysTableFile which are not in GeneEssentialitiesFile
path_unmapped_cellLines = find(~ismember(path_Lines,ach_Lines));
path_Lines(path_unmapped_cellLines) = [];
path_Values(:,path_unmapped_cellLines) = [];
path_Labels(:,path_unmapped_cellLines) = [];
%mutation.labels(path_unmapped_cellLines) = [];
%CopyNumber.labels(path_unmapped_cellLines) = [];


% hist_unmapped_cellLines = find(~ismember(hist_lines,path_Lines));
% hist_lines(hist_unmapped_cellLines) = [];
% histology(hist_unmapped_cellLines) = [];

% Map cell lines in both dataset
Mapping = StringMatch(path_Lines,ach_Lines);
mapping = cell2mat(Mapping)';
path_Values(:,mapping) = path_Values;
path_Labels(:,mapping) = path_Labels;
%path_Lines(mapping) = path_Lines;

%% Calculation part
[G,N] = size(path_Values);

% new_adjacency = Adjacency;
new_connectivity = zeros(size(Connectivity));
new_connectivity(Connectivity==1)=1;% direct neighbors only 

edge_essens = cell(N,1);
Essentialities = zeros(G,N);

for i = 1:N
    
    proteins = Proteins;
    p_sources = Source;
    p_targets = Target;
    Essentialities(:,i) = path_Values(:,i);
    
    connectivity = new_connectivity;
%     adjacency = new_adjacency;
    
    %Remove unavailable interactions in the pathway
    deleted_edge = find(isnan(path_Values(:,i)));
%     deleted_source = p_sources(deleted_edge);
%     deleted_target = p_targets(deleted_edge);
    p_sources(deleted_edge) = [];
    p_targets(deleted_edge) = [];
    connectivity(:,deleted_edge) = [];

    
    
    % Remove proteins which are not in Achilles
    ach_symbols = UniqueSymbols;
    Unmapped_ach = find(~ismember(ach_symbols,proteins));
    ach_symbols(Unmapped_ach) = [];
    c_scores = scores;
    c_scores(Unmapped_ach,:) = [];
    
    %Remove proteins with zero score from Achilles
    Unmapped_ach = find(c_scores(:,i)==0);
    ach_symbols(Unmapped_ach) = [];
    c_scores(Unmapped_ach,:) = [];
    
    
    %Find proteins with zero score in protien list to remove
    Unmapped_prot = find(~ismember(proteins,ach_symbols));
    %Find corresponding interactions
    unmapped_edges = find(cellfun(@(x)sum(strcmp(x, proteins(Unmapped_prot))),p_sources).' | cellfun(@(x)sum(strcmp(x,proteins(Unmapped_prot) )), p_targets).');
    p_sources(unmapped_edges) = [];
    p_targets(unmapped_edges) = [];
    proteins(Unmapped_prot) = [];
    connectivity(Unmapped_prot, :) = [];
%     adjacency(Unmapped_prot, :) = [];
%     adjacency(:,Unmapped_prot) = [];
    %Remove proteins and interactions from the network
    connectivity(:,unmapped_edges) = [];
    
    
    Mapping = StringMatch(ach_symbols,proteins);
    mapping = cell2mat(Mapping)';
    c_c_scores = zeros(size(c_scores));
%     c_ach_symbols= ach_symbols;
    c_c_scores(mapping,:) = c_scores;
%     c_ach_symbols(mapping) = ach_symbols;
    
    y = c_c_scores(:,i);
    y = - log10(y);
    
    %Remove template-reaction interactions
    Temp_interactions = strcat(source,target);
    p_interactions = strcat(p_sources,p_targets);
    mapped = StringMatch(Temp_interactions,p_interactions);
    unmapped_edges = cell2mat(mapped);
    
    p_sources(unmapped_edges) = [];
    p_targets(unmapped_edges) = [];
    connectivity(:,unmapped_edges) = [];
    
    %Calculate essentialities
    edge_essens{i} = get_edge_essen( connectivity,y,alpha,w);
  
    all_interactions = strcat(Source,Target);
    p_interactions = strcat(p_sources,p_targets);
    mapped = StringMatch(p_interactions,all_interactions);
    mapped_edges = cell2mat(mapped);
    
    Essentialities(mapped_edges,i) = edge_essens{i};   
end

E.Values = Essentialities;
E.Lines = ach_Lines;
E.Histology = histology;
E.Labels = path_Labels;
E.Source = Source;
E.Target = Target;

save (Output,E)

end
