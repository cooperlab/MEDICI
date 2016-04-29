function Essentialities = CalculateEssentialities(SuperpathwayFile, PathwaysTableFile, GeneEssentialitiesFile, TemplateReactionFile, alpha, w, Output)
%Calculates PPIs essentialities in each cell line 
%inputs:
%SuperpathwayFile - filename and path to superpathway file containing
%                   Sources and Targets of interactions.
%PathwaysTableFile - filename and path to the file containing context 
%                    specific pathways which is the output of the
%                    CCLEBuildPathwaysTable.m
%GeneEssentialitiesFile - filename and path to the file contating 
%                         gene-centric essentiality values.
%TemplateReactionFile - filename and path to the file containing
%                       interactions that appear only with the TargetType 
%                       as TemplateReaction.
%alpha - scalar, smoothing parameter used to weight nieghboring essentialities for 
%		each network node when estimating interaction essentialities. A value of '0' means
%		that each node should be handled independently of its neighbors.
%w - scalar, weighting parameter for network inversion. w represents how 
%       much the essentiality of proteins affects the essentiality of interactions
%       in between. Larger values of w show that an interaction is essential 
%       if the corresponding proteins are essential, and smaller values of w show 
%       that no matter how essential the corresponding proteins are the
%       interaction is important.
%Output - filename and path to store a structure containing output variables.
%						This structure contains the following:
%						Values - PPI esseantialitiy scores
%						Lines - Cell lines name
%						Histology - histology
%						Labels - Type of genetic alterations
%						Source - Names of proteins in one side of
%						ineractions
%						Target - Names of proteins in another side of
%						interactions

%Licensed to the Apache Software Foundation (ASF) under one
%or more contributor license agreements.  See the NOTICE file
%distributed with this work for additional information
%regarding copyright ownership.  The ASF licenses this file
%to you under the Apache License, Version 2.0 (the
%"License"); you may not use this file except in compliance
%with the License.  You may obtain a copy of the License at
%
%  http://www.apache.org/licenses/LICENSE-2.0
%
%Unless required by applicable law or agreed to in writing,
%software distributed under the License is distributed on an
%"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
%KIND, either express or implied.  See the License for the
%specific language governing permissions and limitations
%under the License.

load(SuperpathwayFile);

load(PathwaysTableFile)
% load('SuperPathways_v2.mat') %% load pathway(7906*1016) SuperPathway_Lines(1016)
path_Lines = PathwaysTable.Lines;
path_Labels = PathwaysTable.Labels;
path_Values = PathwaysTable.Values;
histology = PathwaysTable.Histology;

load(TemplateReactionFile) % load template reactions source target


load(GeneEssentialitiesFile)
ach_Lines = Lines;
scores = ScoresSecond;

%Remove cell lines in GeneEssentialitiesFile which are not in PathwaysTableFile
ach_unmapped_cellLines = find(~ismember(ach_Lines,path_Lines));
scores(:,ach_unmapped_cellLines) = [];
ach_Lines(ach_unmapped_cellLines) = [];

%Remove cell lines in PathwaysTableFile which are not in GeneEssentialitiesFile
path_unmapped_cellLines = find(~ismember(path_Lines,ach_Lines));
path_Lines(path_unmapped_cellLines) = [];
histology(path_unmapped_cellLines) = [];
path_Values(:,path_unmapped_cellLines) = [];
path_Labels(:,2:end) = [];
for k = 1:size(path_Labels,1)
    path_Labels{k}(path_unmapped_cellLines) = [];
end
% path_Labels(:,path_unmapped_cellLines) = [];

% Map cell lines in both dataset
Mapping = StringMatch(path_Lines,ach_Lines);
mapping = cell2mat(Mapping)';
path_Values(:,mapping) = path_Values;
for k = 1:size(path_Labels,1)
    path_Labels{k} = path_Labels{k}(mapping);
end
% path_Labels(:,mapping) = path_Labels;

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
        
    %Remove unavailable interactions in the pathway
    deleted_edge = find(isnan(path_Values(:,i)));
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
    
    %Remove proteins and interactions from the network
    connectivity(:,unmapped_edges) = [];
    
    Mapping = StringMatch(ach_symbols,proteins);
    mapping = cell2mat(Mapping)';
    c_c_scores = zeros(size(c_scores));
    c_c_scores(mapping,:) = c_scores;    
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
    edge_essens{i} = GetEdgeEssentiality( connectivity,y,alpha,w);
  
    all_interactions = strcat(Source,Target);
    p_interactions = strcat(p_sources,p_targets);
    mapped = StringMatch(p_interactions,all_interactions);
    mapped_edges = cell2mat(mapped);
    
    Essentialities(mapped_edges,i) = edge_essens{i};   
end

Essens = abs(Essentialities);
clear Essentialities; 

for i=1:size(Essens,2)
    temp = Essens(:,i);
    temp(find(isnan(temp)))=[];
    %temp(find(temp==0))=[];
    temp = abs(temp);
    Essens(:,i) = abs(Essens(:,i));
    Essens(:,i) = (Essens(:,i)-min(temp))/(max(temp)-min(temp));
end

Essentialities.Values = Essens;
Essentialities.Lines = ach_Lines;
Essentialities.Histology = histology;
Essentialities.Labels = path_Labels;
Essentialities.Source = Source;
Essentialities.Target = Target;

save(Output,'Essentialities')
end
