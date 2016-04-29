function GenerateCytoscape(OutputPath, EdgeFolder, OncogeneTSPFile)
%Generate cytoscape visualization files to view node and interaction essentiality. This function
%generates 1 node table and 1 edge table for each cell line. Each pair of files can be visualized
%in Cytoscape using the visualization style 'EdgeStyle.xml' located in the Visualization folder.
%inputs:
%OutputPath (string) - path to folder for storing node and edge table output files for Cytoscape visualization.
%EdgeFolder (string) - path to folder where edge essentiality .mat files are stored (one file for each cell line).
%OncogeneTSPFile - filename and path of a text file containing the names of oncogenes and 
%					tumor suppressors. This file contains two columns, the first contains the gene
%					symbol, and the second contains the designation 'TSG' or 'Oncogene'.
%					(see Data folder for example file 'VogelsteinGenes.txt'.

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

%parameters
EdgeOutputSuffix = '.EdgeCS.txt';
NodeOutputSuffix = '.NodeCS.txt';

%get .mat file list in input folder
MatFiles = dir([EdgeFolder '*.mat']);
MatFiles = {MatFiles.name};
Dot = strfind(MatFiles, '.');
MatLines = cellfun(@(x,y)x(1:y(2)-1), MatFiles, Dot, 'UniformOutput', false);

%load first file, build matching string
load([EdgeFolder MatFiles{1}]);
Matching = cellfun(@(x,y)[x y], Source, Target, 'UniformOutput', false);

%create lead information that will be common to all tables
EdgeHeader = [{'Source', 'Target', 'Novel', 'Screened',...
                'Undirected'};...
                Source Target...
                cellfun(@(x)num2str(x), num2cell(Novel), 'UniformOutput', false)...
                cellfun(@(x)num2str(x), num2cell(Screened), 'UniformOutput', false)...
                cellfun(@(x)num2str(x), num2cell(Undirected), 'UniformOutput', false)];

%build null model
Pooled = cell(1,length(MatLines));
for i = 1:length(MatLines)
    
    %load .mat file
    load([EdgeFolder MatFiles{i}]);
    
    %capture null values
    Pooled{i} = Nulls;
    
end
Pooled = [Pooled{:}];
Pooled = Pooled(:);
NullBins = [min(Pooled):(max(Pooled)-min(Pooled))/1000:max(Pooled)];
Null = histc(Pooled, NullBins);
Null = Null / sum(Null);

%load vogelstein genes
OncogenesTSPs = text2cell(OncogeneTSPFile, '\t');
OncogeneTSPgenes = OncogenesTSPs(:,1);
OncogeneTSPstatus = OncogenesTSPs(:,2);
            
%process each cell line, bulding edge and node tables
for i  = 1:length(MatLines)
    
    %load .mat file
    load([EdgeFolder MatFiles{i}]);
    
    %calculate edge ranks, significance
    Ranks = tiedrank(EdgeEssentiality);
    Significance = ones(length(EdgeEssentiality),1);
    for j = 1:length(EdgeEssentiality)
        Diff = EdgeEssentiality(j) - NullBins;
        Bin = find(Diff <= 0, 1, 'first');
        if(~isempty(Bin))
            Significance(j) = sum(Null(1:Bin));
        end        
    end
    
    %capture table elements
    EdgeEssTable = [{'EdgeEssentiality'};...
                        cellfun(@(x)num2str(x), num2cell(EdgeEssentiality), 'UniformOutput', false)];
    EdgeAbsEssTable = [{'abs(EdgeEssentiality)'};...
                        cellfun(@(x)num2str(x), num2cell(abs(EdgeEssentiality)), 'UniformOutput', false)];    
    EdgeRankTable = [{'EdgeRank'};...
                            cellfun(@(x)num2str(x), num2cell(Ranks), 'UniformOutput', false)];
    SigTable = [{'p-value'};...
                            cellfun(@(x)num2str(x), num2cell(Significance), 'UniformOutput', false)];
    LogSigTable = [{'log10(p-value)'};...
                            cellfun(@(x)num2str(x), num2cell(log(Significance)), 'UniformOutput', false)];
                        
    %build edge table, write to disk
    cell2text([EdgeHeader EdgeEssTable EdgeAbsEssTable EdgeRankTable SigTable],...
            [OutputPath MatLines{i} EdgeOutputSuffix]);
        
    %generate unique list of nodes
    Nodes = union(Source, Target);
    
    %get essentiality scores
    NodeEssentiality = zeros(length(Nodes), 1);
    Mapping = StringMatch(Nodes, Source);
    for j = 1:length(Mapping)
        if(~isempty(Mapping{j}))
            Mapping{j} = Mapping{j}(1);
        end
    end
    NodeEssentiality(~cellfun(@isempty, Mapping)) = SourceProteinEssentiality([Mapping{:}]);
    Mapping = StringMatch(Nodes, Target);
    for j = 1:length(Mapping)
        if(~isempty(Mapping{j}))
            Mapping{j} = Mapping{j}(1);
        end
    end
    NodeEssentiality(~cellfun(@isempty, Mapping)) = TargetProteinEssentiality([Mapping{:}]);

    %map to vogelstein genes
    Vogelstein = cell(length(Nodes), 1);
    Mapping = StringMatch(Nodes, OncogeneTSPgenes);
    Vogelstein(~cellfun(@isempty, Mapping)) = OncogeneTSPstatus([Mapping{:}]);
    Vogelstein(cellfun(@isempty, Vogelstein)) = {'NA'};
    
    %capture table elements
    NodeAbsEssentiality = cellfun(@(x)num2str(x), num2cell(abs(NodeEssentiality)), 'UniformOutput', false);
    NodeEssentiality = cellfun(@(x)num2str(x), num2cell(NodeEssentiality), 'UniformOutput', false);
    cell2text([{'Symbol', 'ProteinEssentiality', 'abs(ProteinEssentiality)', 'OncogeneTSP'}; Nodes NodeEssentiality NodeAbsEssentiality Vogelstein],...
                [OutputPath MatLines{i} NodeOutputSuffix]);
    
end
