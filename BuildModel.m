function [Corrected, Connectivity, Novel, Screened, Source, Target, Undirected] = ...
            BuildModel(Adjacency)
%For use with 'GenerateSuperpathway.m'.
%Converts coded adjacency matrix into linear model format Y = AX, where Y
%is an N x 1 vector of protein essentiality scores from an RNAi screen, A
%is an N x P vector encoding the connectivity of the proteins (which proteins
%are connected to which edges), and X is a P x 1 vector of edge values.
%inputs:
%Adjacency - N x N adjacency matrix with the following codes:
%           1 - interaction from public, prior knowledge resource
%           2 - interaction from public resource, confirmed by PPI screen
%           3 - novel interaction from PPI screen
%           Directed interactions are asymetric (Adjacency(i,j) ~=
%           Adjacency(j,i)) and undirected interaction from PPI screens are
%           symmetric.
%outputs:
%Corrected - N x N adjacency matrix, updated to remove self and duplicate
%    interactions.
%Connectivity - N x P matrix representation of the model.
%Novel - P x 1 vector indicating which edges are novel (code 3).
%Screened - P x 1 vector indicating which edges are screened (code 2 or 3).
%Source - P x 1 vector showing source node index of edges (for debugging).
%Target - P x 1 vector showing target node index of edges (for debugging).
%Undirected - P x 1 vector indicating which edges are undirected (for debugging).

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

%get 'N'
N = size(Adjacency,1);

%enumerate edges - label each edge with a 1, 2, 3,...
[m n] = find(Adjacency > 0);
Novel = (Adjacency(sub2ind(size(Adjacency), m, n)) == 3);
Screened = (Adjacency(sub2ind(size(Adjacency), m, n)) == 2);

%initialize 'Undirected'
Undirected = false(size(m));

%generate set of duplicate edge pairs
Forward = cellfun(@(x,y)[num2str(x) '.' num2str(y)], num2cell(m), num2cell(n),...
                    'UniformOutput', false);
Reverse = cellfun(@(x,y)[num2str(x) '.' num2str(y)], num2cell(n), num2cell(m),...
                    'UniformOutput', false);
Mapping = StringMatch(Forward, Reverse);
Duplicates = [find(~cellfun(@isempty, Mapping)).' cell2mat(Mapping).'];
Duplicates(Duplicates(:,1) >= Duplicates(:,2),:) = [];

%set these edges to 'Undirected'
Undirected(Duplicates(:)) = true;

%delete entries below diagonal
%column 1 contains edges below diagonal, column 2 above diagonal
Index = sub2ind([size(Adjacency,1) size(Adjacency,2)],...
                            m(Duplicates(:,1)), n(Duplicates(:,1)));
Corrected = Adjacency;              
Corrected(sub2ind([size(Adjacency,1) size(Adjacency,2)],...
                            m(Duplicates(:,1)), n(Duplicates(:,1)))) = 0;
m(Duplicates(:,1)) = [];
n(Duplicates(:,1)) = [];
Novel(Duplicates(:,1)) = [];
Screened(Duplicates(:,1)) = [];
Undirected(Duplicates(:,1)) = [];

%remove self-interactions
SelfInteraction = find(m == n);
m(SelfInteraction) = [];
n(SelfInteraction) = [];
Novel(SelfInteraction) = [];
Screened(SelfInteraction) = [];
Undirected(SelfInteraction) = [];

%build edge sets for each node, up to layer (length(W))
%add immediate incoming edges last to simplify tracing nodes
EdgeSets = cell(N, 1);
for i = 1:N
    
    %initialize edge set
    EdgeSets{i} = cell(1);
    
    %record outgoing, undirected edges
    Outgoing = find(m == i); %outgoing edges
    Und = find(n == i & Undirected); %undirected edges
    EdgeSets{i}{1} = unique([Outgoing; Und]);
    
    %add incoming edges
    EdgeSets{i}{1} = union(EdgeSets{i}{1}, find(n == i));
    
end

%build connectivity matrix
Connectivity = zeros(N, length(m));
for i = 1:N
	Connectivity(i, EdgeSets{i}{1}) = 1;
end

%capture protein IDs
Source = m;
Target = n;
