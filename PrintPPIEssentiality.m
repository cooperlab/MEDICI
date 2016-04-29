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
%Lines - N-length cell array of strings describing cell line name identifiers.
%Histology - N-length cell array of strings containing cell line histology identifiers.
%outputs:
%Table - M+2 x N+4 cell array of strings that can be written to disk
%           using PrintTable.m.

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

essen = cell(M,N);
for i=1:M
    
   for j=1:N
     essen{i,j} =  Essentialities{i}{j}{1};
   end
end

Table(3:end, 5:end) = essen;

end

