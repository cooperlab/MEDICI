function [IDs iX iY] = IntersectSamples(IDx, X, IDy, Y)
%Takes intersection of two sets of samples, and returns associated data
%with corresponding samples in intersected set.
%inputs:
%IDx - N1-length cell array of strings containing sample IDs for sample set 1.
%X - N1 x M1 or M1 x N1 matrix or cell array of descriptions for set 1.
%IDy - N2-length cell array of strings containing sample IDs for sample set 2.
%Y - N2 x M2 or M2 x N2 matrix or cell array of descriptions for set 2.
%outputs:
%IDs - N-length cell array of strings containing sample IDs of intersection
%      set.
%iX - M x N or N x M matrix or cell array of descriptions for matched set 1.
%iY - M x N or N x M matrix or cell array of descriptions for matched set 2.

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

%match sample IDs
Mapping = StringMatch(IDx, IDy);

%check for multiple mappings
Lengths = cellfun(@(x)length(x), Mapping);
if(any(Lengths > 2))
    error('Multiple mapping violation.');
end

%find indices of corresponding mapped samples
Mapped1 = find(~cellfun(@isempty, Mapping));
Mapped2 = [Mapping{:}];

%output IDs
IDs = IDx(Mapped1);

%trim 'X'
if(length(IDx) == size(X,1))
    iX = X(Mapped1,:);
elseif(length(IDx) == size(X,2))
    iX = X(:,Mapped1);
else
    error('Size of X inconsistent with length(IDx).');
end

%trim 'Y'
if(length(IDy) == size(Y,1))
    iY = Y(Mapped2,:);
elseif(length(IDy) == size(Y,2))
    iY = Y(:,Mapped2);
else
    error('Size of Y inconsistent with length(IDy).');
end
