function Hash = GenerateHash(A)

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

Hash = java.util.Hashtable;

%load items from 'A' into 'Hash'
for i = 1:length(A)
    
    %check if collision
    item = Hash.get(A{i});
    
    %insert 'A{i}' into hash
    if(isempty(item))
        Hash.put(A{i}, i);
    else
        Hash.put(A{i}, [item.' i]);
    end
    
end
