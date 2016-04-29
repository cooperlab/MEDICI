function DeLimited = DeLimit(String, Char)
%separate comma delimited list into cell arrray of strings.
%inputs:
%String - string, list of symbols.
%Char - character separator.
%outputs:
%DeLimited - cell array of strings.

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

%check if input is empty
if(isempty(String))
    
    %return empty cell array
    DeLimited = {};
    
else   
    
    %get number of elements
    Elements = length(strfind(String, Char)) + 1;
    
    %initialize output
    DeLimited = cell(1, Elements);
    
    %if more than one symbols
    if(Elements > 1)
        
        %strip out each element
        for i = 1:Elements
            
            %get leading token
            [Single, String] = strtok(String, Char);
            
            %take off comma and record
            DeLimited{i} = strtrim(Single);
            
        end
        
    else %only one symbol, return
        
        DeLimited{1} = strtrim(String);
        
    end
    
end

end
