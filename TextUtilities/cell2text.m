function result = cell2text(A, filename)

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

try
    fid = java.io.FileOutputStream(filename);
catch exception
    disp(exception.message);
    result = 0;
    return;
end
try
    streamBuffer = java.io.BufferedOutputStream(fid, 65536);
catch exception
    disp(exception.message);
    result = 0;
    return;
end

%determine number of rows, columns in output
columns = size(A,2);
rows = size(A,1);

for i = 1:size(A,1)
    
%     for j = 1:size(A,2)-1
%         streamBuffer.write(uint8([A{i,j} sprintf('\t')]), 0, length(A{i,j})+1);
%     end
%     streamBuffer.write(uint8([A{i,end} sprintf('\n')]), 0, length(A{i,j})+1);
    
    row = A(i,:);
    row(1:end-1) = cellfun(@(x) [x sprintf('\t')], row(1:end-1), 'UniformOutput', false);
    row{end} = [row{end} sprintf('\n')];
    string = [row{:}];
    streamBuffer.write(uint8(string), 0, length(string));
    
    %[i toc]
end

%flush buffer
streamBuffer.flush();

%close file
fid.close();

result = 1;
