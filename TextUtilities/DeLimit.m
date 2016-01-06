function DeLimited = DeLimit(String, Char)
%separate comma delimited list into cell arrray of strings.
%inputs:
%String - string, list of symbols.
%Char - character separator.
%outputs:
%DeLimited - cell array of strings.

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