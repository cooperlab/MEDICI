function DeLimited = ReLimit(Strings, Char)
%generate delimited string from cell arrray of strings.
%inputs:
%Strings - cell array of strings.
%Char - character separator.
%outputs:
%DeLimited - delimited string.

%check if input is empty
if(isempty(Strings))
    
    %return empty string
    DeLimited = '';
    
else   
    
    %append 
    Strings(1:end-1) = cellfun(@(x)[x Char ' '], Strings(1:end-1), 'UniformOutput', false);
    
    %concatenate
    DeLimited = [Strings{:}];
    
end

end