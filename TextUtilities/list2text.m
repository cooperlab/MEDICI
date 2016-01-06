function success = list2text(A, filename)
%Write contents of N-length 1D cell array to text file.
%inputs:
%A - 1D cell array of strings
%filename - string indicating path and filename for output

fid = fopen(filename, 'w');

if(fid > 0)
    for i = 1:length(A)
        fprintf(fid, '%s\n', A{i});
    end
    
    if(fclose(fid) > 0);
        success = 1;
    else
        success = 0;
    end
else
    success = 0;
end