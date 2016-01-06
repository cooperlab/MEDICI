function Hash = GenerateHash(A)

%
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