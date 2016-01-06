function Mapping = StringMatch(A, B, Hash)

%create java hashtable
if(nargin == 2)
    Hash = GenerateHash(B);
end

%initialize output
Mapping = cell(1, length(A));

%Get mapping from 'A' items to 'B' items
for i = 1:length(A)
    
    %recover mapped indices
    Item = Hash.get(A{i});
    if(size(Item,1) > size(Item,2))
        Item = Item.';
    end
    Mapping{i} = Item;
    
end