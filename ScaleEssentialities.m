function ScaleEssentialities(EssentialitiesFile, Output)
%Normalizes interaction essentiality values by scaling to preserve significance of positive
%and negative signs on values from gene essentiality experiments. A negative gene essentiality
%indicates a gene that inhibits proliferation when silenced - measured as the fluorescence
%readout of the cells log(after / before) silencing.
%inputs:
%EssentialitiesFile - filename and path of the 
%Output - desired filename and path for generating the scaled essentialities file that will
%			be consumed by PrintTable.

load(EssentialitiesFile);
Essens = E.Values;
Essens = abs(Essens);

for i=1:size(Essens,2)
    temp = Essens(:,i);
    temp(find(isnan(temp)))=[];
    %temp(find(temp==0))=[];
    temp = abs(temp);
    Essens(:,i) = abs(Essens(:,i));
    Essens(:,i) = (Essens(:,i)-min(temp))/(max(temp)-min(temp));
end

Essentialities.Values = Essens;
Essentialities.Lines = E.Lines;
Essentialities.Histology = E.Histology;
Essentialities.Labels = E.Labels;
Essentialities.Source = E.Source;
Essentialities.Target = E.Target;

save(Output,Essentialities)

end

