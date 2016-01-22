function ScaleEssentialities(EssentialitiesFile,Output)

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

