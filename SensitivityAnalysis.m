function SensitivityAnalysis(ScaledEssenFile,DrugRespondFile,Output)

load(ScaledEssenFile);
lines = Essentialities.Line;

load(DrugRespondFile)
ccle_lines = CCLE.Lines;

mapped = StringMatch(lines,ccle_lines);
ccle_AUC = CCLE.AUC(:,cell2mat(mapped));


[G,N] = size(Essentialities.Values);
correl = zeros(G,length(CCLE.Compounds));
pval = zeros(G,length(CCLE.Compounds));
n_cellLines = zeros(G,length(CCLE.Compounds));

for i = 1:G
    for j=1: length(CCLE.Compounds)
        e = Essentialities.Values(i,:)';
        auc = ccle_AUC(j,:)';
        nans = find(isnan(e));
        e(nans) = [];
        auc(nans) = [];
        
        nans = find(isnan(auc));
        e(nans) = [];
        auc(nans) = [];
               
        [RHO,PVAL] = corr(e,auc);
        correl(i,j) = RHO;
        pval(i,j) = PVAL;
        n_cellLines(i,j) = length(e);
    end
end

siz = size(ccle_AUC,1) * length(Target);
r_n_cellLines = reshape(n_cellLines,[1,siz]);
r_correl = reshape(correl,[1,siz]);
r_pval = reshape(pval,[1,siz]);

r_Target = Target';
r_Source = Source';

for i=1:(length(CCLE.Compounds)-1)
   r_Target = [r_Target,Target'];
   r_Source = [r_Source,Source'];
end

ppi = strcat(r_Source,'-',r_Target);

r_drug = {};

for i = 1:length(CCLE.Compounds)
    d = cell(length(Target),1);
    [d{:}] = deal(CCLE.Compounds{i});
    r_drug = [r_drug,d'];
end

T = table(ppi' ,r_drug',r_correl',r_pval',r_n_cellLines');
writetable(T,Output,'Delimiter','\t')

end
