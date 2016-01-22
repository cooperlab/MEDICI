function  ExtractTemplateReactions(Superpathway_HGNC_EntrezFile,HUGOFile,Output )

% Superpathway_HGNC_EntrezFile ='Superpathway.HGNC-Entrez.xlsx';
[num,txt,raw] = xlsread(Superpathway_HGNC_EntrezFile,1,'I2:T15716');
text = txt(:,[1 4 12]);

TargetType = text(:,2);
TemplateReaction = StringMatch({'TemplateReaction'}, TargetType);
other = text;
other(cell2mat(TemplateReaction),:) = [];
text = text(cell2mat(TemplateReaction),:);

[E,c] = size(text);

source = cell(E,1);
target = cell(E,1);
TRs = cell(E,1);

j = 1;

for i=1:E
    m_c = strsplit(char(text(i,1)),{'/',','});
    m_t = strsplit(char(text(i,3)),{'/',','});
    
    for t = 1:length(m_t)
        for n = 1:length(m_c)
            source{j} = m_c{n};
            target{j} = m_t{t};
            TRs{j} = [source{j} target{j}];
            j = j + 1;
        end
    end
    
end



[E,c] = size(other);

other_source = cell(E,1);
other_target = cell(E,1);
others = cell(E,1);

j = 1;

for i=1:E
    m_c = strsplit(char(other(i,1)),{'/',','});
    m_t = strsplit(char(other(i,3)),{'/',','});
    
    for t = 1:length(m_t)
        for n = 1:length(m_c)
            other_source{j} = m_c{n};
            other_target{j} = m_t{t};
            others{j} = [other_source{j} other_target{j}];
            j = j + 1;
        end
    end
    
end

% HUGOFile = '/Users/sharati/Google Drive/Bioinformatics/HUGO/hgnc_complete_set.txt';
%build HUGO structure
HUGO = ParseHUGO(HUGOFile);

source = HUGOLookup(source, HUGO, 'Symbol');
Missing = cellfun(@isempty, source);
source(Missing) = [];
target(Missing) = [];
TRs(Missing) = [];


target = HUGOLookup(target, HUGO, 'Symbol');
Missing = cellfun(@isempty, target);
source(Missing) = [];
target(Missing) = [];
TRs(Missing) = [];


Missing = find(cellfun(@(x)sum(strcmp(x, '')),others));
others(Missing) = [];
other_source(Missing) = [];
other_target(Missing) = [];


other_source = HUGOLookup(other_source, HUGO, 'Symbol');
Missing = cellfun(@isempty, other_source);
other_source(Missing) = [];
other_target(Missing) = [];
others(Missing) = [];


other_target = HUGOLookup(other_target, HUGO, 'Symbol');
Missing = cellfun(@isempty, other_target);
other_source(Missing) = [];
other_target(Missing) = [];
others(Missing) = [];

mapped = StringMatch(others,TRs);
m = cell2mat(mapped);
TRs(m) = [];
source(m) = [];
target(m) = [];


save(Output, 'source','target')

xlswrite('TemplateReaction',{'Source' 'Target'});
xlswrite('TemplateReaction',[source target],1,'A2');


end

