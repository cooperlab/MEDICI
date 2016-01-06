function HUGO = ParseHUGO(HUGOFile)
%Parses a HUGO gene symbol file describing HUGO gene symbols, synonyms, past
%symbols, and accession numbers. The HUGO file is obtained from the
%downloads section of http://genenames.org
%inputs: 
%HUGOFile - file linking HUGO gene symbols with synonyms, past symbols, and
%           accession numbers.
%outputs:
%HUGO - structure containing list of approved gene symbols and indexed
%       lists of previous symbols, synonyms, accessions linked to these
%       approved symbols.  This structure is used with the HUGOLookup
%       function to determine if a given string is a gene symbol or associated 
%       with a gene symbol.

%open HUGO file
Contents = text2cell('hgnc_complete_set.txt', '\t');

%copy header, erase
Header = Contents(1,:);
Contents(1,:) = [];

%identify approved symbol column
Hit = find(strcmp(Header, 'Approved Symbol'));
if(isempty(Hit))
    error('HUGO File Error: cannot locate column "Approved Symbol"');
else
    Approved = Contents(:,Hit);
end

%identify Locus Group column - protein coding gene flag
Hit = find(strcmp(Header, 'Locus Group'));
if(isempty(Hit))
    error('HUGO File Error: cannot locate column "Locus Group"');
else
    Group = strcmp(Contents(:,Hit), 'protein-coding gene');
end

%identify status symbol column
Hit = find(strcmp(Header, 'Status'));
if(isempty(Hit))
    error('HUGO File Error: cannot locate column "Status"');
else
    Status = Contents(:,Hit);
end

%discard unapproved symbols
Legitimate = find(strcmp(Status, 'Approved'));
Approved = Approved(Legitimate);
Group = Group(Legitimate);
Contents = Contents(Legitimate, :);

%identify previous symbol column
Hit = find(strcmp(Header, 'Previous Symbols'));
if(isempty(Hit))
    error('HUGO File Error: cannot locate column "Previous Symbols"');
else
    Previous = Contents(:,Hit);
end

%identify synonyms column
Hit = find(strcmp(Header, 'Synonyms'));
if(isempty(Hit))
    error('HUGO File Error: cannot locate column "Synonyms"');
else
    Synonyms = Contents(:,Hit);
end

%identify accession numbers column
Hit = find(strcmp(Header, 'Accession Numbers'));
if(isempty(Hit))
    error('HUGO File Error: cannot locate column "Accession Numbers"');
else
    Accession = Contents(:,Hit);
end

%identify Entrez Gene ID column
Hit = find(strcmp(Header, 'Entrez Gene ID'));
if(isempty(Hit))
    error('HUGO File Error: cannot locate column "Entrez Gene ID"');
else
    Entrez = Contents(:,Hit);
end


%create linked list between previous symbols and approved symbols
PreviousSymbols = cell(1,length(Approved));
PreviousIndices = cell(1,length(Approved));
for i = 1:length(Approved)
    PreviousSymbols{i} = DeLimit(Previous{i}, ',');
    PreviousIndices{i} = repmat(i, 1, length(PreviousSymbols{i}));
end
PreviousSymbols = [PreviousSymbols{:}].';
PreviousIndices = [PreviousIndices{:}].';

%create linked list between synonyms and approved symbols
SynonymSymbols = cell(1,length(Approved));
SynonymIndices = cell(1,length(Approved));
for i = 1:length(Approved)
    SynonymSymbols{i} = DeLimit(Synonyms{i}, ',');
    SynonymIndices{i} = repmat(i, 1, length(SynonymSymbols{i}));
end
SynonymSymbols = [SynonymSymbols{:}].';
SynonymIndices = [SynonymIndices{:}].';

%create linked list between accession numbers and approved symbols
AccessionSymbols = cell(1,length(Approved));
AccessionIndices = cell(1,length(Approved));
for i = 1:length(Approved)
    AccessionSymbols{i} = DeLimit(Accession{i}, ',');
    AccessionIndices{i} = repmat(i, 1, length(AccessionSymbols{i}));
end
AccessionSymbols = [AccessionSymbols{:}].';
AccessionIndices = [AccessionIndices{:}].';

%create linked list between Entrez Gene IDs and approved symbols
EntrezSymbols = cell(1,length(Approved));
EntrezIndices = cell(1,length(Approved));
for i = 1:length(Approved)
    EntrezSymbols{i} = DeLimit(Entrez{i}, ',');
    EntrezIndices{i} = repmat(i, 1, length(EntrezSymbols{i}));
end
EntrezSymbols = [EntrezSymbols{:}].';
EntrezIndices = [EntrezIndices{:}].';

%create structure list
ListStruct = struct('Symbols', {}, 'Indices', {});

%put lists into structure
HUGO = struct('Approved', {}, 'Group', {}, 'Previous', {}, 'Synonyms', {},...
                'Accessions', {}, 'Entrez', {});
HUGO(1).Approved = Approved;
HUGO(1).Group = Group;
ListStruct(1).Symbols = PreviousSymbols;
ListStruct(1).Indices = PreviousIndices;
HUGO(1).Previous = ListStruct;
ListStruct(1).Symbols = SynonymSymbols;
ListStruct(1).Indices = SynonymIndices;
HUGO(1).Synonyms = ListStruct;
ListStruct(1).Symbols = AccessionSymbols;
ListStruct(1).Indices = AccessionIndices;
HUGO(1).Accessions = ListStruct;
ListStruct(1).Symbols = EntrezSymbols;
ListStruct(1).Indices = EntrezIndices;
HUGO(1).Entrez = ListStruct;
end