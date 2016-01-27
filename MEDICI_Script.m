MAFFile = './Data/CCLE_hybrid_capture1650_hg19_damaging.maf';
LineFile = './Data/CCLE_sample_info_file_2012-10-18.txt';
CaptureOutputFile = './Data/CaptureFile.mat';

CNOutputFile = './Data/CNFile.mat';

SymbolFile =  './Data/PID.MSigDB.gmt.txt';
InteractionFile = './Data/PID.edge-attributes.txt';
PathwayListFile = './Data/PID.Pathways.txt';
HUGOFile = './Data/hgnc_complete_set.txt';
ScreenFile = './Data/ScreenedInteractions.mat';

PathwaysTableFile = './Data/PathwaysTableFile.mat';
superpathwayFile = './Data/FullSuperpathway.Strict.mat';
GeneEssentialitiesFile = './Data/Achilles.mat';
EssentialiesFile = './Data/Essentialities.mat';
Scaled_EssentialiesFile = './Data/Scaled_Essentialities.mat';
Final_EssentialiesFile = './Data/Final_Essentialies.mat';

GCTFile = './Data/Achilles.QCv2.4.3.rnai.gct';
CellDescriptionFile = './Data/Achilles.CellLineData.txt';
ESNullFile = './Data/ES.NullValues.mat';



TemplateReactionFile = './Data/TemplateReaction.mat';



alpha = 0.5;
w = 0.5;

CCLEHybridCapture(MAFFile, LineFile, CaptureOutputFile);

CN = CCLEBuildCNV('./Data/CCLE_copynumber_byGene_2013-12-03.txt', CNOutputFile);

GenerateSuperpathway(SymbolFile, InteractionFile, PathwayListFile, ...
    HUGOFile, ScreenFile, superpathwayFile);


PathwaysTable = CCLEBuildPathwaysTable(superpathwayFile,CNOutputFile,...
    CaptureOutputFile, PathwaysTableFile);

GeneCentricEssentiality(GCTFile, CellDescriptionFile, HUGOFile, [], ESNullFile, GeneEssentialitiesFile)

E = CalculateEssentialities(superpathwayFile,PathwaysTableFile,...
    GeneEssentialitiesFile,TemplateReactionFile,alpha,w, EssentialiesFile);

PrintTable(superpathwayFile,Scaled_EssentialiesFile,Final_EssentialiesFile);



