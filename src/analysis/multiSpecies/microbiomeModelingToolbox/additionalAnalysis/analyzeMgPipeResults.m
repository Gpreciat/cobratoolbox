function analyzeMgPipeResults(infoFilePath,resPath,varargin)
% This function takes simulation results generated by mgPipe as the input
% and determines which computed fluxes and reaction abundances are
% significantly different between groups. Also creates violin plots of the 
% simulation results. Requires a file with sample information (e.g., 
% disease group, age) for the microbiome models that were generated and 
% interrogated through mgPipe.
%
% USAGE:
%
%    analyzeMgPipeResults(infoFilePath,resPath,varargin)
%
% REQUIRED INPUTS:
% infoFilePath:         Path to text file or spreadsheet with information 
%                       on analyzed samples including group classification
%                       with sample IDs as rows
% resPath:              char with path of directory where simulation 
%                       results are saved
%
% OPTIONAL INPUTS
% statPath:             char with path of directory where results of
%                       statistical analysis are saved
% violinPath:           char with path of directory where violin plots are
%                       saved
% sampleGroupHeaders    list of one or more column headers in file with the
%                       sample information that should be analyzed 
%                       (e.g., disease status, age). If not provided, the
%                       second column will be used.
%
% .. Author: Almut Heinken, 12/2020

parser = inputParser();
parser.addRequired('infoFilePath', @ischar);
parser.addRequired('resPath', @ischar);
parser.addParameter('statPath', [pwd filesep 'Statistics'], @ischar);
parser.addParameter('violinPath', [pwd filesep 'ViolinPlots'], @ischar);
parser.addParameter('sampleGroupHeaders', '', @iscellstr);

parser.parse(infoFilePath, resPath, varargin{:});

infoFilePath = parser.Results.infoFilePath;
resPath = parser.Results.resPath;
statPath = parser.Results.statPath;
violinPath = parser.Results.violinPath;
sampleGroupHeaders = parser.Results.sampleGroupHeaders;

% create the folders
mkdir(statPath)
mkdir(violinPath)

% Read in the file with sample information
infoFile = readInputTableForPipeline(infoFilePath);

% get all spreadsheet files in results folder
dInfo = dir(resPath);
fileList={dInfo.name};
fileList=fileList';
fileList(~contains(fileList(:,1),{'.csv','.txt'}))=[];
fileList(contains(fileList(:,1),{'ModelStat'}))=[];

% analyze data in spreadsheets
for i=1:length(fileList)
    sampleData = readInputTableForPipeline([resPath filesep fileList{i}]);
    
    % merge columns for shadow price results
    if strcmp(sampleData{1,2},'Source')
        for j=2:size(sampleData,1)
            sampleData{j,1}=[sampleData{j,1} '_' sampleData{j,2}];
        end
        sampleData(:,2)=[];
    end
    if strcmp(sampleData{1,3},'Source')
        for j=2:size(sampleData,1)
            sampleData{j,2}=[sampleData{j,2} '_' sampleData{j,3}];
        end
        sampleData(:,2:3)=[];
    end
    if strcmp(sampleData{1,2},'Objective')
        for j=2:size(sampleData,1)
            sampleData{j,1}=[sampleData{j,1} '_' sampleData{j,2}];
        end
        sampleData(:,2)=[];
    end
    
    sampleData(1,2:end)=strrep(sampleData(1,2:end),'microbiota_model_samp_','');
    sampleData(1,2:end)=strrep(sampleData(1,2:end),'microbiota_model_diet_','');
    sampleData(1,2:end)=strrep(sampleData(1,2:end),'microbiota_model_pDiet_','');
    sampleData(1,2:end)=strrep(sampleData(1,2:end),'microbiota_model_rich_','');
    sampleData(1,2:end)=strrep(sampleData(1,2:end),'host_microbiota_model_samp_','');
    
    % remove entries not in data
    [C,IA]=intersect(infoFile(2:end,1),sampleData(1,2:end));
    if length(C)<length(sampleData(1,2:end))
        error('Some sample IDs are not found in the file with sample information!')
    end
    
    % perform statistical analysis
    if ~isempty(sampleGroupHeaders)
        for j=1:length(sampleGroupHeaders)
            [Statistics,significantFeatures] = performStatisticalAnalysis(sampleData,infoFile,'stratification',sampleGroupHeaders{j});
            
            % Print the results as a text file
            filename = strrep(fileList{i},'.csv','');
            filename = strrep(filename,'.txt','');
            writetable(cell2table(Statistics),[statPath filesep filename '_' sampleGroupHeaders{j} '_Statistics'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
            if size(significantFeatures,2)>1
                writetable(cell2table(significantFeatures),[statPath filesep filename '_' sampleGroupHeaders{j} '_SignificantFeatures'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
            end
            
            % create violin plots
            currentDir=pwd;
            cd(violinPath)
            
            % create violin plots for net uptake and secretion files
            if any(contains(fileList{i,1},{'net_uptake_fluxes.csv','net_secretion_fluxes.csv'}))
                makeViolinPlots(sampleData, infoFile, 'stratification',sampleGroupHeaders{j}, 'plottedFeature', filename, 'unit', 'mmol/person/day')
            end
            cd(currentDir)
        end
    else
        [Statistics,significantFeatures] = performStatisticalAnalysis(sampleData,infoFile);
        
        % Print the results as a text file
        filename = strrep(fileList{i},'.csv','');
        filename = strrep(filename,'.txt','');
        writetable(cell2table(Statistics),[statPath filesep filename '_Statistics'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
        if size(significantFeatures,2)>1
            writetable(cell2table(significantFeatures),[statPath filesep filename '_SignificantFeatures'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');
        end
        
        % create violin plots
        currentDir=pwd;
        cd(violinPath)
        
        % create violin plots for net uptake and secretion files
        if any(strcmp(fileList{i,1},{'net_uptake_fluxes.csv','net_secretion_fluxes.csv'}))
            makeViolinPlots(sampleData, infoFile, 'plottedFeature', filename, 'unit', 'mmol/person/day')
        end
        cd(currentDir)
    end
    
end
