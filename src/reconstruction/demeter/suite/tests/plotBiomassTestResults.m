function [notGrowing,Biomass_fluxes] = plotBiomassTestResults(refinedFolder, reconVersion, varargin)
% This function plots the growth of refined reconstructions and
% reports whether any reconstructions are incapable of biomass production.
% Optionally, draft reconstructions can be included.
%
% USAGE:
%
%    [notGrowing,Biomass_fluxes] = plotBiomassTestResults(refinedFolder, reconVersion, varargin)
%
%
% REQUIRED INPUTS
% refinedFolder             Folder with refined COBRA models generated by
%                           the refinement pipeline
% reconVersion              Name of the refined reconstruction resource
%
% OPTIONAL INPUTS
% testResultsFolder         Folder where the test results should be saved
% numWorkers                Number of workers in parallel pool (default: 0)
% translatedDraftsFolder    Folder with  translated draft COBRA models
%                           generated by KBase pipeline to analyze (will
%                           only be analyzed if folder is provided)
%
% OUTPUT
% notGrowing                List of IDs for refined reconstructions that
%                           cannot produce biomass on at least one condition
% Biomass_fluxes             Computed biomass production fluxes for each model
%
% .. Author:
%       - Almut Heinken, 09/2020

% Define default input parameters if not specified
parser = inputParser();
parser.addRequired('refinedFolder', @ischar);
parser.addRequired('reconVersion', @ischar);
parser.addParameter('testResultsFolder', [pwd filesep 'TestResults'], @ischar);
parser.addParameter('numWorkers', 0, @isnumeric);
parser.addParameter('translatedDraftsFolder', '', @ischar);

parser.parse(refinedFolder, reconVersion, varargin{:});

refinedFolder = parser.Results.refinedFolder;
testResultsFolder = parser.Results.testResultsFolder;
numWorkers = parser.Results.numWorkers;
reconVersion = parser.Results.reconVersion;
translatedDraftsFolder = parser.Results.translatedDraftsFolder;

mkdir(testResultsFolder)

tol=0.0000001;

notGrowing = {};
cnt=1;

% initialize COBRA Toolbox and parallel pool
global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER)
    initCobraToolbox
end
solver = CBT_LP_SOLVER;

if numWorkers > 0
    % with parallelization
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers)
    end
end
environment = getEnvironment();

if ~isempty(translatedDraftsFolder)
    % test draft and refined reconstructions
    folders={
        translatedDraftsFolder
        refinedFolder
        };
else
    % only refined reconstructions
    folders={
        refinedFolder
        };
end

for f=1:length(folders)
    dInfo = dir(folders{f});
    modelList={dInfo.name};
    modelList=modelList';
    modelList(~contains(modelList(:,1),'.mat'),:)=[];
    
   parfor i=1:length(modelList)
     %    for i=1:length(modelList)
        restoreEnvironment(environment);
        changeCobraSolver(solver, 'LP', 0, -1);
        try
            model=readCbModel([folders{f} filesep modelList{i}]);
        catch % circumvent the verifyModel for the moment
            model=  load([folders{f} filesep modelList{i}]);
            model = model.model;
        end
        biomassID=find(strncmp(model.rxns,'bio',3));
        [AerobicGrowth, AnaerobicGrowth] = testGrowth(model, model.rxns(biomassID));
        aerRes{i,f}=AerobicGrowth;
        anaerRes{i,f}=AnaerobicGrowth;
    end
    
    for i=1:length(modelList)
        growth{f}(i,1)=aerRes{i,f}(1,1);
        growth{f}(i,2)=anaerRes{i,f}(1,1);
        growth{f}(i,3)=aerRes{i,f}(1,2);
        growth{f}(i,4)=anaerRes{i,f}(1,2);
    end
end

data=[];
for f=1:length(folders)
    data(:,size(data,2)+1:size(data,2)+2)=growth{f}(:,1:2);
end

if ~isempty(translatedDraftsFolder)
    % draft and refined reconstructions
    % workaround if growth for drafts on anaerobic is all zeros
    try
        
        if sum(data(:,2))< tol
            data(1,2)=tol;
        end
        
        figure;
        hold on
        violinplot(data, {'Aerobic, Draft','Anaerobic, Draft','Aerobic, Refined','Anaerobic, Refined'});
        set(gca, 'FontSize', 12)
        maxval=max(data,[],'all');
        ylim([0 maxval + maxval/10])
        ylabel('mmol *g dry weight-1 * hr-1')
        h=title(['Growth on rich medium, ' reconVersion]);
        set(h,'interpreter','none')
        set(gca,'TickLabelInterpreter','none')
        print([testResultsFolder filesep 'Growth_rates_Rich_medium_' reconVersion],'-dpng','-r300')
        
        data=[];
        for f=1:length(folders)
            data(:,size(data,2)+1:size(data,2)+2)=growth{f}(:,3:4);
        end
        % workaround if growth for drafts is all zeros
        if sum(data(:,1))< tol
            data(1,1)=tol;
        end
        if sum(data(:,2))< tol
            data(1,2)=tol;
        end
        
        figure;
        hold on
        violinplot(data, {'Aerobic, Draft','Anaerobic, Draft','Aerobic, Refined','Anaerobic, Refined'});
        set(gca, 'FontSize', 12)
        maxval=max(data,[],'all');
        ylim([0 maxval + maxval/10])
        ylabel('mmol *g dry weight-1 * hr-1')
        h=title(['Growth on complex medium, ' reconVersion]);
        set(h,'interpreter','none')
        set(gca,'TickLabelInterpreter','none')
        print([testResultsFolder filesep 'Growth_rates_complex_medium_' reconVersion],'-dpng','-r300')
        
    end
    
    % report draft models that are unable to grow
    fprintf('Report for draft models:\n')
    noGrowth=growth{1}(:,1) < tol;
    if sum(noGrowth) > 0
        fprintf([num2str(sum(noGrowth)) ' models are unable to produce biomass on rich medium.\n'])
    else
        fprintf('All models are able to produce biomass on rich medium.\n')
    end
    
    noGrowth=growth{1}(:,2) < tol;
    if sum(noGrowth) > 0
        fprintf([num2str(sum(noGrowth)) ' models are unable to produce biomass on rich medium under anaerobic conditions.\n'])
    else
        fprintf('All models are able to produce biomass on rich medium under anaerobic conditions.\n')
    end
    
    noGrowth=growth{1}(:,3) < tol;
    if sum(noGrowth) > 0
        fprintf([num2str(sum(noGrowth)) ' models are unable to produce biomass on complex medium.\n'])
    else
        fprintf('All models are able to produce biomass on complex medium.\n')
    end
    
    noGrowth=growth{1}(:,4) < tol;
    if sum(noGrowth) > 0
        fprintf([num2str(sum(noGrowth)) ' models are unable to produce biomass on complex medium under anaerobic conditions.\n'])
    else
        fprintf('All models are able to produce biomass on complex medium under anaerobic conditions.\n')
    end
    
    % report refined models that are unable to grow
    fprintf('Report for refined models:\n')
    noGrowth=growth{2}(:,1) < tol;
    if sum(noGrowth) > 0
        fprintf([num2str(sum(noGrowth)) ' models are unable to produce biomass on rich medium.\n'])
        for i=1:length(noGrowth)
            if noGrowth(i)
                notGrowing{cnt,1}=modelList{i,1};
                cnt=cnt+1;
            end
        end
    else
        fprintf('All models are able to produce biomass on rich medium.\n')
    end
    
    noGrowth=growth{2}(:,2) < tol;
    if sum(noGrowth) > 0
        fprintf([num2str(sum(noGrowth)) ' models are unable to produce biomass on rich medium under anaerobic conditions.\n'])
        for i=1:length(noGrowth)
            if noGrowth(i)
                notGrowing{cnt,1}=modelList{i,1};
                cnt=cnt+1;
            end
        end
    else
        fprintf('All models are able to produce biomass on rich medium under anaerobic conditions.\n')
    end
    
    noGrowth=growth{2}(:,3) < tol;
    if sum(noGrowth) > 0
        fprintf([num2str(sum(noGrowth)) ' models are unable to produce biomass on complex medium.\n'])
        for i=1:length(noGrowth)
            if noGrowth(i)
                notGrowing{cnt,1}=modelList{i,1};
                cnt=cnt+1;
            end
        end
    else
        fprintf('All models are able to produce biomass on complex medium.\n')
    end
    
else
    % only refined reconstructions
    try
        figure;
        hold on
        violinplot(data, {'Aerobic','Anaerobic'});
        set(gca, 'FontSize', 12)
        maxval=max(data,[],'all');
        ylim([0 maxval + maxval/10])
        ylabel('mmol *g dry weight-1 * hr-1')
        h=title(['Growth on rich medium, ' reconVersion]);
        set(h,'interpreter','none')
        set(gca,'TickLabelInterpreter','none')
        print([testResultsFolder filesep 'Growth_rates_Rich_medium_' reconVersion],'-dpng','-r300')
    end
    
    data=[];
    for f=1:length(folders)
        data(:,size(data,2)+1:size(data,2)+2)=growth{f}(:,3:4);
    end
    
    if size(data,1)>5
        figure;
        hold on
        violinplot(data, {'Aerobic','Anaerobic'});
        set(gca, 'FontSize', 12)
        maxval=max(data,[],'all');
        ylim([0 maxval + maxval/10])
        ylabel('mmol *g dry weight-1 * hr-1')
        h=title(['Growth on complex medium, ' reconVersion]);
        set(h,'interpreter','none')
        set(gca,'TickLabelInterpreter','none')
        print([testResultsFolder filesep 'Growth_rates_complex_medium_' reconVersion],'-dpng','-r300')
    end
    
    % report refined models that are unable to grow
    fprintf('Report for refined models:\n')
    noGrowth=growth{1}(:,1) < tol;
    if sum(noGrowth) > 0
        fprintf([num2str(sum(noGrowth)) ' models are unable to produce biomass on rich medium.\n'])
        for i=1:length(noGrowth)
            if noGrowth(i)
                notGrowing{cnt,1}=modelList{i,1};
                cnt=cnt+1;
            end
        end
    else
        fprintf('All models are able to produce biomass on rich medium.\n')
    end
    
    noGrowth=growth{1}(:,2) < tol;
    if sum(noGrowth) > 0
        fprintf([num2str(sum(noGrowth)) ' models are unable to produce biomass on rich medium under anaerobic conditions.\n'])
        for i=1:length(noGrowth)
            if noGrowth(i)
                notGrowing{cnt,1}=modelList{i,1};
                cnt=cnt+1;
            end
        end
    else
        fprintf('All models are able to produce biomass on rich medium under anaerobic conditions.\n')
    end
    
    noGrowth=growth{1}(:,3) < tol;
    if sum(noGrowth) > 0
        fprintf([num2str(sum(noGrowth)) ' models are unable to produce biomass on complex medium.\n'])
        for i=1:length(noGrowth)
            if noGrowth(i)
                notGrowing{cnt,1}=modelList{i,1};
                cnt=cnt+1;
            end
        end
    else
        fprintf('All models are able to produce biomass on complex medium.\n')
    end
end

% export models that cannot grow on at least one condition
notGrowing=unique(notGrowing);
notGrowing=strrep(notGrowing,'.mat','');
if size(notGrowing,1)>0
    save([testResultsFolder filesep 'notGrowing.mat'],'notGrowing');
end

% export computed biomass  fluxes
if ~isempty(translatedDraftsFolder)
    data=growth{1}(:,1:4);
    data=[data,growth{2}(:,1:4)];
    Biomass_fluxes = {'','Unlimited aerobic, Draft','Unlimited anaerobic, Draft','Complex medium aerobic, Draft','Complex medium anaerobic, Draft','Unlimited aerobic, Refined','Unlimited anaerobic, Refined','Complex medium aerobic, Refined','Complex medium anaerobic, Refined'};
    Biomass_fluxes(2:length(modelList)+1,1) = strrep(modelList,'.mat','');
    Biomass_fluxes(2:end,2:9) = num2cell(data);
else
    data=growth{1}(:,1:4);
    Biomass_fluxes = {'','Unlimited aerobic, Refined','Unlimited anaerobic, Refined','Complex medium aerobic, Refined','Complex medium anaerobic, Refined'};
    Biomass_fluxes(2:length(modelList)+1,1) = strrep(modelList,'.mat','');
    Biomass_fluxes(2:end,2:5) = num2cell(data);
end

writetable(cell2table(Biomass_fluxes),[testResultsFolder filesep 'Biomass_fluxes.csv'],'WriteVariableNames',false)

end
