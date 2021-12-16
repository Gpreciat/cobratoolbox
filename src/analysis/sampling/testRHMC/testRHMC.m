
global CBTDIR
if isempty(CBTDIR)
    initCobraToolbox
end

% Select a path
path = [CBTDIR filesep 'src' filesep 'analysis' filesep 'sampling' filesep 'testRHMC'];

% Load model
load([path filesep 'modelRHMC.mat'])

optionsSampling.samplerName = 'RHMC';
optionsSampling.nPointsReturned = 1000;
[~, S_model] =  sampleCbModel(model, [], optionsSampling.samplerName, optionsSampling);

if isempty(S_model)
    error('The sampling was not completed')
end
