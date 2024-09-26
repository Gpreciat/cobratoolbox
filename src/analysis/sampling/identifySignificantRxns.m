function rankedNorms = identifySignificantRxns(model, samples, noOfRxnToCheck, exRxns, printLevel)
% Identifies significant reactions in genome-scale models by evaluating the
% impact of applying constraints through a three-step process: random
% flux space sampling, covariance calculation, and Euclidean distance
% computation. This function will help identify reactions with the greatest
% impact when constraints are applied to a genome-scale metabolic model.
%
% USAGE:
%
%    rankedNorms = identifySignificantRxns(model, samples, noOfRxnToCheck, exRxns, printLevel)
%
% INPUTS:
%    model:    A genome-scale COBRA model
%
%               * .S - Stoichiometric matrix
%               * .mets - Metabolite ID vector
%               * .rxns - Reaction ID vector
%               * .lb - Lower bound vector
%               * .ub - Upper bound vector
%
% OPTIONAL INPUTS:
%    samples:      `n x numSamples` matrix of flux vectors (calculate if it
%                  doesn't exist)
%    noOfRxnToCheck: The number of rxns to rank order (default: 20)
%    exRxns:       A flag to specify whether to consider only exchange
%                  reactions (default: false)
%    printLevel:   Verbosity level for printing (default: 2)
%
% OUTPUTS:
%    rankedNorms:    Table with ranked norms based on Euclidean distances
%
% EXAMPLE USAGE:
%
%    rankedNorms = identifySignificantRxns(model, samples, 20, true, 2)

if nargin < 2  || isempty(samples)
    sample = true;
else
    sample = false;
end

if nargin < 3 || isempty(noOfRxnToCheck)
    noOfRxnToCheck = 20;
end

if nargin < 4
    exRxns = false;    
end

if nargin < 5
    printLevel = 2;
end

% Sample the flux space using the CHRR algorithm
if sample
    try
        options.nStepsPerPoint = 1;
        options.nPointsReturned = 100;
        options.toRound = 1;
        [P, ~] =  sampleCbModel(model, [], [], options);
        options.nStepsPerPoint = 8 * size(P.A, 2);
        options.nPointsReturned = size(P.A, 2)^2;
        options.toRound = 0;
        [~, samples] =  sampleCbModel(model, [], [], options);
    catch
    end
end

% Calculate the covariance of the sampled space
covS = cov(samples');

% Check if focus in exRxn or not
if exRxns
    model = findSExRxnInd(model);
    flags = model.SIntRxnBool;
else
    flags = zeros(size(model.rxns));
end

% Calculate the norms
warning_state = warning('off', 'all');
if printLevel > 1
    display('Computing norms...')
    display(' ')
end

if printLevel > 1
    [S_rxns, S_adjusted_norms, S_norms] = getTopKNorms_adaptive(covS, noOfRxnToCheck, flags, printLevel);
else
    [S_rxns, S_adjusted_norms, S_norms] = getTopKNorms_adaptive(covS, noOfRxnToCheck, flags, 0);
end
plot_data = [S_adjusted_norms S_norms - S_adjusted_norms];
warning(warning_state);

% Create empty table
if length(unique(model.rxns(S_rxns))) < noOfRxnToCheck
    if printLevel > 1
        display(' ')
        display(['noOfRxnToCheck changed to ' num2str(length(unique(model.rxns(S_rxns))))])
    end
    noOfRxnToCheck = length(unique(model.rxns(S_rxns)));
end
nRows = noOfRxnToCheck;
varTypes = {'string', 'string', 'double', 'double'};
if exRxns
    varNames = {'mets', 'metNames', 'adjustedNorm', 'norm'};
else
    varNames = {'rxns', 'rxnNames', 'adjustedNorm', 'norm'};
end
rankedNorms = table('Size', [nRows length(varTypes)], 'VariableTypes', varTypes,...
    'VariableNames', varNames);

% Fill the table
if exRxns
    met_names = {};
    for j = 1:noOfRxnToCheck
        met_names{j} = model.metNames{find(full(model.S(:, strmatch(model.rxns{S_rxns(j)}, model.rxns))))};
        metID{j} = model.mets{find(full(model.S(:, strmatch(model.rxns{S_rxns(j)}, model.rxns))))};
    end
    [~, ia, ~] = unique(met_names);
    met_names = met_names(sort(ia));
    metID = metID(sort(ia));
    plot_data = plot_data(sort(ia), :);
    rankedNorms.mets = metID';
    rankedNorms.metNames = met_names';
else
    rankedNorms.rxns = model.rxns(S_rxns);
    rankedNorms.rxnNames = model.rxnNames(S_rxns);
end
rankedNorms.adjustedNorm = S_adjusted_norms(1:noOfRxnToCheck);
rankedNorms.norm = S_norms(1:noOfRxnToCheck);

% Make the figure
if printLevel > 0

    figure
    
    % Check if there are negative values use normal data istead log data
    logPlot_data = [log10(plot_data(:, 1))  log10(sum(plot_data, 2)) - log10(plot_data(:, 1))];
    if ~any(logPlot_data(:) < 0)
        plot_data = logPlot_data;
    end

    barh(plot_data,'stacked');
    if ~any(logPlot_data(:) < 0)
        xlabel('log Euclidean Norm', 'interpreter', 'latex')
    else
        xlabel('Euclidean Norm', 'interpreter', 'latex')
    end
    set(gca, 'yticklabels', S_rxns);
    yticks(1:noOfRxnToCheck)
    if exRxns
        title('Priority metabolites', 'FontSize', 14)
        set(gca, 'yticklabel', regexprep(rankedNorms.metNames,'\_', '\\_'), 'YDir', 'reverse', 'FontSize', 14);
    else
        title('Priority reactions', 'FontSize', 14)
        set(gca, 'yticklabel', regexprep(rankedNorms.rxnNames,'\_', '\\_'), 'YDir', 'reverse', 'FontSize', 14);
    end
end
