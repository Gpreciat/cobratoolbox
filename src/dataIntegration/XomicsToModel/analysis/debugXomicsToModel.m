function debugXomicsToModel(genericModel, directory, modelGenerationReport, coreData)
% Examine the debug files generated by the XomicsToModel function for
% metabolites, reactions, and genes of interest to determine whether or not
% they were removed during the XomicsToModel process.
%
% USAGE:
%
%    debugXomicsToModel(model, directory, contextSpecificData)
%
% INPUT:
%    genericModel:     The generic input COBRA model used in XomicsToModel
%
%       * .S - Stoichiometric matrix
%       * .mets - Metabolite ID vector
%       * .rxns - Reaction ID vector
%       * .lb - Lower bound vector
%       * .ub - Upper bound vector
%       * .genes - Upper bound vector
%
%   directory:  Folder where the debug files are located

%   coreData: List of metabolites, reactions and genes to examine
%   in order to determine whether or not they were removed during the
%   XomicToModel process:
%
%       * .genes - List of Entrez gene IDs to check
%       * .rxns - List of reactions IDs to check
%       * .mets - List of metabolites IDs to check
%

if isempty(modelGenerationReport)
    if isempty(coreData)
        error('both modelGenerationReport and coreData must not be empty')
    else
        if isfield(coreData, 'genes')
            % Genes in correct format
            if isnumeric(coreData.genes)
                format long g
                tmp = cell(length(coreData.genes), 1);
                for i = 1:length(coreData.genes)
                    tmp{i,1} = num2str(coreData.genes(i));
                end
                coreData.genes = tmp;
            end
            coreData.genes = coreData.genes;
        else
            coreData.genes = [];
        end
        if ~isfield(coreData, 'rxns')
            coreData.rxns = [];
        end
        if ~isfield(coreData, 'mets')
            coreData.mets = [];
        end
    end
else
    coreData.rxns = modelGenerationReport.coreRxnAbbr0;
    coreData.mets = modelGenerationReport.coreMetAbbr0;
    coreData.genes = modelGenerationReport.activeGeneID0;
end

% Find all debug files in specified directory
if ~exist(directory,'dir')
    msg = ['Directory ' directory ' does not exist'];
    error(msg)
end
debugFiles = dir([directory filesep '*.mat']);
debugFiles = {debugFiles.name}';
debugFiles(ismember(debugFiles, 'Model.mat')) = [];
debugFiles(~contains(debugFiles, 'debug')) = [];
[~, ia] = sort(str2double(regexp(debugFiles, '\d+', 'match', 'once')));
debugFiles = debugFiles(ia);

% Generic data
genericGenes = regexprep(genericModel.genes, '\.\d', '');
genericRxns = genericModel.rxns;
genericMets = genericModel.mets;

% Prit
% Active list
fprintf('%13s %s %13s %s %13s %s %s %s\n', '#Active_genes', '            ', ...
    '#Active_rxns', '           ',  '#Active_mets', '           ',  '  ', 'Stage')
fprintf('%13u %12s %13u %11s %13u %11s %s %s\n', length(coreData.genes), ...
    ' ', length(coreData.rxns), ' ', length(coreData.mets),  ' ',  '  ', 'Active list')

% Model data
if ~isempty(coreData.genes)
    boolGenes = ismember(coreData.genes, genericGenes);
else
    boolGenes = zeros(size(genericGenes));
end
if ~isempty(coreData.rxns)
    boolRxns = ismember(coreData.rxns, genericRxns);
else
    boolRxns = zeros(size(genericRxns));
end
if ~isempty(coreData.mets)
    boolMets = ismember(coreData.mets, genericMets);
else
    boolMets = zeros(size(genericMets));
end

% Check feasibility
sol = optimizeCbModel(genericModel);
if  sol.stat ~= 1
    message = 'Infeasible';
else
    message = 'Feasible';
end
if nnz(genericModel.c)==1
    lb_obj = genericModel.lb(genericModel.c~=0);
else
    lb_obj = NaN;
end
% Draft model data
draftGenes = regexprep(genericModel.genes, '\.\d', '');
draftRxns = genericModel.rxns;
draftMets = genericModel.mets;
    
fprintf('%13s %s %13s %s %13s %s %s %s %s %s %s\n', '#Active_genes', '#Model_genes', ...
    '#Active_rxns', '#Model_rxns', '#Active_mets', '#Model_mets',  '   lb_obj', '   Obj', ' Message', '  ', 'Stage')
fprintf('%13u %12u %13u %11u %13u %11u %8g %8g %s %s %s\n', nnz(boolGenes), length(draftGenes), ...
    nnz(boolRxns), length(draftRxns), nnz(boolMets), length(draftMets),  lb_obj, sol.f, message, '  ', 'Generic model')

% Draft data
for i = 1:length(debugFiles)
    
    % Load the model at each stage
    load([directory filesep debugFiles{i}], 'model');
    
    % Draft model data
    draftGenes = regexprep(model.genes, '\.\d', '');
    draftRxns = model.rxns;
    draftMets = model.mets;
    
    if ~isempty(coreData.genes)
        boolGenes = ismember(coreData.genes, draftGenes);
    else
        boolGenes = zeros(size(draftGenes));
    end
    
    if ~isempty(coreData.rxns)
        boolRxns = ismember(coreData.rxns, draftRxns);
    else
        boolRxns = zeros(size(genericRxns));
    end
    if ~isempty(coreData.mets)
        boolMets = ismember(coreData.mets, draftMets);
    else
        boolMets = zeros(size(genericMets));
    end
    
    % Check feasibility
    sol = optimizeCbModel(model);
    if  sol.stat ~= 1
        message = 'Infeasible';
    else
        message = 'Feasible';
    end
    
    if nnz(model.c)==1
        lb_obj = model.lb(model.c~=0);
    else
        lb_obj = NaN;
    end
    
    %printConstraints(model,-inf,inf,ismember(model.rxns,'biomass_reaction'))
   
    
    fprintf('%13u %12u %13u %11u %13u %11u %11g %11g %s %s %s\n', nnz(boolGenes), length(draftGenes), ...
        nnz(boolRxns), length(draftRxns), nnz(boolMets), length(draftMets),  lb_obj, sol.f, message, '  ', debugFiles{i})
    
end
