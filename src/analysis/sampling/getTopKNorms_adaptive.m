function [A_rxns, adjustedA_norms, A_norms] = getTopKNorms_adaptive(A, k, flags, printLevel)
% This function provides a heuristic approach to select a subset of 
% reactions with high covariance. This function aims to identify a subset 
% of reactions with the highest covariance in a given covariance matrix, 
% biasing the selection towards "independent" reactions. It is particularly 
% useful for genome-scale metabolic models when considering the impact of 
% applying constraints.
%
% USAGE:
%
%    [A_rxns, adjustedA_norms, A_norms] = getTopKNorms_adaptive(A, k, flags, printLevel)
%
% INPUTS:
%    A:          The covariance matrix of samples
%    k:          The number of reactions to select
%    flags:      A binary array indicating exchange reactions (1: 
%                non-exchange, 0: exchange), so we will not select any 
%                reaction with a 1.
%    printLevel:   Verbosity level for printing
%
% OUTPUTS:
%    A_rxns:     A list of indices corresponding to the top k reactions
%    adjustedA_norms: The corresponding adjusted norms of A_rxns
%    A_norms:    The corresponding norms of A_rxns
%
% EXAMPLE USAGE:
%
%    [A_rxns, A_norms] = selectTopIndependentReactions(covarianceMatrix, 10, exchangeFlags, 0)
%

if nargin < 3
    flags = zeros(size(A, 2), 1);
end
if nargin < 4
    printLevel = 0;
end

indices = [];
adjustedA_norms = [];
for i = 1:k

    A_modified = remove_indices(A, indices);
    A_norms = sqrt(sum(A_modified.^2, 2));
    [A_norms_sorted, A_index] = sort(A_norms,1, 'descend');
%     A_index_flagged = A_index(find(~flags));
    % A_norms_sorted_flagged = A_norms_sorted(find(~flags));
    
    j = 1;
    while flags(A_index(j)) == 1
       j = j + 1; 
    end
    A_rxn = A_index(j);
    % A_norms_sorted = A_norms_sorted_flagged(1);
    
    indices = [indices; A_rxn];
    adjustedA_norms = [adjustedA_norms; A_norms_sorted(j)];
    if printLevel > 0
        fprintf('Reaction %d/%d selected\n', i, k);
    end
end

A_rxns = indices;
A_norms = sqrt(sum(A.^2, 2));
A_norms = A_norms(indices);

end

function [A] = remove_indices(A, indices)
% This function removes rows specified by indices from matrix A.

if isempty(indices)
    return;
end
    
% Extract out the rows;
B = A(indices, :);

for i = 1:size(A, 2)
    x = (B' \ A(i, :)')';
    A(i, :) = A(i, :) - x * B;
end
    
end
