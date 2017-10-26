function charge = reactionCharge(rxnFile)
% Calculates a reaction charge based on its RXN file
%
% USAGE:
%
%    charge = reactionCharge(rxnFile)
%
% INPUTS:
%    rxnFile:        Name and path of the RXN file to calculate charge.
%
% OUTPUTS:
%    charge:         Charge of the reaction.
%
% EXAMPLE:
%
%    charge = reactionCharge(rxnFile)
%
% .. Author: - German A. Preciat Gonzalez 25/05/2017

if exist(rxnFile, 'file') == 2
    
    rxnFile = regexp( fileread(rxnFile), '\n', 'split')';
    charges = strmatch('M  CHG', rxnFile);
    
    begmol = strmatch('$MOL', rxnFile);
    substrates = str2double(rxnFile{5}(1:3));
    products = str2double(rxnFile{5}(4:6));
    
    if ~isempty(charges)
        chargeS = 0;
        chargeP = 0;
        for i = 1:length(charges)
            chargeL = strsplit(strtrim(rxnFile{charges(i)}));
            if charges(i) < begmol(substrates + 1)
                for j = 5:2:length(chargeL)
                    chargeS = chargeS + str2double(chargeL(j));
                end
            else
                for j = 5:2:length(chargeL)
                    chargeP = chargeP + str2double(chargeL(j));
                end
            end
        end
        charge = chargeS - chargeP;
    else
        charge = 0;
    end
    
else
    charge = NaN;
end
