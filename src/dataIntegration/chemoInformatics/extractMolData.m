function molData = extractMolData(filename, molFilePath)
% Read a mol file and generate MOLdata structure with the molecule information
%
% USAGE:
%
%    molData = extractMOLdata(filename, molFilePath)
%
% INPUTS:
%    filename:      Name of the mol file.
%
% OPTIONAL INPUTS:
%    molFilePath    Path of the mol file.
%
% OUTPUTS:
%    MOLdata:       Struct array with the data of MDL MOL file with 
%                   following fields:%
%                       * .mName - Molecule name.
%                       * .InChiKey - InChiKey of the molecule.
%                       * .atoms - # of atoms in the molecule.
%                       * .bonds - # of bonds in the molecule.
%                       * .formula - Formula of the molecule.
%                       * .charge - Charge of the molecule.        
%
% EXAMPLE:
%
%    filename = 'ATP.mol'
%    molData = extractMOLdata(filename)
%
% .. Author: - German A. Preciat Gonzalez 04/04/2016

if nargin < 2 || isempty(molFilePath)
    molFilePath = [pwd filesep];        
else
    molFilePath = [regexprep(molFilePath,'(/|\\)$',''), filesep];
end              

% Read the MDL MOL file
molFile = regexp( fileread([molFilePath filename]), '\n', 'split')';

% The molecule name and the InChIKey
molData.mName = molFile{1};
molData.InChIKey = molFile{3};

% Number of atoms and number of bonds in the molecule
molData.atoms = str2double(molFile{4}(1:3));
molData.bonds = str2double(molFile{4}(4:6));

% Molecule's formula
for i = 1:molData.atoms
    atoms{i} = strtrim(molFile{4 + i}(32:33));
end
rGroups = strmatch('A', atoms, 'exact');
if ~isempty(rGroups)
    for i = 1:length(rGroups)
      atoms{rGroups(i)} = 'R';
    end
end
uniqueAtoms = unique(atoms);
formula = [];
for i = 1:length(uniqueAtoms)
    noOfAtoms = strmatch(uniqueAtoms{i}, atoms, 'exact');
    if length(noOfAtoms) ~= 1
        formula = [formula, uniqueAtoms{i} num2str(length(noOfAtoms))];
    else
        formula = [formula, uniqueAtoms{i}];
    end
end
molData.formula = formula;

% Molecule's charge
charges = strmatch('M  CHG', molFile);
if ~isempty(charges)
    charge = 0;
    for k = 1:length(charges)
        chargeL = strsplit(strtrim(molFile{charges(k)}));
        for j = 5:2:length(chargeL)
            charge = charge + str2double(chargeL(j));
        end
    end
else
    charge = 0;
end
molData.charge = charge;
