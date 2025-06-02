function [E,density] = calcModulus(alloyElements,alloyComposition,varargin)

%% Pre-define options
% This script assumes composition is in weight percent or weight fraction
% unless specified as atomic percent or atomic fraction
flagAtomic = check_option(varargin,'atomic'); 

if isa(alloyElements,'cell')
    alloyElements = alloyElements';

elseif isa(alloyElements,'char')
    % Replace semicolons with commas if necessary to standardise the delimiters
    alloyElements = strrep(alloyElements, '; ', ',');
    alloyElements = strtrim(strsplit(alloyElements, ','))';

else
    error('alloyElements must be of class cell or char.');
end

% Check if composition is a row vector and convert it to a column vector if necessary
if size(alloyComposition, 2) > size(alloyComposition, 1)
    alloyComposition = alloyComposition'; 
end

if size(alloyElements,1) ~= size(alloyComposition,1)
    error('The sizes of the elements and composition arrays do not match.')
end

if sum(alloyComposition / sum(alloyComposition)) ~= 1
    warning('The composition will be normalised as it does not sum to 1.')
end


% Get the reference elemental data
data = getElementData; 


% Initialise arrays to store the elemental properties
modulusList = []; 
atomicNumberList = []; 
massList = []; 
densityList = []; 

% Loop through each alloy symbol to collect the necessary elemental properties
for ii = 1:length(alloyElements)
    idx = find(strcmpi(data.symbol, alloyElements{ii})); 
    if ~isempty(idx)
        modulusList = [modulusList;  data.modulus(idx)]; 
        atomicNumberList = [atomicNumberList;  data.atomicNumber(idx)]; 
        massList = [massList;  data.mass(idx)]; 
        densityList = [densityList;  data.density(idx)]; 
    else
        warning('Element %s not found in the dataset.', alloyElements{ii}); 
    end
end

% Check for and delete nan values from the elemental data
nanCheck = isnan(modulusList) | isnan(densityList);
if any(nanCheck)
    disp('...'); 
    disp('Elements excluded from calculation (NaNs):')
    elements2Exclude = char(alloyElements(nanCheck));
    disp(elements2Exclude);

    % Remove nan elemental data from all variables
    alloyComposition = alloyComposition(~nanCheck);
    modulusList = modulusList(~nanCheck);
    atomicNumberList = atomicNumberList(~nanCheck);
    massList = massList(~nanCheck);
    densityList = densityList(~nanCheck);
end


if flagAtomic
    % Convert from atomic percent (at.%) to weight percent (wt.%)
    atomicFraction = alloyComposition / sum(alloyComposition);  % convert to atomic fraction by normalising
    weightFraction = (atomicFraction .* massList .* densityList) / ...
        sum(atomicFraction .* massList .* densityList); 
else
    % Convert from weight percent (wt.%) to weight fraction by normalising
    weightFraction = alloyComposition / sum(alloyComposition); 
end

% Convert weight fractions to volume fractions
volumeFraction = (weightFraction ./ (densityList .* massList)) ./ sum(weightFraction ./ (densityList .* massList)); 

% Calculate the Voigt (Rule of Mixtures, Reinforced),
% Reuss (Series) and  Voigt-Reuss-Hill (VRH) Average elastic moduli
E.voigt = sum(volumeFraction .* modulusList); 
E.reuss = 1 / sum(volumeFraction ./ modulusList); 
E.average = (E.voigt + E.reuss) / 2; 

% Calculate the alloy density using the inverse rule of mixtures
density.invRM = 1 / sum(weightFraction ./ densityList); 

% Calculate the alloy density using the weighted average
density.wtAvg = sum(volumeFraction .* densityList) / sum(volumeFraction); 

% Display the results
disp('...'); 
fprintf('Estimated theoretical elastic moduli:\n'); 
fprintf('Voigt modulus   = %.2f GPa\n', E.voigt); 
fprintf('Reuss modulus   = %.2f GPa\n', E.reuss); 
fprintf('Average modulus = %.2f GPa\n', E.average); 
disp('----'); 
fprintf('Estimated theoretical density:\n'); 
fprintf('Density (invRM) = %.4f g/cm3\n', density.invRM); 
fprintf('Density (wtAvg) = %.4f g/cm3\n', density.wtAvg); 
disp('...'); 
end





