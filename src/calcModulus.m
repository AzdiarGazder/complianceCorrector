function [E,density] = calcModulus(alloyElements,alloyComposition,varargin)

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
data = referenceElementalData; 


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




function elemental = referenceElementalData
elemental.symbol = {...
    'H'; 'He'; 'Li'; 'Be'; 'B'; 'C'; 'N'; 'O'; 'F'; 'Ne'; 
    'Na'; 'Mg'; 'Al'; 'Si'; 'P'; 'S'; 'Cl'; 'Ar'; 'K'; 'Ca'; 
    'Sc'; 'Ti'; 'V'; 'Cr'; 'Mn'; 'Fe'; 'Co'; 'Ni'; 'Cu'; 'Zn'; 
    'Ga'; 'Ge'; 'As'; 'Se'; 'Br'; 'Kr'; 'Rb'; 'Sr'; 'Y'; 'Zr'; 
    'Nb'; 'Mo'; 'Tc'; 'Ru'; 'Rh'; 'Pd'; ' Ag'; 'Cd'; 'In'; 'Sn'; 
    'Sb'; 'Te'; 'I'; 'Xe'; 'Cs'; 'Ba'; 'La'; 'Ce'; 'Pr'; 'Nd'; 
    'Pm'; 'Sm'; 'Eu'; 'Gd'; 'Tb'; 'Dy'; 'Ho'; 'Er'; 'Tm'; 'Yb'; 
    'Lu'; 'Hf'; 'Ta'; 'W'; 'Re'; 'Os'; 'Ir'; 'Pt'; 'Au'; 'Hg'; 
    'Tl'; 'Pb'; 'Bi'; 'Po'; 'At'; 'Rn'; 'Fr'; 'Ra'; 'Ac'; 'Th'; 
    'Pa'; 'U'; 'Np'; 'Pu'; 'Am'; ...
    }; 


% Note: The elastic modulus of C is highly anisotropic and varies from 18
% to 40 GPa. Consequently, an average value of 29 GPa is used here.
elemental.modulus = [... % in GPa
    nan; nan; 4.9; 287; nan; 29; nan; nan; nan; nan; 
    10; 45; 70; 47; nan; nan; nan; nan; nan; 20; 
    74; 116; 128; 279; 198; 211; 209; 200; 130; 108; 
    nan; nan; 8; 10; nan; nan; 2.4; nan; 64; 67; 
    105; 329; nan; 447; 275; 121; 85; 50; 11; 50; 
    55; 43; nan; nan; 1.7; 13; 37; 34; 37; 41; 
    46; 50; 18; 55; 56; 61; 64; 70; 74; 24; 
    67; 78; 186; 411; 463; nan; 528; 168; 78; nan; 
    8; 16; 32; nan; nan; nan; nan; nan; nan; 79; 
    nan; 208; nan; 96; nan; ...
    ]; 


elemental.atomicNumber = [... % Z
    1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 
    11; 12; 13; 14; 15; 16; 17; 18; 19; 20; 
    21; 22; 23; 24; 25; 26; 27; 28; 29; 30; 
    31; 32; 33; 34; 35; 36; 37; 38; 39; 40; 
    41; 42; 43; 44; 45; 46; 47; 48; 49; 50; 
    51; 52; 53; 54; 55; 56; 57; 58; 59; 60; 
    61; 62; 63; 64; 65; 66; 67; 68; 69; 70; 
    71; 72; 73; 74; 75; 76; 77; 78; 79; 80; 
    81; 82; 83; 84; 85; 86; 87; 88; 89; 90; 
    91; 92; 93; 94; 95; ...
    ]; 


elemental.mass = [... % in g/mol
    1.008;  4.002602;  6.94;  9.0121831;  10.81;  12.011;  14.007;  15.999;  18.998403163;  20.1797; 
    22.98976928;  24.305;  26.9815385;  28.085;  30.973761998;  32.06;  35.45;  39.948;  39.0983;  40.078; 
    44.955908;  47.867;  50.9415;  51.9961;  54.938044;  55.845;  58.933194;  58.6934;  63.546;  65.38; 
    69.723;  72.63;  74.921595;  78.971;  79.904;  83.798;  85.4678;  87.62;  88.90584;  91.224; 
    92.90637;  95.95;  97;  101.07;  102.9055;  106.42;  107.8682;  112.414;  114.818;  118.71; 
    121.76;  127.6;  126.90447;  131.293;  132.90545196;  137.327;  138.90547;  140.116;  140.90766;  144.242; 
    145;  150.36;  151.964;  157.25;  158.92535;  162.5;  164.93033;  167.259;  168.93422;  173.045; 
    174.9668;  178.49;  180.94788;  183.84;  186.207;  190.23;  192.217;  195.084;  196.966569;  200.592; 
    204.38;  207.2;  208.9804;  209;  210;  222;  223;  226;  227;  232.0377; 
    231.03588;  238.02891;  237;  244;  243; ...
    ]; 


elemental.density = [... % in g/cm3
0.0000899; 0.0001785; 0.535; 1.848; 2.46; 2.26; 0.001251; 0.001429; 0.001696; 0.0009; 
0.968; 1.738; 2.7; 2.33; 1.823; 1.96; 0.003214; 0.001784; 0.856; 1.55; 
2.985; 4.507; 6.11; 7.19; 7.47; 7.874; 8.9; 8.908; 8.96; 7.14; 
5.904; 5.323; 5.727; 4.819; 3.12; 0.00375; 1.532; 2.63; 4.472; 6.511; 
8.57; 10.28; 11.5; 12.37; 12.45; 12.023; 10.49; 8.65; 7.31; 7.31; 
6.697; 6.24; 4.94; 0.0059; 1.879; 3.51; 6.146; 6.689; 6.64; 7.01; 
7.264; 7.353; 5.244; 7.901; 8.219; 8.551; 8.795; 9.066; 9.32; 6.57; 
9.841; 13.31; 16.65; 19.25; 21.02; 22.59; 22.56; 21.45; 19.3; 13.534; 
11.85; 11.34; 9.78; 9.196; nan; 0.00973; nan; 5; 10.07; 11.724; 
15.37; 19.05; 20.45; 19.816; 13.67; ...
];


% % Reference list:
% % [1] Ashby, Shercliff, Cebon, "Materials: Engineering, Science, Processing".
% % [2] Callister, "Materials Science and Engineering: An Introduction".
% % [3] ASM International Handbook.
% % [4] CRC Handbook of Chemistry and Physics.
end
