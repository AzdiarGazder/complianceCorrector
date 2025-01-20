function [E,density] = calcModulus(alloyElements,alloyComposition,varargin)

%% Pre-define options
% This script assumes composition is in weight percent or weight fraction
% unless specified as atomic percent or atomic fraction
flagAtomic = check_option(varargin,'atomic');

% Replace semicolons with commas if necessary to standardise the delimiters
alloyElements = strrep(alloyElements, ';', ',');
alloyElements = strtrim(strsplit(alloyElements, ','))';

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
for i = 1:length(alloyElements)
    idx = find(strcmp(data.symbol, alloyElements{i}));
    if ~isempty(idx)
        modulusList = [modulusList; data.modulus(idx)];
        atomicNumberList = [atomicNumberList; data.atomicNumber(idx)];
        massList = [massList; data.mass(idx)];
        densityList = [densityList; data.density(idx)];
    else
        warning('Element %s not found in the dataset.', alloyElements{i});
    end
end

if flagAtomic
    % Convert from atomic percent (at.%) to weight percent (wt.%)
    atomicFraction = alloyComposition / sum(alloyComposition); % convert to atomic fraction by normalising
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
    'Al'; 'Sb';
    'Ba'; 'Be'; 'Bi'; 'B';
    'Cd'; 'Ca'; 'C'; 'C_diamond'; 'Ce'; 'Cr'; 'Co'; 'Cu';
    'Dy';
    'Er'; 'Eu';
    'Fr';
    'Gd'; 'Ga'; 'Ge'; 'Au';
    'Hf'; 'He'; 'Ho'; 'H';
    'In'; 'I'; 'Ir';
    'Fe';
    'Kr';
    'La'; 'Pb'; 'Li'; 'Lu';
    'Mg'; 'Mn'; 'Hg'; 'Mo';
    'Nd'; 'Ne'; 'Ni'; 'Nb'; 'N';
    'Os';
    'Pd'; 'P'; 'Pt'; 'Pu'; 'Po'; 'K';
    'Ra'; 'Rn'; 'Re'; 'Rh'; 'Rb'; 'Ru';
    'Sm'; 'Sc'; 'Se'; 'Si'; 'Ag'; 'Na'; 'Sr'; 'S';
    'Ta'; 'Te'; 'Tb'; 'Tl'; 'Th'; 'Tm'; 'Sn'; 'Ti';
    'W';
    'U';
    'V';
    'Xe';
    'Yb'; 'Y';
    'Zn'; 'Zr'};

% Note: The elastic modulus of C is highly anisotropic and varies from 18
% to 40 GPa. Consequently, an average value of 29 GPa is used here.
elemental.modulus = [... % in GPa
    69; 55;
    13; 287; 32; 380;
    50; 20; 29; 1050; 33; 279; 209; 110;
    69;
    69; 18;
    nan;
    78; 9.8; 103; 78;
    78; nan; 64; nan;
    11; 11; 528;
    211;
    nan;
    36; 16; 5.5; 70;
    45; 198; nan; 329;
    41; nan; 200; 105; nan;
    550;
    121; 11; 168; 96; nan; 3.5;
    nan; nan; 463; 380; 2.4; 447;
    49; 74; 10; 130; 83; 10; 15; 7;
    186; 43; 55; 8; 79; 74; 50; 116;
    411;
    208;
    128;
    nan;
    24; 63;
    108; 68];

elemental.atomicNumber = [...
    13; 51;
    56; 4; 83; 5;
    48; 20; 6; 6; 58; 24; 27; 29;
    66;
    68; 63;
    87;
    64; 31; 32; 79;
    72; 2; 67; 1;
    49; 53; 77;
    26;
    36;
    57; 82; 3; 71;
    12; 25; 80; 42;
    60; 10; 28; 41; 7;
    76;
    46; 15; 78; 94; 84; 19;
    88; 86; 75; 45; 37; 44;
    62; 21; 34; 14; 47; 11; 38; 16;
    73; 52; 65; 81; 90; 69; 50; 22;
    74;
    92;
    23;
    54;
    70; 39;
    30; 40];

elemental.mass = [... % in g/mol
    26.98; 121.76;
    137.33; 9.0122; 208.98; 10.81;
    112.41; 40.08; 12.01; 12.01; 140.12; 52; 58.93; 63.55;
    162.5;
    167.26; 151.98;
    223;
    157.25; 69.72; 72.63; 196.97;
    178.49; 4.0026; 164.93; 1.008;
    114.82; 126.9; 192.22;
    55.85;
    83.8;
    138.91; 207.2; 6.94; 175;
    24.31; 54.94; 200.59; 95.95;
    144.24; 20.18; 58.69; 92.91; 14.01;
    190.23;
    106.42; 30.97; 195.08; 244.06; 209.98; 39.1;
    226.03; 222; 186.21; 102.91; 85.47; 101.07;
    150.36; 44.96; 78.96; 28.09; 107.87; 22.99; 87.62; 32.07;
    180.95; 127.6; 158.93; 204.38; 232.04; 168.93; 118.71; 47.87;
    183.84;
    238.03;
    50.94;
    131.29;
    173.04; 88.91;
    65.38; 91.22];

elemental.density = [... % in g/cm3
    2.7; 6.68;
    3.62; 1.848; 9.78; 2.34;
    8.65; 1.54; 2.267; 3.51; 6.77; 7.19; 8.9; 8.96;
    8.55;
    9.06; 5.24;
    nan;
    7.9; 5.91; 5.32; 19.32;
    13.31; 0.0001786; 8.8; 0.00008988;
    7.31; 4.93; 22.56;
    7.87;
    0.00375;
    6.15; 11.34; 0.534; 9.84;
    1.738; 7.43; 13.534; 10.28;
    7.01; 0.0008999; 8.9; 8.57; 0.0012506;
    22.59;
    12.03; 1.82; 21.45; 19.86; nan; 0.856;
    5.5; 0.00973; 21.02; 12.41; 1.532; 12.37;
    7.52; 2.98; 4.81; 2.33; 10.49; 0.968; 2.64; 2.07;
    16.69; 6.24; 8.23; 11.85; 11.72; 9.32; 7.31; 4.54;
    19.25;
    18.95;
    6.11;
    0.005887;
    6.9; 4.47;
    7.14; 6.52];

% % Reference list:
% % [1] Ashby, Shercliff, Cebon, "Materials: Engineering, Science, Processing".
% % [2] Callister, "Materials Science and Engineering: An Introduction".
% % [3] ASM International Handbook.
% % [4] CRC Handbook of Chemistry and Physics.
end