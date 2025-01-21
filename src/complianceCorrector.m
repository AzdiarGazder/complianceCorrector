function complianceCorrector(alloyElements,alloyComposition,varargin)
%% Function description:
% This script automatically corrects uniaxial tension/compression
% test data for compliance.
%
%% Author:
% Dr. Azdiar Gazder, 2025, azdiaratuowdotedudotau
%
%% Acknowledgements:
% Dr C.B. Finfrock
% For developing the original "correctcompliance.m" script, rev. Sep 2021.
%
% For details, please refer to the following reference:
% CB Finfrock, 'Temperature and strain rate dependence of the martensitic
% transformation and mechanical properties in advanced high strength
% s teels', Colorado School of Mines, PhD thesis, 2022.
% https://www.researchgate.net/publication/360627771_Temperature_and_Strain_Rate_Dependence_of_the_Martensitic_Transformation_and_Mechanical_Properties_in_Advanced_High_Strength_Steels
% https://www.researchgate.net/post/How-to-do-compliance-correction-for-a-given-stress-strain-data
% https://data.niaid.nih.gov/resources?id=mendeley_4p6w99kxg5
%
%% Syntax:
%  complianceCorrector()
%
%% Input:
%  alloyElements    - @char, a comma-separated list of the chemical 
%                     symbols of the elements constituting the alloy
%  alloyComposition - @double, a comma-separated list of the chemical 
%                     composition of the elements constituting the alloy.
%                     Alloy composition is assumed to be in weight percent 
%                     or weight fraction. Use the flag 'atomic' if the 
%                     alloy composition is in atomic percent or atomic 
%                     fraction.
%
%% Output:
%  *.txt            - a text file containing the compliance corrected
%                     time, displacement, force, engineering strain and
%                     engineering stress data
% Figures comprising:
% - Plot: Uncorrected and compliance corrected force vs displacement
% - Plot: Compliance corrected eng. stress vs eng. strain
%
%% Options:
% 'totalLength'     - @double, defines the length of the dog-bone sample.
% 'gageLength'      - @double, defines the parallel gage length of the
%                     dog-bone sample.
% 'width'           - @double, defines the width of the dog-bone sample 
%                     in the parallel gage length region.
% 'thickness'       - @double, defines the thickness of the dog-bone sample 
%                     in the parallel gage length region.
%
%%



%% Pre-define options
% If these values are specified, they are subsequently used for the 
% calculations. However, if these variables are empty, then the values in 
% the *.txt file are used.
totalLength_mm = get_option(varargin,'totalLength',[]); % total length (in mm)
gageLength_mm = get_option(varargin,'gageLength',[]); % gage length (in mm)
width_mm = get_option(varargin,'width',[]); % width (in mm)
thickness_mm = get_option(varargin,'thickness',[]); % thickness (in mm)
%% 



%% Default directories - Do not modify
iniDir = pwd;
Ini.dataPath = [strrep(iniDir,'\','/'),'/data/'];
Ini.inputPath = [Ini.dataPath,'input/'];
Ini.outputPath = [Ini.dataPath,'output/'];
%%



%% Load the stage data
[fileName, pathName] = uigetfile([Ini.inputPath, '*.txt'], 'Load the stage data');
if fileName == 0
    error('The program was terminated by the user');
    return;

else
    tic
    disp('...');
    disp('Loading stage data...');
    pfName = [pathName fileName];
    disp(pfName);

    %% Read the width, thickness, gage, and total lengths
    % Read the sample dimension data from the tab delimited text file
    opts1 = detectImportOptions(pfName, 'Delimiter', ':');
    % Extract only the necessary rows by reading between lines 8 to 11
    opts1.DataLines = [8 11];   %

    % Load the data as a table
    warning off;
    dataTable1 = readtable(pfName, opts1);
    warning on;

    % Extract the variables from the appropriate rows
    if isempty(width_mm)
        width_mm = dataTable1{1, 2}; % width (in mm)
    end
    if isempty(totalLength_mm)
        totalLength_mm = dataTable1{2, 2}; % total length (in mm)
    end
    if isempty(gageLength_mm)
        gageLength_mm = dataTable1{3, 2}; % gage length (in mm)
    end
    if isempty(thickness_mm)
        thickness_mm = dataTable1{4, 2}; % thickness (in mm)
    end


    % Read the uniaxial tension test data from row 14 onwards
    opts2 = detectImportOptions(pfName, 'Delimiter', '\t');
    % Start reading from the 14th row onwards
    opts2.DataLines = [14 Inf];
    % 13th row contains variable names
    opts2.VariableNamesLine = 13;
    % Load the data as a table
    warning off;
    dataTable2 = readtable(pfName, opts2);
    warning on;
    % Convert the table to a matrix
    dataMatrix2 = table2array(dataTable2);

    % Convert the data to the units as specified
    t = dataMatrix2(:, 1);            % time (in seconds)
    d = dataMatrix2(:, 4) * 10^-3;    % crosshead displacement (in mm)
    f = dataMatrix2(:, 5) * 10^-3;    % force (in kN)

    disp('Finished loading stage data...');
    toc
    disp('...');

    % Display the extracted data
    disp('...');
    disp('Extracted sample dimensions:');
    disp(['Width        = ', num2str(width_mm), ' mm']);
    disp(['Total length = ', num2str(totalLength_mm), ' mm']);
    disp(['Gage length  = ', num2str(gageLength_mm), ' mm']);
    disp(['Thickness    = ', num2str(thickness_mm), ' mm']);
    disp('...');

    % Calculate the initial cross-sectional area
    csArea_mm2 = width_mm * thickness_mm; % width * thickness, in mm2
end
%%



%% Calculate the theoretical elastic modulus of the alloy
[E,~] = calcModulus(alloyElements,alloyComposition,varargin{:});
elasticModulus_GPa = E.average; % in GPa
%%


%% Define the region-of-interest from test start to just at failure
uiwait(helpdlg({'LEFT-click, drag & release = Select a ROI from test start to just at failure';...
    'ENTER = when selection completed'}));

figure;
plot(t,d,'.k');
xlabel('Time (s)');
ylabel('Displacement (mm)');
roi1 = drawrectangle('color','r','lineWidth',0.5);
% Allow the user to resize and reposition the roi rectangle by forcing
% the pressing of any key to continue
pause;

% Select the data within the region of interest
tf1 = inROI(roi1,t,d);
close all;
t = t(tf1 == 1);
d = d(tf1 == 1);
f = f(tf1 == 1);

t_s = t - t(1); % zeroed time
d_mm = d - d(1); % zeroed displacement
f_kN = f; % force
%%


%% Perform least squares-based piecewise linear regression modelling on the
% time-displacement data to reduce noise in the displacement data
disp('...');
disp('Performing PLRM on stage displacement vs. time data...');
tic
[~,d_mm] = calcPLRM(t,d_mm,3,...
    'xLabel','Time (s)',...
    'yLabel','Displacement (mm)');
disp('Finished PLRM on stage displacement vs. time data...');
toc
disp('...')
pause;
close all;
%%


%% Define the region-of-interest from actual test start to just before
% failure. This step is needed to remove any slack at test start.
uiwait(helpdlg({'LEFT-click, drag & release = Select a ROI without the slack at test start to failure';...
    'ENTER = when selection completed'}));

figure;
plot(d_mm,f_kN,'.k');
xlabel('Displacement (mm)');
ylabel('Force (kN)');
roi2 = drawrectangle('color','g','lineWidth',0.5);
% Allow the user to resize and reposition the roi rectangle by forcing
% the pressing of any key to continue
pause;

% Select the data within the region of interest
tf2 = inROI(roi2,d_mm,f_kN);
close all;
t_s = t_s(tf2 == 1);
d_mm = d_mm(tf2 == 1);
f_kN = f_kN(tf2 == 1);
%%


%% Define the region-of-interest for the elastic region
uiwait(helpdlg({'LEFT-click, drag & release = Select a ROI defining the elastic region';...
    'ENTER = when selection completed'}));

figure;
plot(d_mm, f_kN,'.k');
xlabel('Displacement (mm)');
ylabel('Force (kN)');
roi3 = drawrectangle('color','b','lineWidth',0.5);
% Allow the user to resize and reposition the roi rectangle by forcing
% the pressing of any key to continue
pause;

% Select the data within the region of interest
tf3 = inROI(roi3,d_mm,f_kN);
close all;
d_mm_elastic = d_mm(tf3 == 1);
f_kN_elastic = f_kN(tf3 == 1);

% Fit the effective compliance of the system (sample + load frame)
% The unit for compliance is kN/mm
k_effective = polyfit(d_mm_elastic, f_kN_elastic, 1);

% Check the quality of the linear fit for the elastic region
f_kN_elasticFit = polyval(k_effective, d_mm_elastic);
figure;
plot(d_mm_elastic, f_kN_elastic, '.k')
hold all;
plot(d_mm_elastic, f_kN_elasticFit, '.-r')
xlabel('Displacement (mm)');
ylabel('Force (kN)');
legend('Elastic', 'Fitted', 'Location','southeast');
legend('boxoff');
hold off;

% Re-calculate the force data values based on the linear fit
fPlastic = f_kN(tf3 ~= 1);
fPlastic = fPlastic - fPlastic(1); % zero the plastic force
f_kN_plastic = f_kN_elasticFit(end) + fPlastic;
f_kN_corrected = [f_kN_elasticFit; f_kN_plastic];

% Re-calculate the displacement data values
dPlastic = d_mm(tf3 ~= 1);
dPlastic = dPlastic - dPlastic(1); % zero the plastic displacement
d_mm_plastic = d_mm_elastic(end) + dPlastic;
d_mm_corrected = [d_mm_elastic; d_mm_plastic];

% Use the springs-in-series reciprocal addition relationship to correct
% for compliance
% 1/k_effective = 1/k_sample + 1/k_loadFrame % (in  kN/mm)
k_effective(:,2) = [];
k_sample = (elasticModulus_GPa * csArea_mm2) / gageLength_mm;
k_loadFrame = 1 / ((1/k_effective) - (1/k_sample));

% Display the calculated compliance
disp('...');
disp('Calculated compliance:');
disp(['Sample compliance    = ', num2str(k_sample), ' kN/mm']);
disp(['Loadframe compliance = ', num2str(k_loadFrame), ' kN/mm']);
disp(['Effective compliance = ', num2str(k_effective), ' kN/mm']);
disp('...');

% Correct the displacement based on the calculated compliance values
d_mm_corrected = d_mm_corrected - (f_kN_corrected ./ k_loadFrame);
d_mm_corrected = d_mm_corrected - d_mm_corrected(1,1); % zero corrected displacement
% Center the displacement start at zero for the corrected data
startForce = f_kN_corrected(1,1);
startDisplacement = startForce / k_sample;
d_mm_corrected = d_mm_corrected + startDisplacement;

d_mm_corrected(1,1) = 0;
f_kN_corrected(1,1) = 0;
t_s = t_s - t_s(1);
%%


%% Plot the uncorrected and corrected data
figure;
plot(d_mm, f_kN,'.-k');
hold all;
plot(d_mm_corrected, f_kN_corrected, '.-r');
xlabel('Displacement (mm)');
ylabel('Force (kN)');
legend('Uncorrected', 'Corrected', 'Location','southeast');
legend('boxoff');
hold off;


%% Plot the corrected engineering stress-strain data
engStrain = d_mm_corrected./gageLength_mm;
engStress = (f_kN_corrected.*10^3)./csArea_mm2; % (in MPa)
figure;
plot(engStrain, engStress, '.-r');
xlabel('Eng. strain');
ylabel('Eng. stress (MPa)');
legend('Corrected', 'Location','southeast');
legend('boxoff');
hold off;


%% Define the paths to save data
outputSubfolder = fullfile(Ini.outputPath,fileName(1:end-4));
if ~exist(outputSubfolder, 'dir')
    mkdir(outputSubfolder);
end
pfName_dataOutput = [outputSubfolder,'/',fileName(1:end-4),'_corrected.txt'];

%% Save the *.txt data file
tic
disp('...');
disp('Saving compliance corrected test data...')

fileHeader1 = {'Alloy elements = ',  alloyElements,...
    'Alloy composition = ',  strrep(num2str(alloyComposition), '  ', ','),...
    'Theoretical elastic modulus (GPa) = ',  elasticModulus_GPa,...
    '----------',...
    'Length (mm) = ', num2str(gageLength_mm),...
    'Width (mm) = ', num2str(width_mm),...
    'Thickness (mm) = ', num2str(thickness_mm),...
    '----------'};

fileHeader2 = {'Time (s)',...
    'Displacement (mm)',...
    'Force (N)',...
    'Eng. strain',...
    'Eng. stress (MPa)'};

fileData = [t_s,...
    d_mm_corrected,...
    (f_kN_corrected.*10^3),...
    engStrain,...
    engStress]';

fid = fopen(pfName_dataOutput,'wt');
fprintf(fid,'%s%s\t\n%s%s\t\n%s%s\t\n%s\t\n%s%s\t\n%s%s\t\n%s%s\t\n%s\t\n',fileHeader1{:});
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t\n',fileHeader2{:});
fprintf(fid,'%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t\n',fileData);
fclose(fid);
toc
disp('...');
%%
