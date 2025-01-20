% *********************************************************************
%            complianceCorrector - Demonstration Example
% *********************************************************************
% This script demonstrates how to invoke the complianceCorrector function. 
% Save the compliance corrected uniaxial data in a sub-folder within the 
% ~/data/input directory.
% All output will be saved in a sub-folder of the same name in the 
% ~/data/output directory.
%
% *********************************************************************
% Dr. Azdiar Gazder, 2025, azdiaratuowdotedudotau
% (Remove "dot" and "at" to make this email address valid)
% *********************************************************************

%% Clear variables
home; clc; clear all; clear hidden; close all;
currentFolder;
warning off MATLAB:subscripting:noSubscriptsSpecified
set(0,'DefaultFigureWindowStyle','normal');


%% Call the complianceCorrector
% Use the sample dimensions in the *.txt file
% Here alloy composition is assumed to be in weight percent or weight
% fraction
complianceCorrector('Zr, Ti, Nb',[35, 40, 25]);

% % Specify a user-defined set of sample dimensions
% % when a sample dimension is not specified, the value in the *.txt file 
% % is used
% complianceCorrector('Zr, Ti, Nb',[35, 40, 25],...
%     'width',2.2,...
%     'thickness', 1.0,...
%     'gageLength',5.0);

% % Specify the alloy composition using atomic percent
% complianceCorrector('Zr, Ti, Nb',[35, 40, 25],...
%     'atomic',...
%     'width',2.2);

