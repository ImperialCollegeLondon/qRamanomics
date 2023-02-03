%% the calibrationCurve function performs data preprocessing steps for all input Raman spectrums

function [result] = calibrationCurve(matName)

% This function loads WITec .mat files of single spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the file into variable S, matName is a string that does not have .mat
filename = strcat(matName, '.mat');
S = load(filename);

%names will store all the field names for the variable S
names = fieldnames(S);

%xAxis will store the wavenumbers for the file
xAxis = getXAxis(S, names);

% Pull Raman intensity values matrix from exported WITec .mat file
yData = S.(names{1}).data(1, :);

% ModPoly contains all preprocessing steps for the input Raman spectrums
[waveNumber,Raman] = ModPoly(xAxis',yData');
intensity = max(Raman);

% result = matrix to be exported to workspace
result = intensity;

end
