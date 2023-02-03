%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % LEAST SQUARES REGRESSION FOR RAMAN IMAGE DATASETS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Vernon LaLone
% Date: 23 March 2021

% Save Raman image files into same folder with copies of regression
% analysis.m, reference.mat, and Xaxis.mat

% If all has been saved/copied accurately as described, the user may simply
% insert name of sample dataset in line below and Run this
% regression_analysis.m script.

%Insert name of sample dataset here (do not include ".mat"):

% Desired results are output into coeffMatrix file which may be saved
% and/or exported for plotting in alternative software.

% More detailed description of steps outlined below:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Step 1: Export Data from WITec Project Software 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Export all files in matlab format (right click > export > Matlab) as Matlab Export Type DSO 6.0

% Ensure the XUnits are in Wavenumbers (rel. 1/cm)

% Save all files in appropriate folders accessible thru Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Step 2: Import Data to Matlab Workspace 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = final(matName)
%this function creates heatmaps of the tissue sections 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loads the specified file structure into the workspace as "S"
filename = strcat(matName, '.mat');
S = load(filename);

% names will store all the field names for the variable S
names = fieldnames(S);

% xAxis will store the wavenumbers for the file
load('xAxis.mat');
xAxis = getXAxis(S, names);

%%% Load data for model
load('reference.mat');
refMat=double(reference');

%%% Set xaxis
waveNumbers_ref=double(xAxis');

%acquire the dimensions of the image so that we can reshape it later
numRows = S.(names{1}).imagesize(1);
numCols = S.(names{1}).imagesize(2);

% removes column if the entire column is zero
% sampleMat(:,all(sampleMat == 0))=[];

% normalize all spectra to one
% sampleMat = normalize(sampleMat, 'range');

%get the number of pixels that are in the image
length = getLength(S, names);

%initialize an array of 0s for the data
spectra = zeros (1122, length(1));

%run through every pixel and preprocess spectral data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(1)
%this function aquires the spectrum data for the current pixel
    yData = S.(names{1}).data(i, :);
    [waveNumber,Raman] = ModPoly_images(xAxis',yData'); % ModPoly_images perform the preprocessing steps for the raw Raman spectrum
    [spectra(:,i)] = Raman';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Step 3: Create Reference Spectral Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% indicate wavenumber range, here: 400 - 1800cm^-1 and 2700 - 3200cm^-1
ind_ref=((waveNumbers_ref<3200)+(waveNumbers_ref>400)+(waveNumbers_ref>2700)+(waveNumbers_ref<1800))==3;

% indicate each reference, here: 7 reference spectra (column 1 to 7)
% Reference spectrum for: protein ; cyt; dna; glycogen; lipid_sat; lipid_mix; lipid_unsat
R=[refMat(ind_ref,1), refMat(ind_ref,2), refMat(ind_ref,3), refMat(ind_ref,4), refMat(ind_ref,5) refMat(ind_ref,6), refMat(ind_ref,7), refMat(ind_ref,8)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Step 4: Run the Regression Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computes number of spectra in sampleMat
count_columns=width(sampleMat);

%%% Solution using lsqnonneg for 50000 samples
[lin,col]=size(R);
coeffMatrix=zeros(col,count_columns);
contribMatrix=zeros(col,count_columns);
estimateMatrix=zeros(lin,count_columns);
residualMatrix=zeros(lin,count_columns);
count=1;

for i= 1:count_columns
    M1=sampleMat(ind_spl,i);
    [W2,resnorm,residual2,exitflag]=lsqnonneg(R,M1);
    coeffMatrix(:,count)=W2;
    contribMatrix(:,count)=W2./sum(W2);
    estimateMatrix(:,count)=S*W2;
    residualMatrix(:,count)=residual2;
    count=count+1;
end

% Create Total Lipid coeffMatrix
coeffMatrix_lipid = sum(coeffMatrix(5:7,:));

%%% Plot local concentration boxplots
figure(1);
boxplot(coeffMatrix_nan');
title('BoxPlot');
xlabel('Component','fontweight','b','fontsize',16);
ylabel('Local Concentration (mg/mL)','fontweight','b','fontsize',16);

%%% Plot protein vs cytochrome concentration correlation
figure(2);
scatter(coeffMatrix(1,:), coeffMatrix(2,:));
R= corrcoef(coeffMatrix(1,:), coeffMatrix(2,:));
R=R(1,2);
L=lsline;
L.LineWidth=3;
L.Color='r';
title('Correlation Coefficient:',num2str(R));
xlabel('Protein (mg/mL)','fontweight','b','fontsize',16);
ylabel('cytochrome (ug/mL)','fontweight','b','fontsize',16);

%%% Plot protein vs glycogen concentration correlation
figure(3);
scatter(coeffMatrix(1,:), coeffMatrix(4,:));
R= corrcoef(coeffMatrix(1,:), coeffMatrix(4,:));
R=R(1,2);
L=lsline;
L.LineWidth=3;
L.Color='r';
title('Correlation Coefficient:',num2str(R));
xlabel('Protein (mg/mL)','fontweight','b','fontsize',16);
ylabel('Glycogen (mg/mL)','fontweight','b','fontsize',16);

%%% Plot protein vs Nucleic Acid concentration correlation
figure(4);
scatter(coeffMatrix(1,:), coeffMatrix(3,:));
R= corrcoef(coeffMatrix(1,:), coeffMatrix(3,:));
R=R(1,2);
L=lsline;
L.LineWidth=3;
L.Color='r';
title('Correlation Coefficient:',num2str(R));
xlabel('Protein (mg/mL)','fontweight','b','fontsize',16);
ylabel('Nucleic Acid (mg/mL)','fontweight','b','fontsize',16);

%%% Plot protein vs Total Lipid concentration correlation
figure(5);
scatter(coeffMatrix(1,:), coeffMatrix_lipid(1,:));
R= corrcoef(coeffMatrix(1,:), coeffMatrix_lipid(1,:));
R=R(1,2);
L=lsline;
L.LineWidth=3;
L.Color='r';
title('Correlation Coefficient:',num2str(R));
xlabel('Protein (mg/mL)','fontweight','b','fontsize',16);
ylabel('Total Lipid (mg/mL)','fontweight','b','fontsize',16);

%%% Plot main outcomes
figure(6);
plot(waveNumber(ind,1),M1,'color','r','LineWidth',1.5);
hold on;
plot(waveNumber(ind,1),estimateMatrix(:,sol),'color','b','LineWidth',1.5);
set(gca,'FontWeight','bold','FontSize',32);
title('Main Outcome');
xlabel('Wavenumber','fontweight','b','fontsize',32);
ylabel('Intensity','fontweight','b','fontsize',32);
legend('Sample','Solution');

%reshape the Areas array into a matrix of the same dimensions of the image
coeff_1 = reshape(coeffMatrix(1,:), numRows, numCols);
coeff_2 = reshape(coeffMatrix(2,:), numRows, numCols);
coeff_3 = reshape(coeffMatrix(3,:), numRows, numCols);
coeff_4 = reshape(coeffMatrix(4,:), numRows, numCols);
coeff_5 = reshape(coeffMatrix(5,:), numRows, numCols);
coeff_6 = reshape(coeffMatrix(6,:), numRows, numCols);
coeff_7 = reshape(coeffMatrix(7,:), numRows, numCols);
coeff_8 = reshape(coeffMatrix(8,:), numRows, numCols);

%get the pixel dimension
pixelWidth = S.(names{1}).imageaxisscale{2,1}(2); % =0.5
pixelHeight = S.(names{1}).imageaxisscale{1,1}(2); % =1

% Save the Protein image in 'tiff' format
FileName=strcat('Protein-',sample_data,'.tiff');
coeff_1=mat2gray(coeff_1,[0 300]); % The range of intensity is limited from 0 to 300
imwrite(coeff_1,FileName);

% Save the Cytochrome image in 'tiff' format
FileName=strcat('Cytochrome-',sample_data,'.tiff');
coeff_2=mat2gray(coeff_2,[0 100]);
imwrite(coeff_2,FileName);

% Save the DNA image in 'tiff' format
FileName=strcat('DNA-',sample_data,'.tiff');
coeff_3=mat2gray(coeff_3,[0 60]);
imwrite(coeff_3,FileName);

% Save the Glycogen image in 'tiff' format
FileName=strcat('Glycogen-',sample_data,'.tiff');
coeff_4=mat2gray(coeff_4,[0 200]);
imwrite(coeff_4,FileName);

% Save the Lipid-saturated image in 'tiff' format
FileName=strcat('LipidSat-',sample_data,'.tiff');
coeff_5=mat2gray(coeff_5,[0 100]);
imwrite(coeff_5,FileName);

% Save the Lipid-mix image in 'tiff' format
FileName=strcat('LipidMix-',sample_data,'.tiff');
coeff_6=mat2gray(coeff_6,[0 100]);
imwrite(coeff_6,FileName);

% Save the Lipid-unsaturated image in 'tiff' format
FileName=strcat('LipidUnsat-',sample_data,'.tiff');
coeff_7=mat2gray(coeff_7,[0 100]);
imwrite(coeff_7,FileName);
