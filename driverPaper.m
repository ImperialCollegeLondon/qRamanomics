%Run the driver to process all Raman hyperspectral images listed in files.txt

imageList = {''};
imageListIndex = 1;
fid = fopen('files.txt', 'r');
while ~feof(fid) % Process all files listed in files.txt
    line = fgetl(fid);
    imageList{imageListIndex} = line;
    imageListIndex = imageListIndex + 1;
end

%% To output images, press 1; to generate quantitative calibration curve, press 2
prompt = 'Enter 1 for Images. Enter 2 for Calibration Curve\n';
decision = input(prompt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%% - Press 1 for images generation (5 images) - %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (decision == 1)
    for i = 1:size(imageList, 2)
        final(char(imageList(i))); % final can generate the images for different components within the reference spectrums
    end  
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%% - Press 2 for calibration curve generation - %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The script will load .mat single spectrum files exported from WITec and  
% run quantitative preprocessing algorithm on each spectrum then fit 
% regression line of intensity at specified shift vs. concentration

elseif (decision == 2)
    intensity = zeros(size(imageList, 2), 1);
   
    for i = 1:size(imageList, 2)
        intensity(i) = calibrationCurve(char(imageList(i)));
    end
   
    disp(intensity);
    %Enter calibration solution concentration values - mg/mL as X values
    X = [0,0,0,10,10,10,50,50,50,100,100,100,200,200,200,300,300,300];
    Y = intensity';
    
figure()
c = [0,0,0,2,2,2,4,4,4,6,6,6,8,8,8,10,10,10];
scatter(X,Y,500,c,'filled'); % scatter plot of the intensity
hold on;
p = polyfit(X,Y,1);
fitPoly = polyval(p,X);
R2 = norm(fitPoly - mean(Y))^2/norm(Y - mean(Y))^2;
plot(X,fitPoly);
F=strcat('y=',num2str(p(1)),'x+ ',num2str(p(1)),' ,R^2=',num2str(R2));
legend(F)
    
else
    errorMessage = 'Error: You requested a command that does not exist.';
    disp(errorMessage);
end
