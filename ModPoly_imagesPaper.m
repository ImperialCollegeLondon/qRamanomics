%% This function is run only when 2 was pressed in driver.m, see explanation of code from ModPoly.m

function [waveNumber, raman]=ModPoly_images(x,y)
y = Min2Zero(y);
[x,y] = Water2Zero(x,y);

[x0,y0]=PBS_Norm('PBS_1.txt','PBS_2.txt','PBS_3.txt');
X=x;
Y=y-y0; % Subtract the input spectrum by PBS spectrum
index=(x>400&x<3100);
Xoriginal=X(index);
Yoriginal=Y(index);
deviation=std(Yoriginal);

% Separate the spectrums into two regions
ind =(x>400 & x<1900);
ind2 = (x>=1900 & x<3100);
X1=X(ind); % 1st region
Y1=Y(ind);
X2=X(ind2); % 2nd region
Y2=Y(ind2);
% Set the parameters for the baseline subtraction algorithm
lambda=1000000;
p=0.02;       
if deviation<0.004
    [ fluo,raman,waveNumber ] = arPLSBaseline(Xoriginal,Yoriginal);
    raman=zeros(size(waveNumber)); % Set all noise spectrums to zero
else
    [fluo1,raman1,waveNumber1 ] = whittakerbln(X1,Y1,lambda,p);
    [fluo2, raman2, waveNumber2] = arPLSBaseline(X2,Y2);
    % Combine the two regions of the spectrum
    fluo=[fluo1;fluo2];
    raman=[raman1';raman2];
    waveNumber=[waveNumber1;waveNumber2];
end
    
%%Remove the Cosmic Ray
[waveNumber,raman]=removepeak(waveNumber,raman);

end
