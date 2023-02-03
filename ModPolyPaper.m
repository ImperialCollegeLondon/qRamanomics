%% The ModPoly function is used to perform the preprocessing steps, only used when pressed one in the driver.m

function [waveNumber, raman]=ModPoly(x,y)
% The function Min2Zero subtracts the y-values of original Raman spectrum by the
% minimum value within the spectrum:

% function y=Min2Zero(y)
%  ymin = min(y);
%  y = y-ymin;
% end
y = Min2Zero(y);

% The function Water2Zero locate the water peak and normalise it to one
% function [x,y]=Water2Zero(x,y)
%     Wmin=3300;
%     Wmax=3400;
%     ind=(x>Wmin&x<Wmax);
%     A=mean(y(ind));
%     y = y/A;
% end
[x,y] = Water2Zero(x,y);

%The function PBS_Norm averages the PBS spectrums
% function [x,y]=PBS_Norm(fa,fb,fc)
%     [x1,y1]=Readfile(fa);
%     [x2,y2]=Readfile(fb);
%     [x3,y3]=Readfile(fc);
%     if (x1~=x2 | x2~=x3 | x3~=x1)
%         print('X is not same in the PBS file!')
%     end
%     y = (y1+y2+y3)/3;
%     y = Min2Zero(y);
%     [x,y] = Water2Zero(x1,y);
% end
[~,y0]=PBS_Norm('PBS_1.txt','PBS_2.txt','PBS_3.txt');

X=x;
Y=y-y0; % Subtract the orginal spectrum by the averaged PBS spectrum

index=(x>400&x<3100);
Xoriginal=X(index);
Yoriginal=Y(index);
deviation=std(Yoriginal);

% Separate the spectrum into two regions in order to perform
% different baseline fitting algorithms based on the spectrum
% characteristics
ind =(x>400 & x<1900); % The first region is from waveNumber = 400 to 1900
ind2 = (x>=1900 & x<3100); % The second region is from waveNumber = 1900 to 3100
X1=X(ind); % X1 and Y1 is the first spectrum region
Y1=Y(ind);

X2=X(ind2); % X2 and Y2 is the second spectrum region
Y2=Y(ind2);

% Set the parameters lambda and p
lambda=1000000;
p=0.02;       
if deviation<0.004
    [ fluo,raman,waveNumber ] = arPLSBaseline(Xoriginal,Yoriginal); % Use the arPLS algorithm is the std of input spectrum is less than 0.004
    raman=zeros(size(waveNumber));%if noise is interpreted as DNA, delete this line, all noise spectrums are removed by setting them to zeroes
else
    [fluo1,raman1,waveNumber1 ] = whittakerbln(X1,Y1,lambda,p); % if deviation >= 0.004, use the combination of whittaker and arPLS algorithms
    [fluo2, raman2, waveNumber2] = arPLSBaseline(X2,Y2);
    fluo=[fluo1;fluo2]; % Combine the two regions into one spectrum for the output
    raman=[raman1';raman2];
    waveNumber=[waveNumber1;waveNumber2];
end

%% Plot the input, output spectrums and the baseline generated
figure();
set(gcf, 'color', 'w');
plot(waveNumber, Yoriginal, waveNumber, fluo, waveNumber, raman)
legend({'Original Raman Signal' 'Fluorescence Background' 'Pure Raman Spectrum' }, 'FontSize', 12)
xlabel('Raman Shift (cm^{-1})', 'FontSize', 12)
ylabel('Intensity (a.u.)', 'FontSize', 12)
title('Vancouver Raman Algorithm', 'FontSize', 12)
set(gca, 'FontSize', 12)

end
