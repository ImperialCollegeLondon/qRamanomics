%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % arPLS Baseline Subtraction FOR RAMAN IMAGE DATASETS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arPLS: adaptive iteratively reweighted penalized least squares
% fluo is the y values for baseline
% raman is the y values for the output Raman spectrum with baseline
% subtracted
% waveNumber if the x values(waveNumber) for the output Raman spectrum

function [fluo, raman, waveNumber] = arPLSBaseline(x,y)
%%This function has two parameters that can be changed according to the specific Raman Spectrums: lambda and ratio
%The input x and y must be double

finalY = y; % initialise the final y axis to be the same as the input y-axis
x = double(x);% make sure the x and y are double
y = double(y);

%% Set the two parameters: lambda and ratio

% You can chose to add an additional Sovisky-Golay filter to smooth out the
% initial signal, then do baseline subtraction:
% y = sgolayfilt(y,3,91); lambda=100000; ratio=0.01;

%Set the lambda to be 10,000,000, and ratio to be 0,1

% lambda is the smoothness parameter. If lambda is too large, the fitted
% baseline would not catch and follow the curved baseline. If lambda is too
% small, the fitted baseline will follow peaks within the Raman spectrum
lambda = 10000000; 

% ratio sets the termination condition of the final spectrum, which is used
% to decide whether the baseline is good enough or not. The smaller the
% ratio, the better the baseline fits the input Raman spectrum
ratio = 0.1; 

% Estimate baseline with arPLS in Matlab
N = length(y);
D = diff(speye(N), 2); % speye(N) returns a N-by-N identity matrix, with ones on the main diagonal and zeros elsewhere
H = lambda*D'*D;
w = ones(N,1);
while true
    W = spdiags(w, 0, N, N); % creates an N-by-N sparse matrix W by taking the columns of w and placing them along the diagonals specified by 0.
    C = chol(W + H); % Cholesky decomposition of (W + H)
    z = C \ ( C'\(w.*y) ); % z is the baseline of the Raman spectrums
    d = y-z; % By subtracting the original y values from Raman spectrum by the y values of baseline, we get the ouput y values: d
    % make d negative, and get w^t with m and s
    dn = d(d<0); % dn is the negative part of d
    m = mean(dn); % mean of negative d
    s = std(dn); % standard deviation of dn
    wt = 1./ ( 1 + exp( 2* (d-(2*s-m))/s) );
    % check exit condition and backup
    if norm(w-wt)/norm(w) < ratio, break; end
    w = wt;
end

raman = finalY-z; % raman is the y values for the output Raman spectrum with baseline subtracted
fluo = z; %fluo is the y values of the baseline generated
waveNumber = x; % waveNumber is the x values for the output Raman spectrum

end
