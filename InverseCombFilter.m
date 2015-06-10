% Kunal Jathal
%
% Inverse Comb Filtering - Pitch Detection
% ========================================

function delayN = InverseCombFilter(inputSignal, fs, lowerBound, upperBound)

% Pre-emphasis filter (which is basically a simple high pass filter)
b1 = [1 -0.99];
a1 = 1;

x1 = filter(b1, a1, inputSignal);

% Inverse comb filter and energy computation
startN = floor(fs/upperBound);
stopN = ceil(fs/lowerBound);
powerMean = zeros(1,stopN-startN);

a = 1;
r = -0.8;
b = [1 zeros(1,startN-1) r]; % set filter FIR inverse comb filter
i = 1;

for N=startN:stopN % sweep from delay that corresponds from min to max freq (in Hz)
    b = [1 zeros(1,N-1) r]; % set filter IIR inverse comb filter
    y = filter(b, a, x1);    % filter signal 

    powerMean(i) = mean(y.^2);
    
    i = i + 1;
end

[amp delayIndex] = min(powerMean);

% The delay index is what we need to return, since fs/index will give us
% the fundamental frequency

delayN = delayIndex + startN - 1;

return

