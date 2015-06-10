% Kunal Jathal
%
% Zero Crossing Rate - Pitch Detection
% ====================================

function [zcr, zc] = ZeroCrossing(x, windowLength)

% Get variables ready
xMean = zeros(1,length(x));
lenX1 = length(x);
windowIndex = 1;
meanIndex = 1;
zc = [];

% Moving Average Filter to smooth out the signal
while(1)
    if windowIndex + windowLength > lenX1
        break; % reached end
    else
        xMean(meanIndex) = mean(x(windowIndex:(windowIndex+windowLength-1)));
        meanIndex = meanIndex + 1;
        windowIndex = windowIndex + windowLength;
    end
end

% interpolate back to original signal size
xMean = interp(xMean, windowLength);

% initialize sign
if xMean(1) < 0
    theSign = -1;
else
    theSign = 1;
end

% step through waveform
idx = 1;
for i=2:length(xMean)
    if (theSign < 0)
        if xMean(i) >= 0
            % zero crossing occured from - to +
            theSign = 1;    % flip sign
            zc(idx) = i;    % save sample number
            idx = idx + 1;  % increment zc array
        end
    else
        if xMean(i) < 0
            % zero crossing occured from + to -
            theSign = -1;   % flip sign
            zc(idx) = i;    % save sample number                    
            idx = idx + 1;  % increment zc array
        end        
    end
end

% ZCR is basically the number of samples between successive zero crossings.
zcr = mean(diff(zc)); 

return
