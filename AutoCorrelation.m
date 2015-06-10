% Kunal Jathal
%
% Auto Correlation - Pitch Detection
% ==================================

function peakPeriod = AutoCorrelation(inputSignal, thold)

% Get the autocorrelation vector of the signal
autoCorr = xcorr(inputSignal);

% Now we need to get the period between peaks to get the fundamental freq.

% initialize sign
if ((autoCorr(2) - autoCorr(1))/2) < 0
    theSign = -1;
else
    theSign = 1;
end

% step through waveform
idx = 1;
allPeaks = [];
for i=3:length(autoCorr)
    if (theSign < 0)
        if ((autoCorr(i) - autoCorr(i-1))/2) >= 0
            % slope change occured from - to +
            theSign = 1;    % flip sign
        end
    else
        if ((autoCorr(i) - autoCorr(i-1))/2) < 0
            % zero crossing occured from + to -
            theSign = -1;   % flip sign
            allPeaks(idx) = i - 1;    % save sample number                    
            idx = idx + 1;  % increment all Peaks array
        end        
    end
end


% We only care about peaks whose value are above the threshold...
% For the threshold, we will use the % of the max entered by the user
threshold = thold * max(autoCorr);

peaks = [];
peakCounter = 1;

% Build a new array that contains only the peaks of interest
for i=1:length(allPeaks)
    if (autoCorr(allPeaks(i)) > threshold)
        peaks(peakCounter) = allPeaks(i);
        peakCounter = peakCounter + 1;
    end
end

% Now, peaks contains the list of peaks that are above the threshold. All
% we need to do is compute the mean distance between them, which is P

peakPeriod = mean(diff(peaks)); 

return
