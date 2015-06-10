% Kunal Jathal
%
% Fundamental Frequency Computation
% =================================

%% USAGE

% FundamentalFreqComputation(input, method)

% input - input signal (egs: 'bass_clarinet_fhorn.wav')

% method - Method by which you want to calculate the fundamental frequency

% 'zcr'                 - Zero Crossing
% 'auto'                - Auto Correlation
% 'icf'                 - Inverse Comb Filtering
% 'hps'                 - Harmonic Product Spectrum
% 'cepstrum'            - Cepstrum Analysis
% 'chroma'              - Chroma

%%
function FundamentalFreqComputation(inputWaveForm, method)

% Read the original signal
[inputSignal, fs] = wavread(inputWaveForm);

% Remove DC
inputSignal = inputSignal - mean(inputSignal);

% Get input signal length
signalLength = length(inputSignal);

% Get variables ready
frameLength = 256;
frameIndex = 1;
rmsIndex = 1;
xSignalStartIndex = 1;
xSignalEndIndex = 1;

% Get number of frames
numFrames = ceil(signalLength / frameLength);

xRMS = zeros(numFrames, 1);
xSignalStart = [];
xSignalEnd = [];

%% SILENCE DETECTION

% We want to break the signal into frames. Then, we compute the RMS of each
% frame. Store all the RMS values in an array. Finally, find the indices
% where the RMS is less than 5% of the max RMS. These are the indices that
% are the boundaries of the frames, so we can split the final signal into
% the required frames.

% Get the RMS array
while(1)
    if frameIndex + frameLength > signalLength
        break; % reached end
    else
        xRMS(rmsIndex) = sqrt(mean(sum((inputSignal(frameIndex:(frameIndex + frameLength - 1))).^2)));
        rmsIndex = rmsIndex + 1;
        frameIndex = frameIndex + frameLength;
    end        
end

% Now, we have an array of all the RMS values. Let's get the indices for
% where the signal portions of the entire wave start and stop. The
% threshold we use to gauge this is 5% of the max RMS value.
maxRMS = max(xRMS);
threshold = 0.05 * maxRMS;
signalFlag = 0;

% initialize RMS processing
if xRMS(1) < threshold
    signalFlag = -1;
else
    signalFlag = 1;
    xSignalStart(1) = 1;
    xSignalStartIndex = xSignalStartIndex + 1;
end

% step through RMS array
for rmsIndex = 2:length(xRMS)
    if (signalFlag < 0)
        if xRMS(rmsIndex) >= threshold
            % threshold change occured from silence to no silence
            signalFlag = 1;    % flip flag
            xSignalStart(xSignalStartIndex) = rmsIndex;    % save sample number
            xSignalStartIndex = xSignalStartIndex + 1;  % increment array
        end
    else
        if xRMS(rmsIndex) < threshold
            % threshold change occured from no silence to silence
            signalFlag = -1;   % flip flag
            xSignalEnd(xSignalEndIndex) = rmsIndex - 1;    % save sample number                    
            xSignalEndIndex = xSignalEndIndex + 1;  % increment array
        end        
    end
end

if xSignalStartIndex > xSignalEndIndex
    % The signal ends in non-silence
    xSignalEnd(xSignalEndIndex) = length(xRMS);
end

% Now we have rid the signal off any silence. Pass each of these signals to
% the corresponding fundamental frequency computation algorithms.

%% FUNDAMENTAL FREQUENCY DETECTION

% Now we send each part of the signal to be analyzed for fundamental
% frequency estimation.

% For each method, there are a set of variables that the user needs to
% enter, that are CRUCIAL in determining the accuracy of the fundamental
% frequency detection. These variables are initialized here...

windowLength = 0; % Zero Crossing
threshold = 0; % AutoCorrelation
lowerBound = 0; % Inverse Comb Filtering
upperBound = 0; % Inverse Comb Filtering
endFactor = 0; % Harmonic Product Spectrum
minFrequency = 0; % Cepstrum Analysis
maxFrequency = 0; % Cepstrum Analysis

% Ask user to enter values for the required variables. In each case, I
% provide example values that I recommend the user enter.
if (strcmp(method, 'zcr'))    
    windowLength = input(['For the Zero Crossing Method, the window length is crucial in \n' ...
    'determining the fundamental frequency. Please enter the window length you would like \n' ...
    'to use in SAMPLES (for egs, to use 50 samples, enter 50): ']);
elseif (strcmp(method, 'auto'))
    threshold = input(['For the AutoCorrelation method, the threshold value for peaks to be used \n'...
        'in computing the peak period is vital. Please enter a percentage of the MAX peak value \n'...
        'you would like to use as the threshold. For egs, to only use peaks that are at least 50% \n'...
        'of the max peak, enter 0.5. : ']);
elseif (strcmp(method, 'icf'))
    lowerBound = input(['For the Inverse Comb Filtering method, the lower and upper bounds of the \n'...
        'frequency to sweep are vital in determining the fundamental frequency. Please enter the \n'...
        'lower bound in Hz (egs: for 40 Hz, enter 40): ']);
    upperBound = input('Enter the upper bound in Hz (egs: for 1500, enter 1500): ');
elseif (strcmp(method, 'hps'))
    endFactor = input(['For the Harmonic Product Spectrum method, the ending factor by which we \n'...
        'downsample the signal is vital in determining the fundamental frequency. Please enter \n'...
        'enter the ending factor (egs: to stop computation after downsampling by a factor of 4, enter 4): ']);
elseif (strcmp(method, 'cepstrum'))    
    minFrequency = input(['For the Cepstrum method, the lower and upper bounds of the frequencies to show \n'...
        'are vital in determining the fundamental frequency. Please enter the lower bound of the \n'...
        'frequency to show in Hz (egs: for 30 Hz, enter 30): ']);
    maxFrequency = input('Enter the upper bound in Hz (egs: for 1000, enter 1000): ');
end


for counter=1:length(xSignalStart)
    
    startIndex = ((xSignalStart(counter) - 1) *  frameLength) + 1;
    endIndex = (xSignalEnd(counter)) * frameLength;
    sig = inputSignal(startIndex:endIndex);
    
    switch method
        case 'auto'
            peakPeriod = AutoCorrelation(sig, threshold);
            % Fund Freq is sample rate over the peak period
            fundFreq = fs/peakPeriod;
            display(strcat('The pitch of sound number #', num2str(counter), ' using autocorrelation is (in Hz) : ', num2str(fundFreq)));
        case 'zcr'
            [zcr, zc] = ZeroCrossing(sig, windowLength);            
            % Fund Freq is sample rate over the ZCR
            fundFreq = fs/(zcr*2);
            display(strcat('The pitch of sound number #', num2str(counter), ' using ZCR is (in Hz) : ', num2str(fundFreq)));
        case 'icf'
            delayIndex = InverseCombFilter(sig, fs, lowerBound, upperBound);
            % Fund Freq is sample rate over the delay index
            fundFreq = fs/delayIndex;
            display(strcat('The pitch of sound number #', num2str(counter), ' using Inverse Comb Filtering is (in Hz) : ', num2str(fundFreq)));
        case 'cepstrum'
            sample = cepstrum(sig, fs, minFrequency, maxFrequency);
            % Fund Freq is sample frequency divided by the sample index
            fundFreq = fs/sample;
            display(strcat('The pitch of sound number #', num2str(counter), ' using Cepstrum Analysis is (in Hz) : ', num2str(fundFreq)));
        case 'hps'
            fundFreq = HarmonicProductSpectrum(sig, fs, endFactor);
            % Fund Freq is the peak position
            display(strcat('The pitch of sound number #', num2str(counter), ' using Harmonic Product Spectrum is (in Hz) : ', num2str(fundFreq)));
        case 'chroma'
            fundFreq = chroma(sig, fs);
            % Fund Freq is the chroma pitch
            display(['The pitch of sound number #' num2str(counter) ' using Chroma is (in Hz) : ' num2str(fundFreq)]);
        otherwise
            error('Invalid method entered...');
    end
end

end


