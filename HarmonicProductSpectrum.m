% Kunal Jathal
%
% Harmonic Product Spectrum - Pitch Detection
% ===========================================

function peakPosition = HarmonicProductSpectrum(inputSignal, fs, endFactor)

% DFT
fftLength = 1024;
theFFT = abs(fft(inputSignal.*hamming(length(inputSignal)), fftLength));

tempArray = [];
finalArray = theFFT;

startFactor = 2;

% Let's downsample starting from 2 and ending at the factor entered by the
% user
for factor = startFactor:endFactor
    tempArray = downsample(theFFT, factor);

    % Before element-wise multiplication, we need to make both arrays the
    % same size, so we zero pad the temp array (which has been downsized)
    fftLength = length(theFFT);
    
    numZeroes = floor((1 - (1/factor)) * fftLength);
    tempArray = [tempArray; zeros(numZeroes, 1)];
    
    finalArray = finalArray.*tempArray;
end

% Let's get the max peak index
[maxValue, maxIndex] = max(finalArray);

% Return this index
peakPosition = (fs/fftLength) * maxIndex;

return