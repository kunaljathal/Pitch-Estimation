% Kunal Jathal
%
% Cepstrum Analysis - - Pitch Detection
% =====================================

function sample = cepstrum(inputSignal, fs, minFrequency, maxFrequency)

% Get the DFT
fftLength = 1024;
theFFT = fft(inputSignal.*hamming(length(inputSignal)), fftLength);

% Cepstrum Computation
xCepstrum = ifft(log10(abs(theFFT)+eps));

% set min and max frequency to show
sampleStart = floor(1/maxFrequency*fs); % (max in Hz)
sampleEnd   = floor(1/minFrequency*fs);   % (min in Hz)

% Fundamental Frequency Computation
[temp fundamentalSample] = max(xCepstrum(sampleStart:sampleEnd));
sample = sampleStart + fundamentalSample;

return

