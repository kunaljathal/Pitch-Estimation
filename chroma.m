% Kunal Jathal
%
% Chroma - - Pitch Detection
% ==========================

function chromaPitch = chroma(inputSignal, fs)

% Let's first make our chroma 'constant Q' pitch classes. We will go from
% C0 to C7. C0 is 16.35 Hz. C7 is 84 semitones above C0.

constantC0 = 16.35;
chromaPitches = [constantC0];

for semitone=1:84
    newSemitone = constantC0 * (2^(semitone/12));
    chromaPitches = [chromaPitches newSemitone];
end

% This is the array where we will store the cumulative chroma energies
energyBuckets = zeros(size(chromaPitches));

% Let's get the FFT of our input signal now
fftLength = 2^nextpow2(length(inputSignal));
theFFT = abs(fft(inputSignal.*hamming(length(inputSignal)), fftLength));

% We want to get the index range of FFT indices that we want to group
% together that correspond to the chroma pitch classes. Let's go through
% each of the chroma pitches to determine this.

indexMapping = chromaPitches * (fftLength/fs);

% Now we want to use the indexes in the indexMapping array as boundaries
% for the chroma pitch classes.

for i=1:length(indexMapping) - 1
    startIndex = round(indexMapping(i));
    endIndex = round(indexMapping(i+1)) - 1;
    energyBuckets(i) = mean(theFFT(startIndex:endIndex).^2);
end

% Now, all we need to do is go through the energy buckets and see which
% bucket has the most energy:

[maxEnergy, maxIndex] = max(energyBuckets);

% To get the fundamental frequency, get the index. The index is the
% semitone. The fundamental frequency will be the constant Q of it.

chromaPitch = constantC0 * (2^(maxIndex/12));

return
