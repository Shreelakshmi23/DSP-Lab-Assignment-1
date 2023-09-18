%Defining the Chirp Signal
%1.1
Fs=100;  %Sampling Rate
t=0:1/Fs:10-1/Fs;  %Defining Time Period
F=(0.6).*t+4;      %Time varying Frquency

x=sin(2*pi.*F.*t);   %Chirp Signal
plot(t,x)
xlabel('Time (s)');
ylabel('Chirp Signal');
title('Signal vs Time');
%%
%1.2
%FFT of the Chirp Signal
N=1000;   %No. of Samples
DFT=fft(x,N);   %FFT of the Signal
frequencies =F;  %Frequencies
plot(frequencies,abs(DFT))
xlabel('Frequencies(Hz)');
ylabel('abs(FFT) of Signal');
title('Frequency Spectrum of Chirp Signal');
%%
%Q1.3

overlap = 10;     % Overlap size in samples
nfft = 1000;       % Number of FFT points

[S, F, T] = spectrogram(x,hamming(100), overlap,200,100);
figure;
imagesc(T, F, 20*log10(abs(S))); % Convert to dB for better visualization
axis xy; % Flip the y-axis to have lower frequencies at the bottom
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of Chirp Signal');
colorbar;
colormap('jet'); % You can choose a different colormap if desired
%%
%2.1
[y, Fs]=audioread('instru1.wav');

[S, F, T] = spectrogram(y,hamming(10000), 1000,500,Fs);
figure;
imagesc(T, F, 20*log10(abs(S))); % Convert to dB for better visualization
axis xy; % Flip the y-axis to have lower frequencies at the bottom
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of Instrument');
colorbar;
colormap('jet'); % You can choose a different colormap if desired



%%
% Load the .wav file using audioread
filename = 'instru1.wav';
[audio, sampleRate] = audioread(filename);

% Perform FFT on the audio data
n = length(audio); % Length of the audio data
frequencies = (0:n-1)*(sampleRate/n); % Frequency values corresponding to FFT bins
audio_fft = fft(audio); % Compute FFT

% Calculate the magnitude of the FFT and convert to dB scale
magnitude = 20*log10(abs(audio_fft));

% Plot the magnitude versus frequency
plot(frequencies, magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude vs. Frequency of Instrument Wave');
grid on;

% Find the index of the peak magnitude in the DFT
[~, peakIndex] = max(abs(audio_fft));

% Calculate the corresponding frequency value (in Hz)
fundamental_frequency = frequencies(peakIndex);

disp(['Fundamental Frequency: ', num2str(fundamental_frequency), ' Hz']);



%%
%2.2

[z, Fs]=audioread('Opera (2).wav');

[S, F, T] = spectrogram(z,hamming(22000), 11000,512,Fs);
figure;
imagesc(T, F, 20*log10(abs(S))); % Convert to dB for better visualization
axis xy; % Flip the y-axis to have lower frequencies at the bottom
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of Opera wav');
colorbar ;
colormap('jet'); % You can choose a different colormap if desired
%%
% Load the .wav file using audioread
filename = 'Name1.mp4';
[audio, sampleRate] = audioread(filename);
[S, F, T] = spectrogram(audio,hamming(22000), 11000,512,Fs);
figure;
imagesc(T, F, 20*log10(abs(S))); % Convert to dB for better visualization
axis xy; % Flip the y-axis to have lower frequencies at the bottom
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of Name(Shreelakshmi)');
colorbar ;
colormap('jet'); % You can choose a different colormap if desired