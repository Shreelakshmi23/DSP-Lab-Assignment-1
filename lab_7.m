%Q1.1)
x=hamming(100);
subplot(3,1,1)
stem(x)
title('Hamming Window')
y=hanning(100);
subplot(3,1,2)
stem(y)
title('Hanning Window')
z=rectwin(100);
subplot(3,1,3)
stem(z)
title('Rectangular Window')
%%
%Q1.2)
N=1024;
L1=100;
L2=200;
L3=300;
x1=hanning(L1);
x2=hanning(L2);
x3=hanning(L3);
dft1=fftshift(fft(x1,N));
dft2=fftshift(fft(x2,N));
dft3=fftshift(fft(x3,N));
y1=(20*log10(abs(dft1)));
y2=(20*log10(abs(dft2)));
y3=(20*log10(abs(dft3)));
f=((0):(N-1))/N;
figure
plot(f,y1)
hold on 
plot(f,y2)
plot(f,y3)
title('Spectrum of Hanning Window of Different Lengths ')
ylabel('Magnitude in dB')
xlabel('Normalised Frequency')
legend('Length=100','Length=200','Length=300')
%%
figure
b = fir1(20,0.5,'low',hanning(21));
freqz(b)
%%
impz(b)
title('Impulse response of Hanning Window')
%%
figure
b1 = fir1(20,0.5,'low',rectwin(21));
freqz(b1)
%%
impz(b1)
title('Impulse response of Rectangular Window')
%%
%Q3
% Load the audio file
[y_with_dc, fs] = audioread('instru1.wav');
y= y_with_dc - mean(y_with_dc);
% Compute the spectrogram
window_size = 1024;  % Choose an appropriate window size
overlap = window_size / 2;  % Adjust overlap as needed
nfft = 1024;  % Number of FFT points
[S, F, T] = spectrogram(y_with_dc, hamming(window_size), overlap, nfft, fs);

% Plot the spectrogram
figure;
imagesc(T, F, 20*log10(abs(S)));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of Original Signal');
colorbar;

filename = 'instru1.wav';
[audio,fs] = audioread(filename);

% Perform FFT on the audio data
n = length(audio); % Length of the audio data
frequencies = (0:n-1)*(fs/n); % Frequency values corresponding to FFT bins
audio_fft = fft(audio); % Compute FFT
% Find the index of the peak magnitude in the DFT
[~, peakIndex] = max(abs(audio_fft));

% Calculate the corresponding frequency value (in Hz)
f_f = frequencies(peakIndex);


b2 = fir1(1023,[(f_f-10)*2/fs (f_f+10)*2/fs],'bandpass',hamming(1024));
% Apply the filter to extract the fundamental frequency
filtered_signal = filter(b2, 1, y);

% Write the filtered signal to a new WAV file
audiowrite('filtered_instru.wav', filtered_signal, fs);

% Listen to the filtered audio
%sound(filtered_signal, fs);

% Compute and plot the spectrogram of the filtered signal
[S_filtered, F_filtered, T_filtered] = spectrogram(filtered_signal, hamming(window_size), overlap, nfft, fs);

figure;
imagesc(T_filtered, F_filtered, 20*log10(abs(S_filtered)));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of Filtered Signal');
colorbar;