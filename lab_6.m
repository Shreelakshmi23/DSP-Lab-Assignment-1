%Q1- Low Pass Butterworth Filter Design
Fs = 720 ;     %Sampling Frequency in Hz
Wp =1/36 ;   %Wp is passbandedge frequency normalised 
Ws = 1/18;   %Ws is stopband edge frequency normalised
Rp = 1;      %Bandpass Ripple
Rs= 40 ;     %Attenuation in the stopband
[n,Wn] = buttord(Wp,Ws,Rp,Rs);
[b,a] = butter(n,Wn);
%H_1=tf(b,a);
H=tf(b,a,1/Fs)
freqz(b,a);
%%
% Create the pole-zero plot
figure;
pzmap(H);
%%
% Create the Bode plot
figure;
F=0:1/100:50;
bode(H,2*pi*F)
%%
%Impulse and Unit Response
figure;
impulse(H,1)
hold on 
step(H,1)
title('Impulse and Step Response')
legend('Impulse','Step')
%%
%Q2-Filtering ECG Data
data=load('ECG_Data.txt');
[b,a] = butter(n,Wn);
filtered_data = filter(b, a, data);
figure;
plot(data)
hold on
plot(filtered_data)
title('Filtered Signal');
xlabel('Samples')
ylabel('Magnitude of Signal')
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

% Design a Butterworth bandpass filter
low_cutoff = 300;  % Lower cutoff frequency in Hz
high_cutoff = 350;  % Upper cutoff frequency in Hz
order = 4;  % Filter order

[b, a] = butter(order, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');

% Apply the filter to extract the fundamental frequency
filtered_signal = filter(b, a, y);

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


%%

%Q4-Chebyshev Type-1 Filter
[n1,Wp1] = cheb1ord(Wp,Ws,Rp,Rs)
[b1,a1] = cheby1(n1,Rp,Wp1);
H1=tf(b1,a1,1/Fs)


%%
figure;

F=0:1/100:50;
bode(H,2*pi*F)

hold on
bode(H1,2*pi*F)
legend('Butterworth','Cheby1');
hold off


%%
figure;
subplot(2,1,1)
impulse(H,1)
hold on 
step(H,1)
hold off
title('Impulse and Step Response of Butterworth Filter')
legend('Impulse','Step')

subplot(2,1,2)
impulse(H1,1)
hold on
step(H1,1)
hold off
title('Impulse and Step Response of Chebyshev Type 1 Filter')
legend('Impulse','Step')
%%