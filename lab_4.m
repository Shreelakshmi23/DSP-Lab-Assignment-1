%Q1.1
Fs = 120; % Sampling frequency in Hz
T=2;       %Time of Signal
t = 0:1/Fs:T-1/Fs; % Time vector
f = 15; % Frequency of the sine wave in Hz
signal = sin(2*pi*f*t);

% Step 2: Compute the DFT using the fft function
N = 120; % Number of samples to consider
DFT = fft(signal,N);

% Step 3: Calculate the frequency axis
frequencies = (0:N-1) * (Fs / N);

% Step 4: Plot the magnitude of the DFT against frequency
figure;
plot(frequencies, abs(DFT));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude of DFT vs. Frequency');
grid on;

%%
%Q1.2
Fs = 120; % Sampling frequency in Hz
t = 0:1/Fs:2-1/Fs; % Time vector
f = 15; % Frequency of the sine wave in Hz
signal = sin(2*pi*f*t);

% Step 2: Compute the DFTs of 2 different samples using the fft function
N1 = 120; % Number of samples to consider
N2 = 130;
DFT_1 = fft(signal,N1);
DFT_2 = fft(signal,N2);

% Step 3: Calculating the frequency axis of 2 samples
frequencies_1 = (0:N1-1) * (Fs / N1);
frequencies_2 = (0:N2-1) * (Fs / N2);

% Step 4: Plot the magnitude of the DFT against frequency
figure;
plot(frequencies_1, abs(DFT_1),frequencies_2, abs(DFT_2));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude of DFT vs. Frequency for Different Samples');
grid on;
legend('120 Samples','130 Samples')


%%
%Q1.3
Fs = 120; % Sampling frequency in Hz
T=2;       %Time of Signal
t = 0:1/Fs:T-1/Fs; % Time vector
f = 15; % Frequency of the sine wave in Hz
signal = sin(2*pi*f*t);

% Step 2: Compute the DFT using the fft function
N = 160; % Number of samples to consider
DFT = fft(signal,N);

% Step 3: Calculate the frequency axis
frequencies = (0:N-1) * (Fs / N);

% Step 4: Plot the magnitude of the DFT against frequency
figure;
plot(frequencies, abs(DFT));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('c');
grid on;
%%
%Q.2.1
% Define the parameters
Fs = 200;             % Sampling rate (samples per second)
T = 10;               % Duration of the signal (seconds)
N = Fs * T;           % Total number of samples
t = 0:1/Fs:(10-1/Fs);     % Time vector
signal = 0.1*sin(120* pi*t)+cos(126*pi*t);

Y = fft(signal, 215);  % Compute the DFT with 215 points
frequencies = Fs * (0:215-1) / 215;  % Frequency axis for DFT


plot(frequencies, abs(Y));
title('DFT of the Signal for 215 Samples');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%%
%Q.2.2
% Define the parameters
Fs = 200;             % Sampling rate (samples per second)
T = 10;               % Duration of the signal (seconds)
N = Fs * T;           % Total number of samples
t = 0:1/Fs:(10-1/Fs);     % Time vector
signal = 0.1*sin(120* pi*t)+cos(126*pi*t);

Y = fft(signal, 415);  % Compute the DFT with 215 points
frequencies = Fs * (0:415-1) / 415;  % Frequency axis for DFT


plot(frequencies, abs(Y));
title('DFT of the Signal for 415 Samples');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%%
%Q.2.3
% Define the parameters
Fs = 200;             % Sampling rate (samples per second)
T = 10;               % Duration of the signal (seconds)
N = Fs * T;           % Total number of samples
t = 0:1/Fs:(10-1/Fs);     % Time vector
signal = 0.1*sin(120* pi*t)+cos(126*pi*t);

Y = fft(signal, 1115);  % Compute the DFT with 215 points
frequencies = Fs * (0:1115-1) /1115;  % Frequency axis for DFT


plot(frequencies, abs(Y));
title('DFT of the Signal for 1115 Samples');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%%
%Q.2.4
% Define the parameters
Fs = 200;             % Sampling rate (samples per second)
T = 10;               % Duration of the signal (seconds)
N = Fs * T;           % Total number of samples
t = 0:1/Fs:(10-1/Fs);     % Time vector
signal = 0.1*sin(120* pi*t)+cos(126*pi*t);

Y = fft(signal, 1515);  % Compute the DFT with 215 points
frequencies = Fs * (0:1515-1) /1515;  % Frequency axis for DFT


plot(frequencies, abs(Y));
title('DFT of the Signal for 1515 Samples');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%%
%Q.2.5
% Define the parameters
Fs = 200;             % Sampling rate (samples per second)
T = 10;               % Duration of the signal (seconds)
N = Fs * T;           % Total number of samples
t = 0:1/Fs:(10-1/Fs);     % Time vector
signal = 0.1*sin(120* pi*t)+cos(126*pi*t);

Y = fft(signal, 1915);  % Compute the DFT with 215 points
frequencies = Fs * (0:1915-1) /1915;  % Frequency axis for DFT


plot(frequencies, abs(Y));
title('DFT of the Signal for 1915 Samples');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%%
% Define the parameters
Fs = 200;             % Sampling rate (samples per second)
T = 10;               % Duration of the signal (seconds)
N = Fs * T;           % Total number of samples
t = (0:N-1) / Fs;     % Time vector


signal =  0.1*sin(120* pi*t)+cos(126*pi*t);

% Define a Hamming window
window = hamming(N);

% Apply the Hamming window to the signal
windowed_signal = signal .* window';

% Plot the windowed signal
figure;
subplot(3, 1, 1);
plot(t, signal);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(t, window);
title('Hamming Window');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(t, windowed_signal);
title('Windowed Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Compute and plot the DFT of the windowed signal
Y = fft(windowed_signal, 215);  % Compute the DFT with 215 points
frequencies = Fs * (0:215-1) / 215;  % Frequency axis for DFT

figure;
plot(frequencies, abs(Y));
title('DFT of the Windowed Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%%
% Define the parameters
Fs = 200;             % Sampling rate (samples per second)
T = 10;               % Duration of the signal (seconds)
N = Fs * T;           % Total number of samples
t = (0:N-1) / Fs;     % Time vector


signal =  0.1*sin(120* pi*t)+cos(126*pi*t);

% Define a Hamming window
window = hamming(N);

% Apply the Hamming window to the signal
windowed_signal = signal .* window';
% Compute and plot the DFT of the windowed signal
Y = fft(windowed_signal, 415);  % Compute the DFT with 215 points
frequencies = Fs * (0:415-1) / 415;  % Frequency axis for DFT

figure;
plot(frequencies, abs(Y));
title('DFT of the Windowed Signal of 415 samples');
xlabel('Frequency (Hz)');
ylabel('Magnitude')
%%
% Define the parameters
Fs = 200;             % Sampling rate (samples per second)
T = 10;               % Duration of the signal (seconds)
N = Fs * T;           % Total number of samples
t = (0:N-1) / Fs;     % Time vector


signal =  0.1*sin(120* pi*t)+cos(126*pi*t);

% Define a Hamming window
window = hamming(N);

% Apply the Hamming window to the signal
windowed_signal = signal .* window';
% Compute and plot the DFT of the windowed signal
Y = fft(windowed_signal, 1115);  % Compute the DFT with 215 points
frequencies = Fs * (0:1115-1) / 1115;  % Frequency axis for DFT

figure;
plot(frequencies, abs(Y));
title('DFT of the Windowed Signal of 1115 samples');
xlabel('Frequency (Hz)');
ylabel('Magnitude')
%%
% Define the parameters
Fs = 200;             % Sampling rate (samples per second)
T = 10;               % Duration of the signal (seconds)
N = Fs * T;           % Total number of samples
t = (0:N-1) / Fs;     % Time vector


signal =  0.1*sin(120* pi*t)+cos(126*pi*t);

% Define a Hamming window
window = hamming(N);

% Apply the Hamming window to the signal
windowed_signal = signal .* window';
% Compute and plot the DFT of the windowed signal
Y = fft(windowed_signal, 1515);  % Compute the DFT with 215 points
frequencies = Fs * (0:1515-1) / 1515;  % Frequency axis for DFT

figure;
plot(frequencies, abs(Y));
title('DFT of the Windowed Signal of 1515 samples');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%%
% Define the parameters
Fs = 200;             % Sampling rate (samples per second)
T = 10;               % Duration of the signal (seconds)
N = Fs * T;           % Total number of samples
t = (0:N-1) / Fs;     % Time vector


signal =  0.1*sin(120* pi*t)+cos(126*pi*t);

% Define a Hamming window
window = hamming(N);

% Apply the Hamming window to the signal
windowed_signal = signal .* window';
% Compute and plot the DFT of the windowed signal
Y = fft(windowed_signal, 1915);  % Compute the DFT with 215 points
frequencies = Fs * (0:1915-1) / 1915;  % Frequency axis for DFT

figure;
plot(frequencies, abs(Y));
title('DFT of the Windowed Signal of 1915 samples');
xlabel('Frequency (Hz)');
ylabel('Magnitude')
%%
%q4
data = load('Exp4Data1.txt');
Fs = 1; 
T=1/Fs;

N = length(data);

%q4(a)
window = hamming(N);
data_window = data .* window';

n= 10000;
fft_new = fft(data_window, n);
freq = (0:n-1) * (Fs);
freq_normalized = freq / n;
figure;
plot(freq_normalized, abs(fft_new));
title('Magnitude of FFT with Hamming Window');
xlabel('Frequency (Fs)');
ylabel('Magnitude');

[~, max_indices] = maxk(abs(fft_new), 2);
frequencies = freq_normalized(max_indices);
disp(['Estimated Frequencies using Hamming Window: ', num2str(sort(frequencies))]);

%q4(b)
window_rect = rectwin(N);
data_window_rect = data .* window_rect';
fft_rect = fft(data_window_rect, n);
figure;
plot(freq_normalized, abs(fft_rect));
title('Magnitude of FFT with Rectangular Window');
xlabel('Frequency (Fs)');
ylabel('Magnitude');

[~, max_indices] = maxk(abs(fft_rect), 2);
frequencies1 = freq_normalized(max_indices);
disp(['Estimated Frequencies using Rectangular Window: ', num2str(sort(frequencies1))]);