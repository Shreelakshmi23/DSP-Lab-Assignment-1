a = 1+mod(147,3);
%%
%Problem 1
%1)
t = [0:0.0001:1/a];
figure;
y1 = a*cos(2*pi*5*a*t);
y2 = (a/2)*cos(2*pi*6*a*t);
y3 = (a/4)*cos(2*pi*10*a*t);

plot(t, y1, t, y2, t, y3)
xlabel('Time (t) in seconds')
ylabel('y (t)')
title('1.Waveforms')
legend('(a)3cos(2*pi*5*3*t)','(b)(3/2)*cos(2*pi*6*3*t)','(c)(3/4)*cos(2*pi*10*3*t)');
%%
%2)
figure;
y = y1 + y2 + y3;
plot(t,y);
xlabel('Time (t) in seconds')
ylabel('y (t)')
title('2.')
%%
%3)
figure;
x = (3 * 10)-1 %10 considering the max freq
subplot(3,1,1);
%a 42Hz
Fs2 = 14 * a;
n2 = [0:1/Fs2:1/a];
y1 = a *cos(2*pi*5*a*n2);
y2 = a/2 *cos(2*pi*6*a*n2);
y3 = a/4 *cos(2*pi*10*a*n2);
y11 = y1 + y2 + y3;
stem(y11);
subplot(3,1,2);
%b Nyquist rate => 2 * max freq => 2 * 10  *a =>20 * a 60 Hz.
Fs2 = 20 * a;
n2 = [0:1/Fs2:1/a];
y1 = a *cos(2*pi*5*a*n2);
y2 = a/2 *cos(2*pi*6*a*n2);
y3 = a/4 *cos(2*pi*10*a*n2);
y22 = y1 + y2 + y3;
stem(y22);
%plot 0 - len-1 , in samples
subplot(3,1,3);
% 6a aliased with 3a => 18Hz aliased with 9Hz => if Sampling freq is 18 + 9 =
% 27 Hz, 27/2 = 13.5 Hz So upon folding at 13.5, 18 overlaps 9.
Fs3 = 27;
n2 = [0:1/Fs3:1/a];
y1 = a *cos(2*pi*5*a*n2);
y2 = a/2 *cos(2*pi*6*a*n2);
y3 = a/4 *cos(2*pi*10*a*n2);
y33 = y1 + y2 + y3;
stem(y33);
%%
%4)
figure;
subplot(3,1,1);
plot(y11);
subplot(3,1,2);
plot(y22);
subplot(3,1,3);
plot(y33);
%%
%Problem 2)
s1 = sin(2*pi*240*t);
s2 = sin(2*pi*270*t);
s3 = sin(2*pi*300*t);
s4 = sin(2*pi*320*t);
s5 = sin(2*pi*360*t);
s6 = sin(2*pi*400*t);
s7 = sin(2*pi*450*t);
s8 = sin(2*pi*480*t);
s = [s1,s2,s3,s4,s5,s6,s7,s8];
audiowrite("Doremi_48000.wav",s,48000)

% sound(s,48000)
% sound(s,24000)
% sound(s,16000)
% sound(s,8000)
% sound(s,4000)

%increasing sampling rate increases both speed of the audio as well as
%loudness and pitch
%%
%Problem 3
[data,Fs] = audioread("Track003.wav")
sound(data,Fs)
length(data)
data2 = data(1:3:end, 1);
length(data2)
%%
sound(data2,Fs/3)
%%
data3 = data(1:5:end, 1);
length(data3)
sound(data2,Fs/5)
%higher frequencies seem absent as we increase the frequencies