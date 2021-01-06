clear all
[x fs] = audioread('speech.wav');
%load alpha.mat
Y = [zeros(length(x),33)];
time = (1:length(x))/fs;
%%  Gabor quadrature pair to bandpass the signal
clf
k = 1;
clear Y
for i = 2:0.25:10
    fc(k) = 10*(2^(i)); % center frequency in Hz
    Q = 2; % Q-factor = f/df;
    df = fc(k)/Q; % bandwidth in Hz,
    dt = 1/df;
    t = (-3*dt*fs:3*dt*fs)'/fs;
    b = 1/sqrt(pi/2)/fs/dt*exp(-t.^2/2/dt^2).*exp(sqrt(-1)*2*pi*fc(k)*t);
        
    y = filtfilt(real(b),1,x);
              
    Y(:,k) = db(abs(hilbert(y)));
    k = k+1;
end
disp('executed')
%% Plotting the instantenous amplitudes in frequency and time
clf
%set(Y,'YScale','linear')                                                                )
subplot(3,1,1)
imagesc(time,fc,Y')
xlabel('Time(s)')
ylabel('Frequency(Hz)')
axis xy
title('Instantenous Amplitude(Freq,Time) but Log Scaled')


%% interpolation of the logrithmic scale
fclinear = 0:50:8000;
Z= [zeros(length(x),length(fclinear))];

for n = 1: length(x)
     ylinear = interp1(fc,Y(n,:),fclinear);
     Z(n,:) = ylinear;
 end
disp('executed')

subplot(3,1,2)
imagesc(time,fc,Z')
title('Instantenous Amplitude(Freq,Time) Corrected')

xlabel('Time(s)')
ylabel('Frequency(Hz)')
axis xy
disp('executed')
%% creating the hanning window
subplot(3,1,3)
specgram(x,round(fs/20),fs)
title('Spectrogram with Windowing')
%%
disp('Done')

