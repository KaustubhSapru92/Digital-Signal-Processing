%% Sinosoid with 45 Hz Frequency
clear all ; clf
F1 = 45;
Fs = 1000;
t = 0:1/Fs:0.999;
x1 = sin(2*pi*F1*t);
% Subsampling at a 100 Hz sampling rate
x1resampled = decimate(x1,10,1); %Downsampling by th factor 10 to achieve fs = 100Hz from 1000Hz
x1reconstructed = zeros(1,length(t)); %preallocating for speed
% Reconstruction using sinc() function
Ts = 1/100;
ts = 0:Ts:0.999;
samples = length(ts);
for i = 1:1:length(t)
    for n = 1:1:samples
        x1reconstructed(i) = sum(x1resampled .* sinc((t(i) - (1:samples)*Ts) ./ Ts)); %%% CHANGE
    end
end
%plotting
figure(2)
subplot(3,1,1)
plot(t,x1)
title('Original Signal, 45 Hz')
xlabel('Time(s)')
ylabel('Amplitude')
subplot(312)
plot(ts,x1resampled)
legend('Downsampled','original')
xlim([0.4 0.6])
hold on
plot(t,x1)
xlabel('Time(s)')
ylabel('Amplitude')
hold off
subplot(313)
plot(t,x1reconstructed)
xlim([0.4 0.6])
hold on
plot(t,x1)
legend('original','reconstructed')
xlabel('Time(s)')
ylabel('Amplitude')
hold off
%% Calculating the error
for i = 1: length(x1)
    error(i) = x1(i) - x1reconstructed(i);
end
disp('Error as percentage of variance :')
disp(mean(error.^2))