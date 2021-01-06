clear all
[x fs] = audioread('speech.wav');
time = (1:length(x))/fs;
subplot(3,1,1)
specgram(x)
xlabel('Frequency')
ylabel('Amplitude')
title('Original Windowed Signal')
%%
K = 10;
f = 1:500
T = 256;
chunk1 = 1 ;
chunk2 = 256;
for i = 1:123
    if chunk2 > 25455
        sample = x(25246:25455);
        break
    end 
    % applying harmonic model
    for j = 1:length(f)
        phase = [0:T-1]'*[1:K]*2*pi*f(j)/fs;
        S = [sin(phase) cos(phase)];
        sample = x(chunk1:chunk2);
        b = S\sample;
        h = S*b;
        error(i,j) = mean((sample-h).^2);
    end
    [m,idx]=min(error(i,:));
    phase = [0:T-1]'*[1:K]*2*pi*f(idx)/fs;
    S = [sin(phase) cos(phase)];
    b = S\sample;
    h = S*b;
    w_hanning = hanning(length(h));
    x1 = h.*w_hanning;
    X(:,i) = x1;
    pitch_f(:,i) = f(idx);  
    % creating overlapping samples of 16 ms duration
    chunk1 = chunk2-50;
    chunk2 = chunk2+205;
end
%% stiching the signal back
chunk1 =1;
chunk2 =256;
chunk_1 =1;
chunk_2 =256;
for i = 1:123 
    x_windowed(chunk1:chunk2,1) = X(:,i);
    chunk1 = chunk2;
    chunk2 = chunk2+255;
    x_signal(chunk_1:chunk_2,1) = X(:,i);
    chunk_1 = chunk_2-50;
    chunk_2 = chunk_2+205;
end

%%
subplot(3,1,2)
specgram(x_signal)
xlabel('Frequency')
ylabel('Amplitude')
title('Signal Reconstructed from Harmonic Model')
subplot(3,1,3)
plot(pitch_f)
xlim([0 123])
ylabel('Amplitude')
%% 
figure(2); 
plot(x);
hold on; 
plot(x_signal)
%legend('Original','Reconstructed')
hold off
