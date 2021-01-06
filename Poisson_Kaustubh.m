%% Generating random numbers for Poisson Distribution
clf
%clear all
N = 500;
R = rand(N,1);
% Measure
mu_tmp = mean(R);
std_tmp = std(R);
% Generatiing Poission distribution with a scaled mean of the numbers
k = 0:15;
lambda = 10*mu_tmp;
f = (lambda.^k) .* exp(-lambda) ./ factorial(k);
figure(1)
subplot(211)
stem(k,f,'LineWidth',1.5)
title('Poisson Distribution')
%
hold on
x = 0:15;
y = poisspdf(x,lambda);
bar(x,y,1)
xlabel('Observation')
ylabel('Probability')
% for Poisson Distribution Variance = Lambda (mean)
std_desired = sqrt(lambda);
mu_desired  = lambda;
% Normalise and denormalise to Poisson
R1 = (R - mu_tmp) / std_tmp;
R2 = (R1 * std_desired) + mu_desired;

Fano_Factor = var(R2)/mean(R2) % To check for Poisson distribution 

if Fano_Factor ==1
    disp('The newly generated sequence is Poisson distributed with Fano Factor:')
    disp(Fano_Factor)
end
subplot(212)
bar(R2)
title('Random Poisson Distributed Numbers')
xlabel('Samples')
ylabel('Value')
%% Checking the distribution of occurrance of spikes
clear all
load spike.mat
% checking for frequency of events(spike) in a given sample length(chunks
% of 4000).
k=1;
l=1;
X1=zeros(6,10);
X2=zeros(77,10);
chunk1 = 1;
chunk2 = 4000;
for j = 1:10
    for i = chunk1:chunk2
        if x(i,1) == 1
            X1(k,j) = i;
            k = k+1;
        end
        if x(i,2) ==1
            X2(l,j) = i;
            l=l+1;
        end
    end
    A(j) = k-1;
    B(j) = l-1;
    k = 1;
    l = 1;
    chunk1 = chunk2;
    chunk2 = chunk2+4000;
    if chunk2 > 40000
        break
    end
end
% count of the spikes in each chunk is k
%% plotting the probability with the k value for spike train 1
lambda = mean(A);
for i = 1:length(A)
    k(i) = A(i);
    f(i) = (lambda.^k(i)) .* exp(-lambda) ./ factorial(k(i));
end
figure(2)
clf
subplot(211)
hist(A)
% This shows for example : There were 6 spikes in the first chunk and this
% occured twice.
title('Recurrance of Spikes')
xlabel('Number of spikes in each chunk')
ylabel('Recurrance of spike pattern')
subplot(212)
stem(k,f)
title('Probabilty Mass function for channel 1')
xlabel('k')
ylabel('P(X=k)')
% Calculating the Fano factor for spike train 1
Fano_1 = var(A)/mean(A);
% Checking for poisson distributed data
if Fano_1 >0.99 && Fano_1 <1.0
    disp('The occurance of spikes in Spike train 1 is Poisson Distributed')
    disp('The Fano Factor for this data is :')
    disp(Fano_1)
end
%% plotting for the spike train 2
lambda = mean(B);
for i = 1:length(B)
    k(i) = B(i);
    f(i) = (lambda.^k(i)) .* exp(-lambda) ./ factorial(k(i));
end
figure(3)
subplot(211)
hist(B)
title('Recurrance of Spikes')
xlabel('Number of spikes in each chunk')
ylabel('Recurrance of spike pattern')
subplot(212)
stem(k,f)
title('Probabilty Mass function for channel 2')
xlabel('k')
ylabel('P(X=k)')
Fano_2 = var(B)/mean(B);
% Checking for poisson distributed data
if Fano_2 <1
    disp('The occurance of spikes in Spike train 2 is not Poisson Distributed')
    disp('The Fano Factor for this data is :')
    disp(Fano_2)
end
%%
