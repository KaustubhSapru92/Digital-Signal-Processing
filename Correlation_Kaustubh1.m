clear all 
load eeg.mat

x1 = x(1,:);
y1 = x(4,:);
y2 = x(5,:);
y3 = x(6,:);

subplot(221)
plot([ x(1,:); 2000+x(4,:); 4000+x(5,:); 6000+x(6,:)]')
N = length(x1)
xlabel('Samples')
ylabel('x,y1,y2,y3')
title('y,x')
%%
Rxx = x1*x1'/length(x1);
Rxy = x1*y1'/length(x);
A = Rxy\Rxx;
n1 = y1-A*x1;
%n2
Rxy = 1/N * y2*x1';
A = Rxy*Rxx^-1;
n2 = y2 - A*x1;
%n3
Rxy3 = 1/N * y3*x1';
Axy3 = Rxy3*Rxx^-1;
n3 = y3 - A*x1;
%%
subplot(222)
plot([2000+n1; 4000+n2; 6000+n3]')
ylim([0 7000])
xlabel('Samples')
title('n = y-Ax')
legend('n1','n2','n3')
%%
subplot(223)
plot(x1,[y1; y2; y3]','.')
xlabel('x')
ylabel('y')
subplot(224)
plot(x1,[n1; n2; n3]','.')
xlabel('x')
ylabel('n')
%% Multivariate Normal distribution
x1 = -2:0.1:2;
x2 = -2:0.1:2;
n = length(x1)

Rxx = [1 0.5;0.5 0.4]
P = zeros(41,41);
for i = 1:length(x1)
    for k = 1:length(x2)
        M = [x2(k) x1(i)];
        p(i,k) = ((sqrt(((pi*2)^n)*det(Rxx)))^-1)*(exp((-1*(M)*inv(Rxx)*M')/2));
    end

end

P = rescale(p);
figure(2)
mesh(x1,x2,P/2.5)
zlim([0 0.5])
zlabel('Probability Density')
[U,D]= eig(Rxx);
z = randn(2,1000);
W = U*sqrt(D);
x = W*z;

% Contour plot
figure(3);
contour(x1,x2,p)
title('Contour Plot')
hold on
plot(x(1,:),x(2,:),'.')
xlim([-2 2])
hold off
%
%% Slide 18
figure(4);
plot(x(1,:)',x(2,:)','.')
xlim([-4 4])
ylim([-2 2.5])
hold on
d = sqrt(diag(D));
h1 = quiver(0, 0, U(1, 2), -U(2, 2), 7.5*d(1), 'LineWidth', 2);
h2 = quiver(0, 0, U(1, 1), -U(2, 2), 3*d(2), 'LineWidth', 2);
hold off
title('Principle Components(u1,u2)')