% CONVOLVEDIAGRAM Plots for a diagram of convolution.

% MULTIGP

randn('seed', 1e6)
randn('seed', 1e6)

z = linspace(-10, 10, 1000)';
t = linspace(-1, 1, 100)';
kern_u = kernCreate(t, 'rbf');
kern_u.inverseWidth = 400;
Ku = kernCompute(kern_u, z);

u = gsamp(zeros(size(z)), Ku, 1)';

[T, Z] = meshgrid(t, z);
kz1 = 1/(sqrt(2*pi*(1/100)))*exp(-0.5*(T-Z).*(T-Z)/(1/30));
kz2 =  1/(sqrt(2*pi*(1/40)))*exp(-0.5*(T-Z).*(T-Z)/(1/5));

f1 = sum(kz1.*repmat(u, 1, length(t)), 1)'*((max(z)-min(z))/length(z));
f2 = sum(kz2.*repmat(u, 1, length(t)), 1)'*((max(z)-min(z))/length(z));
kern_w1 = kernCreate(t, 'rbf');
kern_w1.inverseWidth = 100;
kern_w1.variance = 0.01;
Kw1 = kernCompute(kern_w1, t);
w1 = gsamp(zeros(size(t)), Kw1, 1)';
y1= w1+f1;

kern_w2 = kernCreate(t, 'rbf');
kern_w2.inverseWidth = 100;
kern_w2.variance = 0.01;
Kw2 = kernCompute(kern_w2, t);
w2 = gsamp(zeros(size(t)), Kw2, 1)';
y2= w2+f2;

subplot(3, 3, 1)
plot(z(451:550), u(451:550));
axis off
subplot(3, 3, 2)
plot(t, kz1(500, :));
axis off
subplot(3, 3, 3)
plot(t, kz2(500, :));
axis off

subplot(3, 3, 4)
plot(t, f1);
axis off
set(gca, 'ylim', [-2 2])
subplot(3, 3, 5)
plot(t, f2);
axis off
set(gca, 'ylim', [-2 2])
subplot(3, 3, 6)
plot(t, w1);
axis off
set(gca, 'ylim', [-2 2])

subplot(3, 3, 7)
plot(t, w2);
axis off
set(gca, 'ylim', [-2 2])
subplot(3, 3, 8)
plot(t, y1);
axis off
set(gca, 'ylim', [-2 2])
subplot(3, 3, 9)
plot(t, y2);
axis off
set(gca, 'ylim', [-2 2])
