close all;
clear all;
clc;
t = [0:0.01:10]
k = 1;
u = 1;
f= zeros(length(t));
for x=length(t)
    f(x,:)=(k.^2*t(x).^2-2*u*t(x)+u^2)/(2*t(x))
    scatterplot(f(x))
end
% f = (k.^2*t.^2-2*u*t+u^2)/(2*t);
figure();
plot(f,t)