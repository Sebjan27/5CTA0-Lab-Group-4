syms A k t mu I_0
I = A*(sqrt((k)/(2*pi*t)))^((-k*(t-mu)^2/(2*t)))+I_0;
didk = diff(I,k)
didmu = diff(I,mu)