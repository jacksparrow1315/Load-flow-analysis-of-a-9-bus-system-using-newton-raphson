function B = line_param(L)
r = 0.00004; x = 0.000168; 
z = r + 1i*x; y = 1i*2e-5;
L = 300;
gamma = sqrt(z*y); Zc = sqrt(z/y);
Z = Zc*sinh(gamma*L);
Y = (2/Zc)*tanh(gamma*L/2);
B=Y;
end