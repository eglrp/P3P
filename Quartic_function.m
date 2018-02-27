function [x1,x2,x3,x4] = Quartic_function(a,b,c,d,e)

p = (8*a*c - 3*b^2)/8*a^2;
q = (b^3 - 4*a*b*c + 8*a^2*d)/8*a^3;

delta0 = c^2 - 3*b*d + 12*a*e;
delta1 = 2*c^3 - 9*b*c*d +27*b^2*e + 27*a*d^2 - 72*a*c*e;

Q = ((delta1 + (delta1^2 - 4*delta0^3)^0.5)/2)^(1/3);
S = 0.5*(-2/3*p + 1/(3*a)*(Q + delta0/Q))^0.5;

x1 = -b/(4*a) - S + 0.5*(-4*S^2 - 2*p + q/S)^0.5;
x2 = -b/(4*a) - S - 0.5*(-4*S^2 - 2*p + q/S)^0.5;
x3 = -b/(4*a) + S + 0.5*(-4*S^2 - 2*p - q/S)^0.5;
x4 = -b/(4*a) + S - 0.5*(-4*S^2 - 2*p - q/S)^0.5;