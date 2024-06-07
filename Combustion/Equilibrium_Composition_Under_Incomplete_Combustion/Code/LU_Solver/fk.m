function [f1,f2,f3,f4,f5] = fk(n1,n2,n3,n4,n5,Kp,lambda,a,b)
f1 = n1 + n3 - a;
f2 = (2 * n2) + (2 * n4) - b;
f3 = (2 * n1) + n2 + n3 - (2 * lambda * (a + b/4));
f4 = (2 * n5) - (2 * 3.76 * lambda * (a + b/4));
f5 = ((n3 * n2) / (n1 * n4)) - Kp;
end