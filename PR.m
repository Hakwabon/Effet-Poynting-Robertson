function dYdt = PR(t,Y)
global c relat beta
dYdt=zeros(size(Y));
x = Y(1); z=Y(2);
vx=Y(3); vz=Y(4);
r = sqrt(x^2+z^2);

dYdt=[vx;vz;...
      (beta*x)/r^3+relat*sqrt(vx^2+vz^2)*z/r^3;(beta*z)/r^3-relat*sqrt(vx^2+vz^2)*x/r^3];

end