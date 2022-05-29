function dWdt = corps(t,W)
  global mu
dWdt=zeros(size(W));
x1 = W(1); z1=W(2);
vx1=W(3); vz1=W(4);


r = sqrt(x1^2+z1^2);

dWdt(1)=vx1;
dWdt(2)=vz1 ;
dWdt(3)= -mu*x1/r^3;
dWdt(4)= -mu*z1/r^3;

end