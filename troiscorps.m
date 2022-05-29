function dYdt = troiscorps (t,Y)
  global mu nu beta_mu beta_nu relat_1 relat_2

  dYdt=zeros(size(Y));
  x_1 = Y(1); x_2 = Y(5); x_p = Y(9); x_h = Y(13);
  y_1 = Y(2); y_2 = Y(6); y_p = Y(10); y_h = Y(14);
  vx_1 = Y(3); vx_2 = Y(7); vx_p = Y(11); vx_h = Y(15);
  vy_1 = Y(4); vy_2 = Y(8); vy_p = Y(12); vy_h = Y(16);
  
  r_12 = sqrt((x_2-x_1).^2+(y_2-y_1).^2);
  r_1p = sqrt((x_p-x_1)^2+(y_p-y_1)^2);
  r_2p = sqrt((x_p-x_2)^2+(y_p-y_2)^2);
  r_1h = sqrt((x_h-x_1)^2+(y_h-y_1)^2);
  r_2h = sqrt((x_h-x_2)^2+(y_h-y_2)^2);
  
  if r_1p < 6.4e6;    #condition stop si la particule tombe sur la Terre
    dYdt = zeros(size(Y));
  else
    dYdt = [vx_1,vy_1,-nu*(x_1-x_2)/r_12^3,-nu*(y_1-y_2)/r_12^3,...
    vx_2,vy_2,-mu*(x_2-x_1)/r_12^3,-mu*(y_2-y_1)/r_12^3,...
    vx_p,vy_p,-nu*(x_p-x_2)/r_2p^3 + -mu*(x_p-x_1)/r_1p^3,-nu*(y_p-y_2)/r_2p^3 + -mu*(y_p-y_1)/r_1p^3,...
    vx_h,vy_h,...
    (-1 + beta_nu)*nu*(x_h-x_2)/r_2h^3 + (-1 + beta_mu)*mu*(x_h-x_1)/r_1h^3 + relat_1*sqrt((vx_h-vx_1)^2+((vy_h-vy_1)^2))*y_h/r_1h^3 + relat_2*sqrt((vx_h-vx_2)^2+((vy_h-vy_2)^2))*y_h/r_2h^3,...
    (-1 + beta_nu)*nu*(y_h-y_2)/r_2h^3 + (-1 + beta_mu)*mu*(y_h-y_1)/r_1h^3 - relat_1*sqrt((vx_h-vx_2)^2+((vy_h-vy_2)^2))*x_h/r_1h^3 - relat_2*sqrt((vx_h-vx_2)^2+((vy_h-vy_2)^2))*x_h/r_2h^3];
  endif 
  