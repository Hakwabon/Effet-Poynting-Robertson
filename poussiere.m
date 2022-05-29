clear all; close all; clc;
global c relat beta mu nu beta_mu beta_nu relat_1 relat_2 % D�finition des variables globales

%Constantes du syst�me :
c = 3e8/400;                 % Vitesse de la lumi�re et vitesse de la lumi�re normalis�e
Ls = 3.845e26;              % Luminosit� solaire
G = 6.67e-11;

%Param�tres de la poussi�re :
rho = 1000;                 % Masse volumique poussi�re
R = 10e-6;               % Rayon poussi�re  
m = 4/3*pi*R*rho*400;           % Masse poussi�re

%Calculs des forces :
mu = 2e19 ;                 % Masse soleil*G (car masse de la poussi�re est n�gligeable)
alpha = 3.845e26/(4*m*c);
beta = (alpha-mu);
rapport = alpha/mu;
%---------------------------Effet sur 2 corps---------------------------------
%Interface utilisateur :

  %Choix du temps de mod�lisation :
annee = inputdlg({'Temps de mod�lisation en ann�e'});
annee = cell2mat(annee);

  %Choix de type de mod�lisation :
menu1 = menu('Choix du mod�le','Mod�lisation de la pression de radiation seule','Ajout de l''effet relativiste');
if menu1 == 1 
  relat = 0 ;
endif
if menu1 == 2 
  relat = alpha/c;
endif
  
  
%Conditions initiales du syst�me :
x = 0; z = 1.5e11; r = sqrt(x*x+z*z);                  % Distance initiale �toile/poussi�re
v_x = sqrt(mu/r); v_z = 0  ;                                        % Vitesse initiale de la poussi�re

CI = [x,z,v_x,v_z];
t = [1:24*3600*365.25*0.01:24*3600*365.25*annee];   % Temps et pas de de mod�lisation



%Calcul des trajectoires
options=odeset('RelTol',1.e-4);
[t,y]=ode45('PR',t,CI,options);
[t,w]=ode45('corps',t,CI,options);

%Affichage
figure 1
subplot(1,2,1)
plot(y(:,1),y(:,2)),
hold on
plot(0,0,'o');
title('Trajectoire de la poussi�re avec l''effet'),
xlabel('Position en x'),
ylabel('Position en y'),
hold on


subplot(1,2,2)
plot(w(:,1),w(:,2)),
hold on
plot(0,0,'o');
title('Trajectoire de la poussi�re avec seulement la gravitation'),
xlabel('Position en x'),
ylabel('Position en y'),
axis square
hold on

%---------------------------Effet sur 3 corps---------------------------------

menu2 = menu('Choix du mod�le','Syst�me solaire','Binaire d''�toile');


if menu2 == 1;
  
  nu = 2e30*G;      %masse soleil *G
  mu = 6e24*G; 
  
  x_1 = 1.5e11; y_1 = 0;
  x_2 = 0; y_2 = 0; 

  r_12 = sqrt((x_2-x_1)^2+(y_2-y_1)^2);

  vx_1 = 0; vy_1 = sqrt(nu/r_12);
  vx_2 = 0; vy_2 = 0;

  x_p = x_1 + 3.85e8; y_p = 0;

  r_p = sqrt((x_p-x_1)^2+(y_p-y_1)^2);
  vx_p = 0; vy_p =  vy_1 + sqrt(mu/r_p);
  
  relat_1 = 0; relat_2 = alpha/c;
  beta_mu = 0; beta_nu = 0.25;

elseif menu2 == 2;
  nu = 2e30*G;      %masse soleil *G
  mu = nu;
  x_1 = 1.5e11; y_1 = 0;
  x_2 = -x_1; y_2 = 0; 

  r_12 = sqrt((x_2-x_1)^2+(y_2-y_1)^2);
  r_1 = sqrt(x_2^2+y_2^2);

  vx_1 = 0; vy_1 = r_1*sqrt(2*mu/r_12^3)
  vx_2 = 0; vy_2 = -vy_1;

  x_p = x_1*4; y_p = 0;

  r_p = sqrt(x_p*x_p+y_p*y_p);
  vx_p = 0; vy_p = sqrt((mu+nu)/r_p) ;
  
  relat_1 = alpha/c; relat_2 = alpha/c;
  beta_mu = 0.15; beta_nu = 0.15;
  
endif

t = [1:24*3600*365.25*0.01:24*3600*365.25*annee];

CI = [x_1,y_1,vx_1,vy_1,x_2,y_2,vx_2,vy_2,x_p,y_p,vx_p,vy_p,x_p,y_p,vx_p,vy_p];

%Calcul des trajectoires :
options=odeset('RelTol',5.e-6,'AbsTol',5.e-6);
[t,Y]=ode45('troiscorps',t,CI,options);
%Affichage des trajectoires :

if menu2 == 1
  corps_1 = 'Terre';
  corps_2 = 'Soleil';
  t_1 = 150;
  t_2 = 350;
else
  corps_1 = 'Etoile 1';
  corps_2 = 'Etoile 2';
  t_1 = 200;
  t_2 = 200;
endif

figure 2
subplot(1,2,1)
axis square
plot(Y(:,1),Y(:,2));
axis square
hold on
plot(Y(:,5),Y(:,6));
axis square
hold on
plot(Y(:,9),Y(:,10));
axis square
hold on
title("Sans effet relativiste")
xlabel('Coordonnees de x');
ylabel('Coordonnees de y');
axis square
hold on 
legend({corps_1,corps_2,'Poussiere'})
hold on 
subplot(1,2,2)
plot(Y(:,1),Y(:,2));
axis square
hold on
title("Poynting-Robertson")
xlabel('Coordonnees de x');
ylabel('Coordonnees de y');
axis square
plot(Y(:,5),Y(:,6));
hold on
plot(Y(:,13),Y(:,14));
axis square
hold on

a = 1.5*x_p;

figure 3


for i= 1:length(Y(:,1));
  scatter(Y(i,1),Y(i,2),t_1,'filled')
  axis ([-a a -a a])
  axis manual
  hold on 
  scatter(Y(i,5),Y(i,6),t_2,'filled')
  hold on 
  scatter(Y(i,13),Y(i,14),50,'filled')
  
  drawnow
  pause(0.001)
  hold off
endfor
