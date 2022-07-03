% define constants
a=10; b=5; d=3; L=3; D=2.5; U=1; e0=8.85e-12;
lambda=4*pi*e0;
lambda_p=-lambda*d/b;
d_p=b^2/d;
theta0=2*L/d;

% create grid
[y,z]=meshgrid(linspace(-15,15,200));
R = sqrt(y.^2+z.^2);
Rk = sqrt(y.^2+(z-D).^2);

% appropriate "masks"
Area1 = Rk<=b;
Area2 = (Rk>b).*(R<a);
Area3 = R>=a;

% calculate electric potential

factor1 = lambda*d/(4*pi*e0);
factor2 = lambda_p*d_p/(4*pi*e0);
tic;
Phi1a=factor1.*integral(@(theta)potential(y,z,d,D,theta),-theta0/2,theta0/2,'ArrayValued',true,'RelTol',1e-6,'AbsTol',1e-12);
toc;
tic;
Phi1b=factor2.*integral(@(theta)potential(y,z,d_p,D,theta),-theta0/2,theta0/2,'ArrayValued',true,'RelTol',1e-6,'AbsTol',1e-12);
toc;
Phi1=Phi1a+Phi1b+U;
Phi2=U;
Phi3=a*U./R;
Phi=Phi1.*Area1+Phi2.*Area2+Phi3.*Area3;

% design spherical conductor
rcirc1=a;
thetacirc1=linspace(0,2*pi,100);
xcirc1=rcirc1*cos(thetacirc1);
ycirc1=rcirc1*sin(thetacirc1);
zcirc1=linspace(90,90,100);

% design spherical cavity
rcirc2=b;
thetacirc2=linspace(0,2*pi,100);
xcirc2=rcirc2*cos(thetacirc2);
ycirc2=D+rcirc2*sin(thetacirc2);
zcirc2=linspace(90,90,100);

% depiction of electric potential
figure(1); surface(y,z,Phi,'LineStyle','None'); hold on; 
plot3(xcirc1,ycirc1,zcirc1,'Color','Black','LineWidth',2);
plot3(xcirc2,ycirc2,zcirc2,'Color','Black','LineWidth',2); colorbar;
caxis([0 6]); axis equal;
xlabel('y');
ylabel('z');
title('Color depiction of normalized electric potential Phi/Vo');

% equipotential surfaces using contour
figure(2); contour(y,z,Phi,[0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.10, 1.25, 2, 3, 5],'ShowText','On','LineWidth',1.25);
hold on; viscircles([0 0],a,'color','black'); viscircles([0 D],b,'color','black'); axis equal;
xlabel('y');
ylabel('z');
title('Equipotential surface for normalized electric field Phi/Vo');

% calculate electric field components
tic;E1ya=factor1.*integral(@(theta)fieldy(y,z,d,D,theta),-theta0/2,theta0/2,'ArrayValued',true,'RelTol',1e-6,'AbsTol',1e-12);toc;
tic;E1yb=factor2.*integral(@(theta)fieldy(y,z,d_p,D,theta),-theta0/2,theta0/2,'ArrayValued',true,'RelTol',1e-6,'AbsTol',1e-12);toc;
tic;E1za=factor1.*integral(@(theta)fieldz(y,z,d,D,theta),-theta0/2,theta0/2,'ArrayValued',true,'RelTol',1e-6,'AbsTol',1e-12);toc;
tic;E1zb=factor2.*integral(@(theta)fieldz(y,z,d_p,D,theta),-theta0/2,theta0/2,'ArrayValued',true,'RelTol',1e-6,'AbsTol',1e-12);toc;
E1y=E1ya+E1yb;
E1z=E1za+E1zb;
E2y=0;
E2z=0;
R3=R.^3;
E3y=a*U.*y./R3;
E3z=a*U.*z./R3;

% multiply with "masks"
Ey=E1y.*Area1+E2y.*Area2+E3y.*Area3;
Ez=E1z.*Area1+E2z.*Area2+E3z.*Area3;

% field lines using streamslice
figure(3); hold on; viscircles([0 0],a,'Color','Black'); viscircles([0 D],b,'Color','black');
norm = sqrt(Ey.^2+Ez.^2);
streamslice(y,z,Ey,Ez,2); axis equal; 
contour(y,z,Phi,[0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.10, 1.25, 2, 3, 5],'LineWidth',1.25,'Color','b');
xlabel('y','fontweight','bold');
ylabel('z','fontweight','bold');
title('Electric field lines');

% in order to calculate the field at the boundary
theta = linspace(0,2*pi,200);
y = b*sin(theta);
z = b*cos(theta)+D;

tic; E1ya=factor1.*integral(@(theta)fieldy(y,z,d,D,theta),-theta0/2,theta0/2,'ArrayValued',true,'RelTol',1e-6,'AbsTol',1e-12); toc;
tic; E1yb=factor2.*integral(@(theta)fieldy(y,z,d_p,D,theta),-theta0/2,theta0/2,'ArrayValued',true,'RelTol',1e-6,'AbsTol',1e-12); toc;
tic; E1za=factor1.*integral(@(theta)fieldz(y,z,d,D,theta),-theta0/2,theta0/2,'ArrayValued',true,'RelTol',1e-6,'AbsTol',1e-12); toc;
tic; E1zb=factor2.*integral(@(theta)fieldz(y,z,d_p,D,theta),-theta0/2,theta0/2,'ArrayValued',true,'RelTol',1e-6,'AbsTol',1e-12); toc;

E1y=E1ya+E1yb;
E1z=E1za+E1zb;

s=-e0.*(sin(theta).*E1y+cos(theta).*E1z);

% surface charge density
figure(4); plot(theta,s./lambda,'LineWidth',2);
xlabel('??? ??? ? ? ????');
ylabel('Normalized surface density sigma in m^{-1}');
title('Normalized induced surface charge density'); grid on;

function P = potential(y,z,d,D,theta)
    P = 1./sqrt(((y-d*sin(theta)).^2+(z-d*cos(theta)-D).^2));
end

function Ey = fieldy(y,z,d,D,theta)
    Ey = (y-d*sin(theta))./sqrt((y-d*sin(theta)).^2+(z-d*cos(theta)-D).^2).^3;
end

function Ez = fieldz(y,z,d,D,theta)
   Ez = (z-d*cos(theta)-D)./sqrt((y-d*sin(theta)).^2+(z-d*cos(theta)-D).^2).^3;
end