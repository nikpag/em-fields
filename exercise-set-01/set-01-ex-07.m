% Question d

% Define variables and constants
a=0.1; d=1; h=1;
e0 = 8.85e-12;
lambda_p = 4*pi*e0; % for electric potential normalization
lambda_f = 4*pi*e0/a; % for electric field normalization

[x1,y] = meshgrid(linspace(0,2,150));
[x2,z] = meshgrid(linspace(0,2,150),linspace(-2,2,200));

% calculate electric potential

Phi1=lambda_p*a/(4*pi*e0)*integral(@(phi)1./R_function(x1,y,d,h,a,phi),0,2*pi,'ArrayValued',true);
Phi2=-lambda_p*a/(4*pi*e0)*integral(@(phi)1./R_function(x1,y,d,-h,a,phi),0,2*pi,'ArrayValued',true);
Phi3=lambda_p*a/(4*pi*e0)*integral(@(phi)1./R_function(x1,y,-d,-h,a,phi),0,2*pi,'ArrayValued',true);
Phi4=-lambda_p*a/(4*pi*e0)*integral(@(phi)1./R_function(x1,y,-d,h,a,phi),0,2*pi,'ArrayValued',true);
Phi = Phi1+Phi2+Phi3+Phi4;

% equipotential lines using contour

contour_points = [0.01 0.1 0.2 0.4 0.8 1.0 1.5 2 2.5 3 3.5 4 5 6 7.5];

figure(1);
contour(x1,y,Phi,contour_points,'ShowText','On','LineWidth',1); hold on; plot([0.9 1.1],[1 1],'LineWidth',1.5,'Color',[0.8 0 0]);
xlabel('x','fontweight','bold');
ylabel('y','fontweight','bold');
title('Equipotential surfaces for normalized electric potential Phi/Vo');
figure(2);
surface(x1,y,Phi,'LineStyle','None'); hold on; plot3([0.91 1.11],[1.01 1.01],[9 9],'LineWidth',1.5,'Color',[0.8 0 0]); colorbar;
xlabel('x','fontweight','bold');
ylabel('y','fontweight','bold');
title('Color depiction of normalized potential Phi/Vo');

% Electric field lines using quiver and streamslice

[Ex,Ey] = gradient(Phi);
Ex = -Ex;
Ey = -Ey;
norm=sqrt(Ex.^2+Ey.^2);

figure(3);
hold on; plot([0.9 1.1],[1 1],'LineWidth',1.5,'Color',[0.8 0 0]); axis equal;
streamslice(x1,y,Ex./norm,Ey./norm); hold on; plot([0.9 1.1],[1 1],'LineWidth',1.5,'Color',[0.8 0 0]);
contour(x1,y,Phi,contour_points,'LineWidth',0.75,'Color','b');
xlabel('x','fontweight','bold');
ylabel('y','fontweight','bold');
title('Electric field lines');

% surface charge density lines using contour

lambda_p=1;

s1=integral(@(phi)-h./R_function_z(x2,z,d,h,a,phi).^3,0,2*pi,'ArrayValued',true);
s2=integral(@(phi)-h./R_function_z(x2,z,d,h,a,phi).^3,0,2*pi,'ArrayValued',true);
s3=integral(@(phi)h./R_function_z(x2,z,-d,h,a,phi).^3,0,2*pi,'ArrayValued',true);
s4=integral(@(phi)h./R_function_z(x2,z,-d,h,a,phi).^3,0,2*pi,'ArrayValued',true);

s=s1+s2+s3+s4;

contourpoints = [-0.1 -0.25 -0.5 -1 -2 -3 -4 -5 -7.5 -10];
figure(4); contour(x2,z,s,contourpoints,'ShowText','On','LineWidth',1.25); hold on; viscircles([h 0],a);
xlabel('x','fontweight','bold');
ylabel('z','fontweight','bold');
title('Normalized surface charge density for y = 0');

axis equal;
function R = R_function(x,y,d,h,a,phi)
    R = sqrt((x-d).^2+(y-h).^2+a^2-2*a.*(x-d).*cos(phi));
end

function R = R_function_z(x,z,d,h,a,phi)
    R = sqrt((x-d).^2+h.^2+z.^2+a^2-2*a.*((x-d).*cos(phi)+z.*sin(phi)));
end
