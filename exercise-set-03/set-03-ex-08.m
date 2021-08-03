%---------------question d1----------
nfig = 0;
d = 2;
h = 1;
a = 0.1;
I = 1;
mu0 = 4*pi*10^-7;
m0 = I*pi*a^2;
factor = mu0*m0/(4*pi);
x_len = 2*d; % x-axis length in meters
y_len = 2*h; % y-axis length in meters
accuracy = 200; % number of data points per axis

z_len = x_len/2; % z-axis length in meters (-z_len < z < z_len)
[x,z] = meshgrid(linspace(0,x_len,accuracy),linspace(-z_len,z_len,accuracy));
y = 1;

% calculating distances
R1 = sqrt((x-d).^2+(y-h).^2+z.^2);
R2 = sqrt((x-d).^2+(y+h).^2+z.^2);
R3 = sqrt((x+d).^2+(y-h).^2+z.^2);
R4 = sqrt((x+d).^2+(y+h).^2+z.^2);

% calculating magnetic vector potential
A1x = z./R1.^3;
A2x = -z./R2.^3;
A3x = z./R3.^3;
A4x = -z./R4.^3;

Ax = factor*(A1x+A2x+A3x+A4x);

A1z = -(x-d)./R1.^3;
A2z = (x-d)./R2.^3;
A3z = -(x+d)./R3.^3;
A4z = (x+d)./R4.^3;

Az = factor*(A1z+A2z+A3z+A4z);

norm = sqrt(Ax.^2+Az.^2);

% nfig = nfig+1; figure(nfig);
% quiver(x,z,Ax./norm,Az./norm);
% axis equal;


% depiction of A
nfig = nfig+1; figure(nfig);
streamslice(x,z,Ax,Az);
viscircles([d 0],[a]);
xlabel('x (m)');
ylabel('z (m)');
title('Depiction of vector potential A on xz-plane for y = 1m');
axis equal;
xlim([0 x_len]);
ylim([-z_len z_len]);



%---------question ä2---------------

% same as above but on the xy-plane
[x,y] = meshgrid(linspace(0,x_len,accuracy),linspace(0,y_len,accuracy));
z = 2;

R1 = sqrt((x-d).^2+(y-h).^2+z.^2);
R2 = sqrt((x-d).^2+(y+h).^2+z.^2);
R3 = sqrt((x+d).^2+(y-h).^2+z.^2);
R4 = sqrt((x+d).^2+(y+h).^2+z.^2);

A1x = z./R1.^3;
A2x = -z./R2.^3;
A3x = z./R3.^3;
A4x = -z./R4.^3;

Ax = factor*(A1x+A2x+A3x+A4x);

Ay = zeros(accuracy,accuracy);

norm = sqrt(Ax.^2+Ay.^2);

% nfig = nfig+1; figure(nfig);
% quiver(x,y,Ax./norm,Ay./norm,0.5);
% axis equal;

nfig = nfig+1; figure(nfig);
streamslice(x,y,Ax,Ay);
line([d-a d+a],[h h],'Color','r','LineWidth',2);
axis equal;
xlim([0 x_len]);
ylim([0 y_len]);
xlabel('x (m)');
ylabel('y (m)');
title('Depiction of vector potential A on xy-plane for z = 1m');


%--------question e-----------------

[x,y] = meshgrid(linspace(0,x_len,accuracy),linspace(0,y_len,accuracy));

z = 0;

R1 = sqrt((x-d).^2+(y-h).^2+z.^2);
R2 = sqrt((x-d).^2+(y+h).^2+z.^2);
R3 = sqrt((x+d).^2+(y-h).^2+z.^2);
R4 = sqrt((x+d).^2+(y+h).^2+z.^2);

% calculating magnetic field

nom = 3*m0*(y-h).*(x-d);
den = 4*pi*R1.^5;
H1x = nom./den;

nom = -3*m0*(y+h).*(x-d);
den = 4*pi*R2.^5;
H2x = nom./den;

nom = 3*m0*(y-h).*(x+d);
den = 4*pi*R3.^5;
H3x = nom./den;

nom = -3*m0*(y+h).*(x+d);
den = 4*pi*R4.^5;
H4x = nom./den;

factor1 = 1./(4*pi*R1.^3);
factor2 = 3*m0*(y-h).^2./R1.^2 - m0;
H1y = factor1.*factor2;

factor1 = 1./(4*pi*R2.^3);
factor2 = -3*m0*(y+h).^2./R2.^2 + m0;
H2y = factor1.*factor2;

factor1 = 1./(4*pi*R3.^3);
factor2 = 3*m0*(y-h).^2./R3.^2 - m0;
H3y = factor1.*factor2;

factor1 = 1./(4*pi*R4.^3);
factor2 = -3*m0*(y+h).^2./R4.^2 + m0;
H4y = factor1.*factor2;

nom = 3*m0*(y-h).*z;
den = 4*pi*R1.^5;
H1z = nom./den;

nom = -3*m0*(y+h).*z;
den = 4*pi*R2.^5;
H2z = nom./den;

nom = 3*m0*(y-h).*z;
den = 4*pi*R3.^5;
H3z = nom./den;

nom = -3*m0*(y+h).*z;
den = 4*pi*R4.^5;
H4z = nom./den;

Hx = H1x + H2x + H3x + H4x;
Hy = H1y + H2y + H3y + H4y;
Hz = H1z + H2z + H3z + H4z;

% depiction of magnetic field H
nfig = nfig+1; figure(nfig);
streamslice(x,y,Hx,Hy);
line([d-a d+a],[h h],'Color','red','LineWidth',2);
axis equal;
xlim([0 x_len]);
ylim([0 y_len]);
xlabel('x (m)');
ylabel('y (m)');
title('Depiction of magnetic field H on xy-plane for z = 0');

%---------question f1---------------------

[y,z] = meshgrid(linspace(0,y_len,accuracy),linspace(-z_len,z_len,accuracy));

x = 0;

R1 = sqrt((x-d).^2+(y-h).^2+z.^2);
R2 = sqrt((x-d).^2+(y+h).^2+z.^2);
R3 = sqrt((x+d).^2+(y-h).^2+z.^2);
R4 = sqrt((x+d).^2+(y+h).^2+z.^2);

% we calculate the magnetic field again
nom = 3*m0*(y-h).*(x-d);
den = 4*pi*R1.^5;
H1x = nom./den;

nom = -3*m0*(y+h).*(x-d);
den = 4*pi*R2.^5;
H2x = nom./den;

nom = 3*m0*(y-h).*(x+d);
den = 4*pi*R3.^5;
H3x = nom./den;

nom = -3*m0*(y+h).*(x+d);
den = 4*pi*R4.^5;
H4x = nom./den;

factor1 = 1./(4*pi*R1.^3);
factor2 = 3*m0*(y-h).^2./R1.^2 - m0;
H1y = factor1.*factor2;

factor1 = 1./(4*pi*R2.^3);
factor2 = -3*m0*(y+h).^2./R2.^2 + m0;
H2y = factor1.*factor2;

factor1 = 1./(4*pi*R3.^3);
factor2 = 3*m0*(y-h).^2./R3.^2 - m0;
H3y = factor1.*factor2;

factor1 = 1./(4*pi*R4.^3);
factor2 = -3*m0*(y+h).^2./R4.^2 + m0;
H4y = factor1.*factor2;

nom = 3*m0*(y-h).*z;
den = 4*pi*R1.^5;
H1z = nom./den;

nom = -3*m0*(y+h).*z;
den = 4*pi*R2.^5;
H2z = nom./den;

nom = 3*m0*(y-h).*z;
den = 4*pi*R3.^5;
H3z = nom./den;

nom = -3*m0*(y+h).*z;
den = 4*pi*R4.^5;
H4z = nom./den;

Hx = H1x + H2x + H3x + H4x;
Hy = H1y + H2y + H3y + H4y;
Hz = H1z + H2z + H3z + H4z;

K1y = -Hz;
K1z = Hy;

% we depict induced surface current on the plane x = 0
nfig = nfig+1; figure(nfig);
streamslice(y,z,K1y,K1z);
line([h h],[-a a],'Color','red','LineWidth',2);
axis equal;
xlim([0 y_len]);
ylim([-z_len z_len]);
xlabel('y (m)');
ylabel('z (m)');
title('Depiction of induced surface current distribution on the plane x = 0');

%---------question f2---------------------

[x,z] = meshgrid(linspace(0,x_len,accuracy),linspace(-z_len,z_len,accuracy));
% the same as above but for y = 0
y = 0;

R1 = sqrt((x-d).^2+(y-h).^2+z.^2);
R2 = sqrt((x-d).^2+(y+h).^2+z.^2);
R3 = sqrt((x+d).^2+(y-h).^2+z.^2);
R4 = sqrt((x+d).^2+(y+h).^2+z.^2);

nom = 3*m0*(y-h).*(x-d);
den = 4*pi*R1.^5;
H1x = nom./den;

nom = -3*m0*(y+h).*(x-d);
den = 4*pi*R2.^5;
H2x = nom./den;

nom = 3*m0*(y-h).*(x+d);
den = 4*pi*R3.^5;
H3x = nom./den;

nom = -3*m0*(y+h).*(x+d);
den = 4*pi*R4.^5;
H4x = nom./den;

factor1 = 1./(4*pi*R1.^3);
factor2 = 3*m0*(y-h).^2./R1.^2 - m0;
H1y = factor1.*factor2;

factor1 = 1./(4*pi*R2.^3);
factor2 = -3*m0*(y+h).^2./R2.^2 + m0;
H2y = factor1.*factor2;

factor1 = 1./(4*pi*R3.^3);
factor2 = 3*m0*(y-h).^2./R3.^2 - m0;
H3y = factor1.*factor2;

factor1 = 1./(4*pi*R4.^3);
factor2 = -3*m0*(y+h).^2./R4.^2 + m0;
H4y = factor1.*factor2;

nom = 3*m0*(y-h).*z;
den = 4*pi*R1.^5;
H1z = nom./den;

nom = -3*m0*(y+h).*z;
den = 4*pi*R2.^5;
H2z = nom./den;

nom = 3*m0*(y-h).*z;
den = 4*pi*R3.^5;
H3z = nom./den;

nom = -3*m0*(y+h).*z;
den = 4*pi*R4.^5;
H4z = nom./den;

Hx = H1x + H2x + H3x + H4x;
Hy = H1y + H2y + H3y + H4y;
Hz = H1z + H2z + H3z + H4z;

K2x = Hz;
K2z = -Hx;

nfig=nfig+1; figure(nfig);
streamslice(x,z,K2x,K2z);
viscircles([d 0],a);
axis equal;
xlim([0 x_len]);
ylim([-z_len z_len]);
xlabel('x (m)');
ylabel('z (m)');
title('Depiction of induced surface current distribution on the plane y = 0');
