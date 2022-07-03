%-----DISCLAIMER: 2h=2a, a, a/2 has to be changed manually----------

global a
global h
global I1
a=0.1;
I1=1;
I2=1;
h=a/4;

nfig=0;
%---------------------question d----------------------------------------
zmin = -2*h;
zmax = +2*h;

syms z

HZ=@(z)(I1.*(a.^2)/2).*((1./((a.^2+(z-h).^2).^(3/2)))+(1./((a.^2+(z+h).^2).^(3/2))));
H1=diff(HZ,z,1);
H2=diff(HZ,z,2);
H3=diff(HZ,z,3);
H4=diff(HZ,z,4);

nfig = nfig + 1; figure(nfig);
ezplot(HZ,[zmin zmax]);
set(gca,'Fontsize',12)
xlabel('z (m)','Fontsize',12)
ylabel('H (A/m)','Fontsize',12)
title('Magnetic field along z-axis','Fontsize',10)
grid on;

nfig=nfig+1; figure(nfig);
ezplot(H1,[zmin zmax]);
set(gca,'Fontsize',12)
xlabel('z (m)','Fontsize',12)
ylabel('dH/dz','Fontsize',12)
title('First derivative','Fontsize',10)
grid on

nfig=nfig+1;
figure(nfig);
ezplot(H2,[zmin zmax]);
set(gca,'Fontsize',12);
xlabel('z (m)','Fontsize',12);
ylabel('d^2H/dz^2','Fontsize',12);
title('Second derivative','Fontsize',10)
%axis equal
grid on;

nfig=nfig+1;
figure(nfig);
ezplot(H3,[zmin zmax])
set(gca,'Fontsize',12)
xlabel('z (m)','Fontsize',12)
ylabel('d^3H/dz^3','Fontsize',12)
title('Third derivative','Fontsize',10)
grid on;

nfig=nfig+1;
figure(nfig);
ezplot(H4,[zmin zmax]);
set(gca,'Fontsize',12);
xlabel('z (m)','Fontsize',12);
ylabel('d^4H/dz^4','Fontsize',12)
title('Fourth derivative','Fontsize',10);
grid on;

%-------------question e------------------

accuracy=125;

xmin = -3*a;
xmax = +3*a;
zmin1 = -(h+2*a);
zmax1 = +(h+2*a);

xx = xmin:(xmax-xmin)/accuracy:xmax;
zz = zmin1:(zmax1-zmin1)/accuracy:zmax1;

[X, Z] = meshgrid(xx,zz);


tic
for ix = 1:length(xx)
    for iz = 1:length(zz)
       x0 = X(ix,iz);
       z0 = Z(ix,iz);
       AA(ix,iz) = integral(@(f1)potentialA(f1,x0,0,z0),0,2*pi,'RelTol',1e-12,'AbsTol',1e-12);
    end
end
disp('time:')
toc

tic
for ix = 1:length(xx)
    for iz = 1:length(zz)
       x0 = X(ix,iz);
       z0 = Z(ix,iz);
       AB(ix,iz) = integral(@(f2)potentialB(f2,x0,0,z0),0,2*pi,'RelTol',1e-12,'AbsTol',1e-12);
    end
end
disp('time:')
toc

AF=AA+AB;

% plotting potential
Phi_max = max(max(AF));
Phi_min = min(min(AF));
Pm = 0.8*Phi_max;

nfig=nfig+1;
figure(nfig);
hold off
surface(X,Z,AF); shading interp
hold on

set(gca,'Fontsize',12);
xlabel('x','Fontsize',12);
ylabel('z','Fontsize',12);
title('Normalized potential A','Fontsize',10);
axis equal;
caxis([0 Pm]);
colorbar;

% equipotential lines

nfig = nfig + 1;
figure(nfig);
hold off;

%cont = [0.01, 0.025, 0.05, 0.075, 0.1];

[CS,H] = contour(X,Z,AF,'Linewidth',1,'Color','r');
%clabel(CS,H,cont);
%hold on

set(gca,'Fontsize',12);
xlabel('x','Fontsize',12);
ylabel('z','Fontsize',12);
title('Equipotential surfaces (lines) on xz-plane', 'Fontsize',10);
axis equal;
grid on;


%---------------question f---------------
% calculation and depiction of magnetic field

tic
for ix = 1:length(xx)
    for iz = 1:length(zz)
       x0 = X(ix,iz);
       z0 = Z(ix,iz);
       y0 = 0;
       HX1(ix,iz) = ((HR1(x0,y0,z0)).*(sin(theta1(x0,y0,z0))))+((HT1(x0,y0,z0)).*(cos(theta1(x0,y0,z0))));
    end
end
disp('time:')
toc

%tic
%y0=0;
%HX1 = HR1(X,y0,Z).*sin(thita1(X,y0,Z))+HT1(X,y0,Z).*cos(thita1(X,y0,Z));
%disp('time:')
%toc

tic
for ix = 1:length(xx)
    for iz = 1:length(zz)
       x0 = X(ix,iz);
       z0 = Z(ix,iz);
       y0 = 0;
       HX2(ix,iz) = ((+1).*(HR2(x0,y0,z0)).*(sin(thita2(x0,y0,z0))))+(+1).*((HT2(x0,y0,z0)).*(cos(thita2(x0,y0,z0))));
    end
end
disp('time:')
toc

tic
for ix = 1:length(xx)
    for iz = 1:length(zz)
       x0 = X(ix,iz);
       z0 = Z(ix,iz);
       y0 = 0;
       HZ1(ix,iz) = ((HR1(x0,y0,z0)).*(cos(theta1(x0,y0,z0))))-((HT1(x0,y0,z0)).*(sin(theta1(x0,y0,z0))));
    end
end
disp('time:')
toc

tic
for ix = 1:length(xx)
    for iz = 1:length(zz)
       x0 = X(ix,iz);
       z0 = Z(ix,iz);
       y0 = 0;
       HZ2(ix,iz) = ((HR2(x0,y0,z0)).*(cos(thita2(x0,y0,z0))))-((HT2(x0,y0,z0)).*(sin(thita2(x0,y0,z0))));
    end
end
disp('time:')
toc

HX_MID=HX2+HX1;
HZ_MID=HZ2+HZ1;

% we multiply with the appropriate "masks"

AreaA = (X>=0);
AreaB = (X<0);

hx_pos=HX_MID.*(AreaA);
hx_neg_mid=rot90((HX_MID),2);
hx_neg=hx_neg_mid.*(AreaB);

hz_poz=HZ_MID.*(AreaA);
hz_neg_mid=rot90((HZ_MID),2);
hz_neg=((hz_neg_mid)).*(AreaB);

HX= hx_pos+hx_neg;
HZ= hz_poz+hz_neg;

% depiction of magnetic field

nfig = nfig + 1;
figure(nfig)
hold off
norm = sqrt(((HX).^2 + (HZ).^2));
quiver(X,Z,HX./norm,HZ./norm,0.5); % for quiver arrow normalization
hold on

hs = streamslice(X,Z,HX,HZ,2);
set(hs,'Color','r','Linewidth',1.0);

set(gca,'Fontsize',12)
xlabel('x','Fontsize',12)
ylabel('z','Fontsize',12)
title('Normalized magnetic field lines','Fontsize',10)
axis equal;
grid on;

%----------------question g--------------
% calculating mutual inductance

syms x
k=x./(sqrt(h.^2+x.^2));
L=(x.*((((2-(k.^2))./k).*(KS(k))) - ((2./k).*(ES(k)))));

nfig=nfig+1;
figure(nfig);
hold off;
ezplot(L,[0 0.25]);
set(gca,'Fontsize',12);
xlabel('a(m)','Fontsize',12);
ylabel('Henry','Fontsize',12);
title('Mutual inductance as a function of the radius of the two loops','Fontsize',10);
grid on;

%---------The functions that were used

function PhiA = potentialA(fi,x0,y0,z0)

global a
I1=1;
I2=1;
global h
A=I1./4.*pi;

RA=sqrt((x0.^2 + y0.^2 + a.^2 - 2.*a.*(x0.*cos(fi)+y0.*sin(fi))+(z0-h).^2));
PhiA=A./RA;

end


% =====================================================================

function PhiB = potentialB(fi,x0,y0,z0)

global a
I1=1;
I2=1;
global h
B=I1./4.*pi;

RB=sqrt((x0.^2 + y0.^2 + a.^2 - 2.*a.*(x0.*cos(fi)+y0.*sin(fi))+(z0+h).^2));
PhiB=B./RB;

end

function resultE = E(L)

resultE = integral(@(w) (sqrt(1-(L.*((sin(w)).^2)))), 0, pi./2, 'RelTol',1e-12,'AbsTol',1e-12);

end

function resultK = K(L)

resultK = integral(@(w) (1./sqrt(1-(L.*((sin(w)).^2)))), 0, pi./2, 'RelTol',1e-12,'AbsTol',1e-12);

end

function resultES = ES(s)

resultES = (pi./2).*(1-(s./4) - (3./64).*(s.^2));

end

function resultKS = KS(s)

resultKS = (pi./2).*(1+(s./4) + (9./64).*(s.^2));

end


function resultT1 = theta1(x0,y0,z0)

global h
resultT1 = atan((sqrt((x0.^2)+(y0.^2))./(z0-h)));

end

function resultR1 = r1(x0,y0,z0)

global h
resultR1 = sqrt((x0.^2)+(y0.^2)+((z0-h).^2));

end

function resultT2 = thita2(x0,y0,z0)

global h
resultT2 = atan((sqrt((x0.^2)+(y0.^2))./(z0+h)));

end

function resultR2 = r2(x0,y0,z0)

global h
resultR2 = sqrt((x0.^2)+(y0.^2)+((z0+h).^2));

end

function resultA1 = A1(x0,y0,z0)

global a
resultA1 = sqrt((r1(x0,y0,z0)).^2 + (a.^2) - 2.*(a).*(r1(x0,y0,z0)).*(sin(theta1(x0,y0,z0))));

end

function resultA2 = A2(x0,y0,z0)

global a
resultA2 = sqrt((r2(x0,y0,z0)).^2 + (a.^2) - 2.*(a).*(r2(x0,y0,z0)).*(sin(thita2(x0,y0,z0))));

end

function resultB1 = B1(x0,y0,z0)

global a
resultB1 = sqrt((r1(x0,y0,z0)).^2 + (a.^2) + 2.*(a).*(r1(x0,y0,z0)).*(sin(theta1(x0,y0,z0))));

end

function resultB2 = B2(x0,y0,z0)

global a
resultB2 = sqrt((r2(x0,y0,z0)).^2 + (a.^2) + 2.*(a).*(r2(x0,y0,z0)).*(sin(thita2(x0,y0,z0))));

end


function resultL1 = L1(x0,y0,z0)

resultL1= 1 - (((A1(x0,y0,z0)).^2)./((B1(x0,y0,z0)).^2));

end

function resultL2 = L2(x0,y0,z0)

resultL2= 1 - (((A2(x0,y0,z0)).^2)./((B2(x0,y0,z0)).^2));

end

function resultHR1 = HR1(x0,y0,z0)

global a
global I1
C3=(I1).*(a.^2).*(cos(theta1(x0,y0,z0))).*(E(L1(x0,y0,z0)));
D3=((pi).*((A1(x0,y0,z0)).^2).*(B1(x0,y0,z0)));
resultHR1=C3./D3;

end

function resultHR2 = HR2(x0,y0,z0)

global a
global I1
C4=(I1).*(a.^2).*(cos(thita2(x0,y0,z0))).*(E(L2(x0,y0,z0)));
D4=(pi).*((A2(x0,y0,z0)).^2).*(B2(x0,y0,z0));
resultHR2= C4./D4;

end

function resultHT1 = HT1(x0,y0,z0)

global a
global I1
C1=(I1).*((r1(x0,y0,z0)).^2 + (a.^2).*(cos(2.*(theta1(x0,y0,z0))))).*(E(L1(x0,y0,z0))) - (((A1(x0,y0,z0)).^2).*(K(L1(x0,y0,z0))));
D1=(2.*pi).*(((A1(x0,y0,z0)).^2).*(B1(x0,y0,z0))).*(sin(theta1(x0,y0,z0)));
resultHT1=C1./D1;

end

function resultHT2 = HT2(x0,y0,z0)

global a
global I1
C2=(I1).*((r2(x0,y0,z0)).^2 + (a.^2).*(cos(2.*(thita2(x0,y0,z0))))).*(E(L2(x0,y0,z0))) - (((A2(x0,y0,z0)).^2).*(K(L2(x0,y0,z0))));
D2=(2.*pi).*((A2(x0,y0,z0)).^2).*(B2(x0,y0,z0)).*(sin(thita2(x0,y0,z0)));
resultHT2=C2./D2;

end
