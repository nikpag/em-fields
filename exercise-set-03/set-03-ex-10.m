%----------------------ερώτημα γ------------------------------------------

nfig = 0;
a=1; % cylinder radius
M0=1;
m0=4*pi*10^-7;

delta_x=6*a/50;
delta_y=6*a/50;
x1=-3*a:delta_x:3*a;
y1=-3*a:delta_y:3*a;
[X,Y]=meshgrid(x1,y1);

% Calculating vector potential A
for i=1:length(x1)
    for j=1:length(y1)
        x=X(i,j);
        y=Y(i,j);
        [TH,R,Z]=cart2pol(x,y,1); % We convert cartesian to polar coordinates in order to calculate the integral
        A(i,j)=-a*integral(@(phi) f1(phi,a,R,TH),0,2*pi);
    end
end

theta=0:2*pi/50:2*pi;
x_circle=a*cos(theta);
y_circle=a*sin(theta);    % we also draw the circle corresponding to the cylindrical magnet

% Depiction of vector potential
nfig = nfig+1; figure(nfig);
surface(X,Y,A); hold on;
plot3(x_circle,y_circle,1000*ones(length(x_circle),1),'LineWidth',2,'Color','k');   % Circle that corresponds to the cylinder
shading interp;
colorbar;
axis equal;
title('Normalized vector potential A_z(x,y)/(m_0M_0/2\pi), a=1m, M_0=1Am^2');
xlabel('x')
ylabel('y')

% Field lines
nfig=nfig+1; figure(nfig);
contour_lines=2.7*[-0.9:0.1:0.9];
[CS,H]=contour(X,Y,A,contour_lines,'r','LineWidth',1);
axis equal;
hold on;
plot(x_circle,y_circle,'LineWidth',2,'Color','k');
title('Vector potential lines A_z(x,y)/(m_0M_0/2\pi), a = 1m, M_0 = 1Am^2');
xlabel('x')
ylabel('y')
clabel(CS,H,contour_lines);

%----------------------question d-----------------------------------------


% Calculating magnetic field (and magnetic induction)
for i=1:length(x1)
    for j=1:length(y1)
        x=X(i,j);
        y=Y(i,j);
        [TH,R,Z]=cart2pol(x,y,1); % converting to polar again
        % magnetic induction
        Bx(i,j)=(m0*M0/(2*pi))*a*integral(@(phi2) f2x(phi2,a,R,TH,y),0,2*pi);
        By(i,j)=-(m0*M0/(2*pi))*a*integral(@(phi2) f2y(phi2,a,R,TH,x),0,2*pi);
        % magnetic field
        if (R>a)
            Hx(i,j)=Bx(i,j)/m0;
            Hy(i,j)=By(i,j)/m0;
        else Hx(i,j)=Bx(i,j)/m0;
            Hy(i,j)=By(i,j)/m0-M0;
        end
    end
end

% depiction of magnetic induction
nfig=nfig+1; figure(nfig);
norm=sqrt(Bx.^2+By.^2); % for quiver arrow normalization
quiv=quiver(X,Y,Bx./norm,By./norm,0.5);
set(quiv,'LineWidth',0.8);
axis equal;
hold on;
plot (x_circle,y_circle,'LineWidth',2,'Color','k')
hs=streamslice(X,Y,Bx,By,1);
set(hs,'color',[0.6,0,0.1],'LineWidth',0.9);
title('Magnetic induction lines B(x,y), a=1m, M0=1Am^2')
xlabel('x');
ylabel('y');
hold off;

% depiction of magnetic field
nfig = nfig+1; figure(nfig);
norm=sqrt(Hx.^2+Hy.^2);  % for quiver arrow normalization
quiv=quiver(X,Y,Hx./norm,Hy./norm,0.5);
set(quiv,'LineWidth',0.8);
axis equal;
hold on;
plot (x_circle,y_circle,'LineWidth',2,'Color','k')
hs=streamslice(X,Y,Hx,Hy,1);
set(hs,'color',[0.6,0,0.1],'LineWidth',1);
title('Magnetic field lines H(x,y), a=1m, M0=1Am^2');
xlabel('x');
ylabel('y');
hold off;

function y=f1(phi2,a,rt,phi)    % for calculating the integrals
y=cos(phi2).*log(a./(rt.^2 +a^2 -2.*rt.*a*cos(phi-phi2)).^(1/2));
end
function y=f2x(phi2,a,rt,phi,yy) % for calculating the integrals
y=((yy-a.*sin(phi2)).*cos(phi2))./(rt.^2+a^2-2.*rt.*a.*cos(phi-phi2));
end
function y=f2y(phi2,a,rt,phi,xx)  % for calculating the integrals
y=((xx-a.*cos(phi2)).*cos(phi2))./(rt.^2+a^2-2.*rt.*a.*cos(phi-phi2));
end
