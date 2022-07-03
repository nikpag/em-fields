nfig=0;
L = 0.5;
h = 0.5;
e0 =  8.8541878128e-12;
lambda0 = e0; % for normalization
e1 = e0;      % in order to have e1 = 5*e0 and e2 = e0,
e2 = 5*e0;    % we must change the values here manually

[x,z] = meshgrid(linspace(-4*L,4*L,200));      % the grid

factor1 = lambda0/(4*pi*e1);
factor1b = (e1-e2)/(e1+e2);

nom1a = L-x+sqrt((L-x).^2+(z-h).^2);
den1a = -L-x+sqrt((L+x).^2+(z-h).^2);

nom1b = L-x+sqrt((L-x).^2+(z+h).^2);
den1b = -L-x+sqrt((L+x).^2+(z+h).^2);

Phi1a = log(nom1a./den1a);
Phi1b = factor1b.*log(nom1b./den1b);
Phi1 = factor1.*(Phi1a+Phi1b);       % potential for z >= 0

factor2 = lambda0/(2*pi*(e1+e2));

nom2 = L-x+sqrt((L-x).^2+(z-h).^2);
den2 = -L-x+sqrt((L+x).^2+(z-h).^2);

Phi2 = factor2.*log(nom2./den2);     % potential for z <= 0

Area1 = z >= 0;
Area2 = z < 0;              % we create appropriate 'masks'

Phi = Phi1.*Area1+Phi2.*Area2;   % ...with which we limit the result to the desired areas

nfig=nfig+1; figure(nfig);
surface(x,z,Phi,'edgecolor','none'); % potential using surface
title('Normalized electric potential');
xlabel('x (m)');
ylabel('z (m)');
colorbar;
axis equal;
hold on;

% draw the ground for reference
line([-4*L 4*L],[0 0],[0.1 0.1],'color','black','linestyle','--','linewidth',2);

contourlines = [0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.2 0.3 0.4 0.5];

nfig=nfig+1; figure(nfig);

% draw the ground for reference
plot(linspace(-4*L,4*L,100),zeros(100),'linestyle','--','color','black','linewidth',2); hold on;

% draw linear charge distribution for reference
plot(linspace(-L,L,10),linspace(h,h,10),'color','r','linewidth',2);

% draw equipotential surfaces using contour
contour(x,z,Phi,contourlines,'showtext','on','linewidth',1.5,'linecolor','b');
axis equal; grid on;
title('Equipotential surfaces of normalized electric potential');
xlabel('x (m)');
ylabel('z (m)');

[Ex,Ez] = gradient(Phi);
Ex = -Ex; Ez = -Ez;   % calculate electric field

nfig=nfig+1; figure(nfig);
% we draw the equipotential lines so we can confirm they are perpendicular
% to the field lines
contour(x,z,Phi,contourlines,'linecolor','b');
axis equal; grid on;

% we prefer streamslice instead of quiver so the lines are easier to
% distinguish
streamslice(x,z,Ex,Ez); hold on;
title('Electric field lines');
xlabel('x (m)');
ylabel('z (m)');

% again, we draw the ground and the linear charge distribution
plot(linspace(-4*L,4*L,10),linspace(0,0,10),'linestyle','--','color','black','linewidth',2);
plot(linspace(-L,L,10),linspace(h,h,10),'color','r','linewidth',2);

Q = e0; % for normalization
N = 50;
a = 0.0025;
Dx = 2*L/N;
chi = linspace(-L,L,N+1);  % edges of conductor pieces
xm = zeros(1,N);           % middle of conductor pieces
R1 = zeros(N,N);           % distance from conductor
R2 = zeros(N,N);           % distance from conductor image
One = ones(N,1);           % matrix full of ones (denoted [1])

for i = 1:N
    xm(i) = 0.5*chi(i) + 0.5*chi(i+1);  % calculating middles
end
for i = 1:N                             % calculating distances
    for j = 1:N
        R1(i,j) = abs(xm(i)-xm(j));
        R2(i,j) = sqrt((xm(i)-xm(j))^2+4*h^2);
    end
end

nom = Dx/2+sqrt(a^2+(Dx/2)^2);
den = -Dx/2+sqrt(a^2+(Dx/2)^2);
term1 = log(nom/den);
term2 = Dx*(e1-e2)/(e1+e2)/(2*h);
Aii = 1/(4*pi*e1)*(term1+term2);        % "problematic" terms

A = zeros(N,N);
for i = 1:N                             % fill A matrix
    for j = 1:N
        if i == j
            A(i,i)=Aii;
        else
            A(i,j) = Dx/(4*pi*e0)*(1/R1(i,j)+(e1-e2)/(e1+e2)/R2(i,j));
        end
    end
end

Phi0 = 1;   % arbitrary potential, we can set this value into anything we want

Lambda0 = A\One*Phi0;

Lambda = A\One*Q/(sum(Lambda0*Dx))*Phi0;   % matrix of linear charge densities
nfig=nfig+1; figure(nfig);

% depiction of linear charge distribution, as calculated by the method of
% moments
plot(xm,Lambda./e0,'linewidth',2); grid on;
title('Linear charge distribution lambda(x)/epsilon_0 (SI)');
xlabel('x (m)');
ylabel('lambda(x)/epsilon_0');

x = linspace(-2,2,N);

factor = 1/(2*pi*(e1+e2))*Dx;
R = zeros(N,N);

for i = 1:N
    for j = 1:N
        R(i,j) = sqrt((x(i)-xm(j))^2+h^2);
    end
end

Phi = zeros(1,N);   % potential on the surface of separation

for i = 1:N
    Phi(i) = factor*sum(Lambda.'./R(i,:));
end

nfig=nfig+1; figure(nfig);
% depiction of potential on the surface of separation
plot(x,Phi,'linewidth',2); grid on;
title('Potential on surface of separation');
xlabel('x (m)');
ylabel('Normalized potential');
