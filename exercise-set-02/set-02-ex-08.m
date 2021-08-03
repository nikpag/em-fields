nfig=0;
N = 3;    % every time we change N (first 3, then 5-75)
a = 0.005;
L = 0.99;
h = 1.00;
D = 1.5;
s0 = 1/160;
I = 500;
Dx = D/N;
Dz = L/N;

xi = zeros(1, 3*N+3);       % coordinates for the edges of the pieces
zi = zeros(1, 3*N+3);

for i = 1:N+1               % left rod bottom-to-top
    xi(i)=-D/2+a;
    zi(i)=-h-a-L+(i-1)*Dz;
end

for i = N+2:2*N+2           % horizontal rod left-to-right
    xi(i)=-D/2+(i-N-2)*Dx;
    zi(i)=-h;
end

for i = 2*N+3:3*N+3         % right rod top-to-bottom
    xi(i)=D/2-a;
    zi(i)=-h-a-(i-2*N-3)*Dz;
end

xm = zeros(1, 3*N);         % coordinates for the center of the pieces
zm = zeros(1, 3*N);

for i = 1:N                 % left rod
    xm(i)=0.5*(xi(i)+xi(i+1));
    zm(i)=0.5*(zi(i)+zi(i+1));
end

for i = N+2:2*N+1           % horizontal rod
    xm(i-1)=0.5*(xi(i)+xi(i+1));
    zm(i-1)=0.5*(zi(i)+zi(i+1));
end

for i = 2*N+3:3*N+2         % right rod
    xm(i-2)=0.5*(xi(i)+xi(i+1));
    zm(i-2)=0.5*(zi(i)+zi(i+1));
end

factor = 1/(4*pi*s0);
Duj = zeros(1,3*N);         % dx for horizontal rod and dz for perpendicular rods

for j = 1:N                 % left rod
    Duj(j) = Dz;
end

for j = N+1:2*N             % horizontal rod
    Duj(j) = Dx;
end

for j = 2*N+1:3*N           % right rod
    Duj(j) = Dz;
end

Rij = zeros(3*N,3*N);       % distances from 'grounder'
Rijprime = zeros(3*N,3*N);  % distances from 'grounder' image

for i = 1:3*N               % calculating distances
    for j = 1:3*N
        Rij(i,j) = sqrt((xm(i)-xm(j))^2+(zm(i)-zm(j))^2);
        Rijprime(i,j) = sqrt((xm(i)-xm(j))^2+(zm(i)+zm(j))^2);
    end
end

VDFij = zeros(3*N,3*N);     % mutual resistance matrix

for i = 1:3*N               % fill matrix, considering 'problematic' terms
    for j = 1:3*N
        if i == j
            nom = Duj(j)/2+sqrt(a^2+(Duj(j)/2)^2);
            den = -Duj(j)/2+sqrt(a^2+(Duj(j)/2)^2);
            term1=log(nom/den);
            term2=Duj(j)/Rijprime(i,i);
            VDFij(i,j) = factor*(term1+term2);
        else
            VDFij(i,j) = factor*Duj(j)*(1/Rij(i,j)+1/Rijprime(i,j));
        end
    end
end

Phi0 = 1;   % arbitrary potential, can be whatever we want
One = ones(3*N,1);  % matrix of ones (denoted [1])

I0 = VDFij\One.*Phi0;      % current matrix corresponding to Phi0

Imatrix = VDFij\One.*I./(Duj*I0)*Phi0;      % matrix of current distribution

% here we depict the charge distribution in every rectilinear conductor
% separately

nfig=nfig+1; figure(nfig);
plot(zm(1:N),Imatrix(1:N),'linewidth',2); % left
grid on;
title('Current distribution of left conductor for N = ' + N);
xlabel('z (m)');
ylabel('Linear current distribution (A/m)');

nfig=nfig+1; figure(nfig);
plot(xm(N+1:2*N),Imatrix(N+1:2*N),'linewidth',2);   % horizontal
grid on;
title('Current distribution of horizontal conductor for N = ' + N);
xlabel('x (m)');
ylabel('Linear current distribution (A/m)');

nfig=nfig+1; figure(nfig);
plot(zm(2*N+1:3*N),Imatrix(2*N+1:3*N),'linewidth',2);     % right
grid on;
title('Current distribution of right conductor for N = ' + N);
xlabel('z (m)');
ylabel('Linear current distribution (A/m)');

Phi = Phi0*I/(Duj*I0);

Rg = Phi0./(Duj*I0);      % ground resistance

size = 100; % length of x vector

x = linspace(-5,5,size);     % x-coordinates of ground

factor = 1/(2*pi*s0);

Rj = zeros(size,3*N);

for i = 1:size
    for j = 1:3*N
        Rj(i,j) = sqrt((x(i)-xm(j))^2+zm(j)^2);
    end
end

Phig = zeros(1,size);      % ground potential

for i = 1:size
    Phig(i) = factor*sum(Imatrix.'.*Duj./Rj(i,:));
end

Rj0 = sqrt(xm.^2+zm.^2);

PhiO = factor*sum(Imatrix.'.*Duj./Rj0);    % potential at O

nfig=nfig+1; figure(nfig);
plot(x,Phig,'linewidth',2);  % depiction of ground potential
grid on;
title('Ground potential for N = ' + N);
xlabel('x (m)');
ylabel('Potential (V)');
