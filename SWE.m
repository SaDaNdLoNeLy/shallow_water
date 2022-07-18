clf;
clear all;
%% define the grid size
n = 100;
dt = 0.02;
dx = 1;
dy = 1;
g = 9.8;

 
H = ones(n+2,n+2); % displacement matrix (+2 because of double step)
U = zeros(n+2,n+2); % x velocity
V = zeros(n+2,n+2); % y velocity

%% draw the mesh
grid = surf(H);
axis([1 n 1 n 0 3]);
hold all;

%% create initial displacement
[x,y] = meshgrid( linspace(-3,3, 20) );
R = sqrt(x.^2 + y.^2);
Z = 1*(sin(R)./R);
Z = max(Z,0);

% second displacement
% Z1 = 1*(sin(R)./R);
% Z1 = max(Z,0);

% add displacement to the height matrix

i = 60:79;
j = 60:79;
H(i,j) = H(i,j) + Z;

% add second displacement
%i1 = 30:49;
%j1 = 40:59;
%H(i1, j1) = H(i1, j1) + Z1;

%% empty matrix for half-step calculations
Hx = zeros(n+1,n+1); 
Hy = zeros(n+1,n+1);
Ux = zeros(n+1,n+1); 
Uy = zeros(n+1,n+1);
Vx = zeros(n+1,n+1);
Vy = zeros(n+1,n+1);
 
%% loop through time
while 1==1
 
 % redraw the mesh
 set(grid, 'zdata', H);
 set(gcf, 'Position',  [500, 200, 700, 700]);
 drawnow
 % blending the edges keeps the function stable
 H(:,1) = H(:,2); 
 H(:,n+2) = H(:,n+1); 
 H(1,:) = H(2,:); 
 H(n+2,:) = H(n+1,:); 
 
 % reverse direction at the x edges
 U(1,:) = -U(2,:);
 U(n+2,:) = -U(n+1,:);
 
 % reverse direction at the y edges
 V(:,1) = -V(:,2);
 V(:,n+2) = -V(:,n+1);
 
 
 % First half step
 i = 1:n+1;
 j = 1:n+1;
 
 % height
 Hx(i,j) = (H(i+1,j+1)+H(i,j+1))/2 - dt/(2*dx)*(U(i+1,j+1)-U(i,j+1)); 
 Hy(i,j) = (H(i+1,j+1)+H(i+1,j))/2 - dt/(2*dy)*(V(i+1,j+1)-V(i+1,j));
 
 % x momentum

 Ux(i,j) = (U(i+1,j+1)+U(i,j+1))/2 - ...
 dt/(2*dx)*( U(i+1,j+1).^2./H(i+1,j+1) - U(i,j+1).^2./H(i,j+1) + ...
 g/2*H(i+1,j+1).^2 - g/2*H(i,j+1).^2 ...
 );
 
 Uy(i,j) = (U(i+1,j+1)+U(i+1,j))/2 - ...
 dt/(2*dy)*( (V(i+1,j+1).*U(i+1,j+1)./H(i+1,j+1)) - (V(i+1,j).*U(i+1,j)./H(i+1,j)) );
 
 % y momentum
 Vx(i,j) = (V(i+1,j+1)+V(i,j+1))/2 - ...
 dt/(2*dx)*((U(i+1,j+1).*V(i+1,j+1)./H(i+1,j+1)) - ...
 (U(i,j+1).*V(i,j+1)./H(i,j+1)));
 
 Vy(i,j) = (V(i+1,j+1)+V(i+1,j))/2 - ...
 dt/(2*dy)*((V(i+1,j+1).^2./H(i+1,j+1) + g/2*H(i+1,j+1).^2) - ...
 (V(i+1,j).^2./H(i+1,j) + g/2*H(i+1,j).^2));
 
 % Second half step
 i = 2:n+1;
 j = 2:n+1;
 
 % height
 H(i,j) = H(i,j) - (dt/dx)*(Ux(i,j-1)-Ux(i-1,j-1)) - (dt/dy)*(Vy(i-1,j)-Vy(i-1,j-1));


 % x momentum
 U(i,j) = U(i,j) - (dt/dx)*((Ux(i,j-1).^2./Hx(i,j-1) + g/2*Hx(i,j-1).^2) - ...
 (Ux(i-1,j-1).^2./Hx(i-1,j-1) + g/2*Hx(i-1,j-1).^2)) ...
 - (dt/dy)*((Vy(i-1,j).*Uy(i-1,j)./Hy(i-1,j)) - ...
 (Vy(i-1,j-1).*Uy(i-1,j-1)./Hy(i-1,j-1)));
 % y momentum
 V(i,j) = V(i,j) - (dt/dx)*((Ux(i,j-1).*Vx(i,j-1)./Hx(i,j-1)) - ...
 (Ux(i-1,j-1).*Vx(i-1,j-1)./Hx(i-1,j-1))) ...
 - (dt/dy)*((Vy(i-1,j).^2./Hy(i-1,j) + g/2*Hy(i-1,j).^2) - ...
 (Vy(i-1,j-1).^2./Hy(i-1,j-1) + g/2*Hy(i-1,j-1).^2));
 
end