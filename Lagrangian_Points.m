clear all, close all, clc;


%% 1 Compute mu for the CRTBP
%all values taken from JPL Horizons

me = 5.597219*(10^24); %kg
mm = 7.349*(10^22); %kg

MU = mm/(me+mm);
Comp = abs(MU-0.012277471);

gme = 398600.435436;

gmm = 4902.800066;

gMU = gmm/(gme+gmm);
gComp = abs(gMU-0.012277471)


%% Initialize system state variable with the initial conditions

mu = 0.012277471; %mass ratio

R = zeros(6,1);
R(1) = 1.061692; %DU
R(5) = 0.403877; %SU

SystemStateVectorLength = length(R);

%% Solve Differential Equation
tf = 25;

opts = odeset('RelTol',1e-13,'AbsTol',1e-100);
[t,r] = ode45(@Syssim4,[0 tf],R,opts,mu);

%% Create position vectors for the individual componenets
rmag = zeros(length(t),1);
vmag = zeros(length(t),1);

rx = r(:,1);
ry = r(:,2);
rz = r(:,3);

rvec = cat(2, rx, ry);
rvec = cat(2,rvec,rz);

vx = r(:,4);
vy = r(:,5);
vz = r(:,6);

vvec = cat(2, vx,vy);
vvec = cat(2, vvec, vz);

%magnitudes of the initial conditions of the vector
for i = 1:length(t)
    
    rmag(i) = sqrt(sum(rvec(i,:).^2));
    vmag(i) = sqrt(sum(vvec(i,:).^2));
end





%% Compute the x coordinates of the L1, L2, and L3 Lagrangian points

%Use of the fzero function in MATLAB is recommended

L1 = @l1;
x0 = 5;
z = fzero(L1,x0);
x1 = 1-mu-z;


L2 = @l2;
x0 = 5;
z = fzero(L2,x0);
x2 = 1-mu+z;

L3 = @l3;
x0 = 5;
z = fzero(L3,x0);
x3 = -mu-z;

x4 = 0.5-mu;
y4 = sqrt(3)/2;
x5 = x4;
y5 = -y4;

%% Play animation of simulation
close all
figure(1)

plot(rx(:,1),ry(:,1), 'k', 0, 0,  '*', 1, 0,  '*')
xlabel('x (DU)')
ylabel('y (DU)')
legend('Orbital Path','Earth','Moon');
grid on;
axis equal;

Lagrange = zeros(5,2);
Lagrange(1,1) = x1;
Lagrange(2,1) = x2;

Lagrange(3,1) = x3;

Lagrange(4,1) = x4;
Lagrange(4,2) = y4;
Lagrange(5,1) = x5;
Lagrange(5,2) = y5;

DU = 384400;

rm = 1737.53 %km
re = 6378.137 %km

rm = rm/DU; %now nondimensionalize
re = re/DU; %now nondimensionalize

theta = linspace(0,2*pi) ;
xe = re*cos(theta);
ye = re*sin(theta);

xm = rm*cos(theta)+1;
ym = rm*sin(theta);


figure(2)
plot( xe, ye, 'b', xm, ym,'k' , x1, 0, '*',x2, 0, '*',x3, 0, '*',x4, y4, '*',x5, y5, '*')
xlabel('x (DU)')
ylabel('y (DU)')

legend('Earth','Moon', 'L1', 'L2', 'L3','L4','L5');
grid on;
axis equal;
