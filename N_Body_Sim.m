%Create a general n-body system simulator that can accept either cartesian
%or spherical coordinas for the specification of initial conditions for an
%arbitrary number if bodies, using the conventions presented in lecture.
%Also allow the value for the universal gravitational constant, G, to be
%specified by the user. Use Planck Units (G = 1) and choose numerical
%values for both the masses and initial conditions of your bodies
%accordingly


%% a. Determine an appropriate set of transformations to convert the spherical
%coordinates initial conditions and present your results. Show diagrams as
%necessary
clear all, close all, clc;

r = 10;
ele = pi/2;
azi = pi/2;

x = r*sin(ele)*cos(azi);
y = r*sin(ele)*sin(azi);
z = r*cos(ele);

%This 3x1 is the vector from rst origin to ijk origin
rstPijk = [-x;-y;-z];

%using euler angles to rotate from the rst frame to the ijk frame first
%about s by the angle of elevation, and then about t the angle of azimuth
Rs = [cos(-ele) 0 sin(-ele); 0 1 0; -sin(-ele) 0 cos(-ele)];
Rt = [cos(-azi) -sin(-azi) 0; sin(-azi) cos(-azi) 0; 0 0 1];
Rr = [1 0 0; 0 1 0; 0 0 1];

%compute the equivalent rotation matrix
rstRijk = Rs*Rt*Rr;

%compute the transform
rstTijk = [rstRijk rstPijk; 0 0 0 1];


%% b. Using experimentation, develop a set of initial conditions for at least
% 3 bodies that yields an interesting set of trajectories for the bodies 
%over time. Ensure that none of the bodies escapes. Present your initial 
%conditions and the masses of the bodies in a table. Provide the following 
%data plots (plot some quantities on the same figure where appropriate;
%use your judgment):
%   1. Positions of each body in 3-space over time
%   2. Conter of mass position in 3-space over time
%   3. Center of mass velocity magnitude vs. time
%   4. Total system kinetic energy vs. time
%   5. Total system potential energy vs. time
%   6. Total system energy vs. time
%   7. Change in total system energy (relative to initial value) vs. time
%   8. Second time derivative of the system's moment of inertia vs. time

%assume that the double variable notation is a derivative with respect to
%time:  r is position, rr is velocity, rrr is acceleration in each
%respective unit vector direction, x, y, & z or in azimuth (az), elevation
%(ele), and distance (r).

%Assume point masses

%% Prompt User

prompt = 'Enter a value for the Universal Gravitational Constant: ';
G = input(prompt);
fprintf('Universal Gravitational Constant (G) = ')
disp(G)

%Take the input of n bodies, and then initialize an array with length n*
%the 6 state variables

prompt2 = 'Enter a value (n) for the number of bodies in the general system simulator: ';
n = input(prompt2);
formatSpec = 'This is a %d-body problem\n';
fprintf(formatSpec,n)

%% Initialize system state variable and an array for mass
R = [];

m = ones(n,1);
totalmass = sum(m);

for i = 1:n
    temparray = zeros(6,1);
    R = cat(1, R, temparray);
    
end
i = 0;

SystemStateVectorLength = length(R);

for i = 1:SystemStateVectorLength
    %This will randomize the initial position and velocity of the state
    %variables in the system state variable. In the current state, the
    %randomizer is -5 to 5
    R(i) = -5 + 10*rand;
end
i = 0;


%% The following code gives us the individual arrays which contain the
%individual position and velocities of the bodies by directional component
rx = [];
ry = [];
rz = [];
vx = [];
vy = [];
vz = [];

for i = 1:n
    temparray = zeros(1,1);
    rx = cat(1, rx, temparray);
    ry = cat(1, ry, temparray);
    rz = cat(1, rz, temparray);
    
    vx = cat(1, rx, temparray);
    vy = cat(1, ry, temparray);
    vz = cat(1, rz, temparray);
    
    rx(i) = R(1 + 6*(i-1));
    ry(i) = R(2 + 6*(i-1));
    rz(i) = R(3 + 6*(i-1));
    vx(i) = R(4 + 6*(i-1));
    vy(i) = R(5 + 6*(i-1));
    vz(i) = R(6 + 6*(i-1));
    
end

rx
example = rx - rx(1)
magnitudes = 
for i = 1:n
    rxmag(i) = sqrt(rx(i) + rx(i+1))
    
i=0;

%% Calculate the center of mass of the system
xcm =0;
ycm = 0;
zcm = 0;

for i = 1:n
    xcm = (xcm + m(i)*rx(i))/totalmass;
    ycm = (ycm + m(i)*ry(i))/totalmass;
    zcm = (zcm + m(i)*rz(i))/totalmass;
end
i=0;
cm = [xcm;ycm;zcm];


%% Solve Differential Equation

% opts = odeset('RelTol',1e-11,'AbsTol',1e-11);
% [t,r] = ode45(@shmsim,[0 tf],R,opts,G,m,rx,ry,rz);


%% Kinetic Energy
% T = .5*sum(m)*sum(v.^2);


%% Potential Energy

% U = -0.5 sumi sumj mi mj / rij


%% Total System Energy

% E = T - U


%% Change in total system energy (relative to initial value) vs. time



%% Second time derivative of the system's moment of inertia vs. time

% The systems total angular momentum is Hvec = sum mi*ri cross rri
% the moment of inertia is its angular momentum divided by its angular
% velocity

