% An orbital tug spacecraft is being considered to perform debris removal
% of spen launch vehicle upper stages. The debris are all 750 km circular
% orbit, with an inclination of 74 degrees

%a) Assume we launch our tug into a 200km circulear orbit. We launch a in
%plane, but out of phase with the debris of the object, trailing it by 120
%deg true anomaly. 
% i)Assuming we use the initial orbit as our phasing orbit,
%how long should we wait to perform our rendezvous initiation burn,
%injecting into a hohmann transfer to the debris orbit?


clear all, close all, clc, clf;



%% Period of Hohmann transfer

re = 6371.01; %km mean volume radius


%% Initialize system state variable with the initial conditions

R = zeros(6,1);

R1 = zeros(6,1);
R1(1) = 6578.1370; %km
R1(5) = -8.7235; %km/s
R1(6) = 4.6384; %km/s

R(1) = 6578.1370; %km
R(5) = 8.7235; %km/s
R(6) = 4.6384; %km/s

mu = 398600.4415; %km^3/s^2

SystemStateVectorLength = length(R);

%% Solve Differential Equation
tf = 25000;

opts = odeset('RelTol',1e-13,'AbsTol',1e-100);
[t,r] = ode45(@Syssim2,[0 tf],R,opts,mu);

[t,r1] = ode45(@Syssim2,[0 tf],R1,opts,mu);


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


%% Create position vectors for the individual componenets
rmag1 = zeros(length(t),1);
vmag1 = zeros(length(t),1);

rx1 = r1(:,1);
ry1 = r1(:,2);
rz1 = r1(:,3);

rvec1 = cat(2, rx1, ry1);
rvec1 = cat(2,rvec1,rz1);

vx1 = r1(:,4);
vy1 = r1(:,5);
vz1 = r1(:,6);

vvec1 = cat(2, vx1,vy1);
vvec1 = cat(2, vvec1, vz1);




%magnitudes of the initial conditions of the vector
for i = 1:length(t)
    
    rmag1(i) = sqrt(sum(rvec1(i,:).^2));
    vmag1(i) = sqrt(sum(vvec1(i,:).^2));
end


%% Play animation of simulation
close all
figure(1)
hold on

plot3(rx(:,1),ry(:,1),rz(:,1), 'r', rx1(:,1),ry1(:,1),rz1(:,1), 'k')
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

[i,j,k]=sphere(25);
earth=surf(i*re,j*re,k*re);
shading interp;
earth.FaceColor = 'b';
earth.EdgeColor = 'g'
earth.FaceAlpha = 0.2;
earth.EdgeAlpha = 0.2;

line([-re,re],[0,0],[0,0], 'LineWidth',1, 'Color', 'k');
line([0,0],[-re,re],[0,0], 'LineWidth',1, 'Color', 'k');
line([0,0],[0,0],[-re,re], 'LineWidth',1, 'Color', 'k');
text(0,0,re*1.1,'$\left| z \right>$', 'Interpreter', 'latex', ...
    'FontSize', 20, 'HorizontalAlignment', 'Center');
text(0,re*1.1,0,'$\left| y \right>$', 'Interpreter', 'latex', ...
    'FontSize', 20, 'HorizontalAlignment', 'Center');
text(re*1.1,0,0,'$\left| x \right>$', 'Interpreter', 'latex', ...
    'FontSize', 20, 'HorizontalAlignment', 'Center');

legend('s1','s2','earth');
grid on;
axis equal;
hold off;


%% Orbital Computations

%Specific Angular Momentum
hvec = cross(rvec,vvec);
hvec = hvec(1,:);
hmag = sqrt(sum(hvec.^2));

% Semi-Latus Rectum
P = (hmag^2)/mu;

% Semi-major axis
a = (2./rmag - ((vmag.^2)/mu)).^-1;
a = a(1);

% Eccentricity
emag = sqrt(1-(P/a(1)));

for i = 1:length(t)
    evec(i,:) = (1/mu)*cross(vvec(i,:),hvec) - rvec(i,:)./rmag(i);
end


%% Specific Mechanical Energy

Eta = zeros(length(t),1);
PE = zeros(length(t),1);
KE = zeros(length(t),1);

% Specific Kinetic and Potential Energy

KE = transpose(0.5.*vmag.*vmag);  %want this get the velocity vector over time

% Specific Potential Energy

PE = -transpose(mu./rmag);

Eta = KE+PE;

figure(2)
A = max(abs(rx));
plot(t(:),KE(:),'--',t(:),PE(:),':',t(:),Eta(:))
xlabel('Time (s)')
ylabel('Specific Energy (J/kg)')
legend('T','U','E')
axis tight
grid on

%% Keplarian Orbital Elements
%Semi-major axis, a. calculated with the apoapsis and periapsis

%need to print the results with their respective units
formatSpeca = 'The semi-major axis is %f kilometers\n';
fprintf(formatSpeca,a);

formatSpece = 'The eccentricity is %f\n';
fprintf(formatSpece,emag);

%Line of nodes:
khat = [0;0;1];

n = cross(khat,hvec);
nmag = sqrt(sum(n.^2));

%inclination
inc = acos(hvec(3)/hmag);

formatSpeci = 'The inclination is %f radians\n';
fprintf(formatSpeci,inc);

%Right ascension of the ascending node
omega = acos(n(1)/nmag);

if n(2) < 0
    omega = 2*pi - omega;
end

formatSpecomega = 'The right ascension is %g radians\n';
fprintf(formatSpecomega,omega);


%argument of perigee, w

% which value of evector should be used?
w = acos(dot(n,evec(2100,:))/(nmag*emag));

if evec(3) < 0
    w = 2*pi- w;
end

formatSpecw = 'The argument of perigee is %f radians\n';
fprintf(formatSpecw,w);


%true anomol, nu
test1 = dot(evec,rvec); %this is a 1x3
test2 = (emag.*rmag); %this is a 1x length(t)

nu = zeros(length(t),1);

for i = 1:length(t)
    nu(i) = real(acos(dot(evec(i,:),rvec(i,:))./(emag.*rmag(i))));
end

if dot(rvec,vvec) < 0
    nu = 2*pi- nu;
end

%% Orbital period

T = 2*pi*sqrt((a(1)^3)/mu);
T_hr = T/3600;

formatSpec4 = 'The orbital period is %f hours, ';
fprintf(formatSpec4,T_hr);
formatSpec5 = 'or %f seconds\n';
fprintf(formatSpec5,T);


%% Quarter Orbital Radius

%loop through the time vector and compare it to half of the period of
%orbit. If they are equal(within certain significant figures),set the half orbital radius at that time stamp.
for i = 1:length(t)
    if round(t(i),-1) == round(T/2,0)
        HalfOrbitalPeriod = i;
    end
end

HalfOrbitalRadius = rmag(HalfOrbitalPeriod);

formatSpecHOR = 'At a half of the orbital period, the spacecraft''s orbital radius %f kilometers\n';
fprintf(formatSpecHOR,HalfOrbitalRadius);


%% Propagate the satellise to one quarter of an orbital period

for i = 1:length(t)
    if round(t(i),-1) == round(T/4,0)
        QuarterOrbitalPeriod = i;
    end
end

QuaterOrbitalRadius = rmag(QuarterOrbitalPeriod);
QuaterOrbitalPosition = rvec(QuarterOrbitalPeriod,:);
QuarterOrbitalVelocity = vvec(QuarterOrbitalPeriod,:);

formatSpecQOR = 'At a quarter of the orbital period, the spacecraft''s orbital radius is %f Kilometers\n';
fprintf(formatSpecQOR,QuaterOrbitalRadius);

formatSpecQOP = 'At a quarter of the orbital period, the spacecraft''s position vector is [%e; %e; %e] kilometers\n';
fprintf(formatSpecQOP,QuaterOrbitalPosition);

formatSpecQOV = 'At a quarter of the orbital period, the spacecraft''s Velocity vector is [%f; %f; %f] in kilometers/second\n';
fprintf(formatSpecQOV,QuarterOrbitalVelocity);

%% Propagate using Kepler's Equation

deltat = (t(QuarterOrbitalPeriod)-t(1));

E0 = atan(tan(nu(QuarterOrbitalPeriod)/2)*(sqrt((1 +emag)/(1-emag))^-1))*2; % Eccentric Anomoly at IC
M0 = E0 - emag*sin(E0); %mean anomoly at IC

%get the initial guess of E at the quarter orbit from vallado book
if -pi <M0 < 0 ||  M0>pi
    Eq(1) = M0 - emag;
else
    Eq(1) = M0 + emag;
end

epsilon = .000001;
err = 1;

for i = 1:100
    Eq(i+1) = Eq(i)-(-M0+ Eq(i) - emag*sin(Eq(i)))/(1-emag*cos(Eq(i)));
    err=abs(Eq(i+1)-Eq(i));
    if err<epsilon
        break
    end
end

%Then we take the that value of E and solve it for Nu
nuq = atan(sqrt((1+emag)/(1-emag))*tan(Eq(length(Eq))/2))*2;

% This will give us the true anomoly at a quarter orbital position, from
% which we can calculate the position and velocity vectors

rpqw(1) = P*cos(nuq)/(1 + emag*cos(nuq));
rpqw(2) = P*sin(nuq)/(1 + emag*cos(nuq));
rpqw(3) = 0;

rpqw = transpose(rpqw);

vpqw(1) = -sqrt(mu/P)*sin(nuq);
vpqw(2) = sqrt(mu/P)*(emag +cos(nuq));
vpqw(3) = 0;

vpqw = transpose(vpqw);

% Setup the 313 rotational transformation
R1 = [1 0 0; 0 cos(-inc) sin(-inc); 0 -sin(-inc) cos(-inc)];
R3o = [cos(-omega) -sin(-omega) 0; sin(-omega) cos(-omega) 0; 0 0 1];
R3w = real([cos(-w) -sin(-w) 0; sin(-w) cos(-w) 0; 0 0 1]);


Transform = R3o*R1*R3w;

rijk = Transform*rpqw;
vijk = Transform*vpqw;


formatSpecKepQOP = 'Propagating using Kepler''s Equation, at a quarter of the orbital period the spacecraft''s position vector is [%e; %e; %e] kilometers\n';
fprintf(formatSpecKepQOP,rijk);

formatSpecKepQOV = 'Propagating using Kepler''s Equation, at a quarter of the orbital period the spacecraft''s Velocity vector is [%f; %f; %f] in kilometers/second\n';
fprintf(formatSpecKepQOV,vijk);
