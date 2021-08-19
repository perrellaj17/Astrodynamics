% An orbital tug spacecraft is being considered to perform debris removal
% of spen launch vehicle upper stages. The debris are all 750 km circular
% orbit, with an inclination of 74 degrees
%Assume point masses

%a) Assume we launch our tug into a 200km circulear orbit. We launch a in
%plane, but out of phase with the debris of the object, trailing it by 120
%deg true anomaly. 
% i)Assuming we use the initial orbit as our phasing orbit,
%how long should we wait to perform our rendezvous initiation burn,
%injecting into a hohmann transfer to the debris orbit?


clear all, close all, clc;


%% 1 Compute mu for the CRTBP
me = 5.597219*(10^24)
mm = 7.349*(10^22)

MU = mm/(me+mm)
Comp = abs(MU-0.012277471)
6.8074e-04

%% Period of Hohmann transfer
mu = 398600.435436; %km^3/s^2 from jpl horizons
re = 6371.01; %km mean volume radius
rtgt = re + 750; %km
Ttgt = 2*pi*sqrt((rtgt^3)/mu);
rint = re + 200; %km
atrans = (rint+rtgt)/2;


Ttrans = pi*sqrt(((atrans)^3)/(mu))


%% Time wait
wtgt = sqrt(mu/(rtgt^3));
wint = sqrt(mu/(rint^3));

closingrate = wint-wtgt;

alphaL = wtgt*Ttrans;

nu = alphaL-pi;
nui = deg2rad(-120);

deltanu = nu-nui;

if deltanu<0
    k=1;
else
    k=0;
end

Twait = (deltanu+2*pi*k)/closingrate %s

% ii) what is the delta v required for this burn, as well as the
% circularizing burn at rendevous?


deltavcircularize = abs(sqrt((2*mu/rtgt)-(mu/atrans))-sqrt(mu/rtgt))
deltavinitialize = abs(sqrt((2*mu/rint)-(mu/atrans))-sqrt(mu/rint));


deltavtot = deltavinitialize+deltavcircularize %km/s


%% b) Assume debris is in perfectly circular orbit, and no external forces imparted by interceptor satellite
%We rendezvous at a point 2km behind the debris in-track. We wish to 
%perform an R-bar approach, setting up a debris capture opportunity
%approximately one and a half orbits later, at a relative position and
%velocity of:

deltar = [-4; 0; 0]; %m
deltardot = [-0.04; 0.05; 0]; %m/s

% in RSW coordinates. We will perform a series of maneuvers between the
% waypoints given in the table below. The time of flight should be
% determined to match the desired relative state given above.



%need to convert the previous delta vectors into RIC

%% DeltaV
%determine the deltaV required for each maneuver using the HCW targeting
%equation

wri = wtgt;
wrivec = [0;0;wri];


%define system state vector:

% dR0 = [drr0, dri0, drc0, drdotr0, drdoti0, drdotc0]

%define state transition matrix
t = zeros(7,1);
t(1) = 0;
t(2) = 1*Ttgt;
t(3) = Ttgt/4;
t(4) = 0.06*Ttgt;
t(5) = 0.06*Ttgt;
t(6) = 0.06*Ttgt;
t(7) = 0.025*Ttgt;

dri = zeros(7,3);
dri(1,1:3) = [0 2000 0];
dri(2,1:3) = [0 -600 0];
dri(3,1:3) = [-300 0 0];
dri(4,1:3) = [-100 0 0];
dri(5,1:3) = [-50 0 0];
dri(6,1:3) = [-10 0 0];
dri(7,1:3) = [-4 0 0];


drdotp=zeros(6,3);


deltaV=zeros(6,3);


for k=1:6
    c = cos(wri*t(k+1));
    s = sin(wri*t(k+1));
    
    PHI = [4-3*c 0 0 s/wri (2/wri)*(1-c) 0;
        6*(s-wri*t(k+1)) 1 0 -(2/wri)*(1-c) (4*s/wri)-3*t(k+1) 0;
        0 0 c 0 0 s/wri;
        3*wri*s 0 0 c 2*s 0;
        -6*wri*(1-c) 0 0 -2*s 4*c-3 0;
        0 0 -wri*s 0 0 c];
    
    for i=1:3
        for j=1:3
            phirr(i,j) = PHI(i,j);
            phirrdot(i,j) = PHI(i,j+3);
            phirdotr(i,j) = PHI(i+3,j);
            phirdotrdot(i,j) = PHI(i+3,j+3);
        end
    end
    
    drdotp(k,1:3)= transpose(inv(phirrdot)*(transpose(dri(k+1,1:3))-phirr*transpose(dri(k,1:3))));
   
    
    if k==1
        dricheck(k,1:3) = (phirr*transpose(dri(k,1:3)));
    end
    
    dricheck(k+1,1:3) = round(phirrdot*transpose(drdotp(k,1:3))+(phirr*transpose(dri(k,1:3))));
    
end


drdotp(7,1:3)=phirr*transpose(dri(6,1:3)) + phirdotrdot*transpose(drdotp(6,1:3)) 


for i=1:6
    
        deltaV(i,1:3)= drdotp(i+1,1:3)-drdotp(i,1:3);
    
end

drdotp
deltaV
dricheck






time = t(2);
dt=t(2)/100;

dr1 = zeros(time/dt,3);
inc=0;
for f = 0:dt:time
    inc = inc+1;
    c = cos(wri*f);
    s = sin(wri*f);
    
    
    PHI = [4-3*c 0 0 s/wri (2/wri)*(1-c) 0;
        6*(s-wri*f) 1 0 -(2/wri)*(1-c) (4*s/wri)-3*f 0;
        0 0 c 0 0 s/wri;
        3*wri*s 0 0 c 2*s 0;
        -6*wri*(1-c) 0 0 -2*s 4*c-3 0;
        0 0 -wri*s 0 0 c];
    
    for i=1:3
        for j=1:3
            phirr(i,j) = PHI(i,j);
            phirrdot(i,j) = PHI(i,j+3);
            phirdotr(i,j) = PHI(i+3,j);
            phirdotrdot(i,j) = PHI(i+3,j+3);
        end
    end
    
    dr1(inc,1:3) = phirr*transpose(dri(1,1:3))+phirrdot*transpose(drdotp(1,1:3));
    
end


time = t(3);
dt=t(3)/100;

dr2 = zeros(time/dt,3);
inc=0;
for f = 0:dt:time
    inc = inc+1;
    c = cos(wri*f);
    s = sin(wri*f);
    
    
    PHI = [4-3*c 0 0 s/wri (2/wri)*(1-c) 0;
        6*(s-wri*f) 1 0 -(2/wri)*(1-c) (4*s/wri)-3*f 0;
        0 0 c 0 0 s/wri;
        3*wri*s 0 0 c 2*s 0;
        -6*wri*(1-c) 0 0 -2*s 4*c-3 0;
        0 0 -wri*s 0 0 c];
    
    for i=1:3
        for j=1:3
            phirr(i,j) = PHI(i,j);
            phirrdot(i,j) = PHI(i,j+3);
            phirdotr(i,j) = PHI(i+3,j);
            phirdotrdot(i,j) = PHI(i+3,j+3);
        end
    end
    
    dr2(inc,1:3) = phirr*transpose(dri(2,1:3))+phirrdot*transpose(drdotp(2,1:3));
    
end

time = t(4);
dt=t(4)/100;

dr3 = zeros(time/dt,3);
inc=0;
for f = 0:dt:time
    inc = inc+1;
    c = cos(wri*f);
    s = sin(wri*f);
    
    
    PHI = [4-3*c 0 0 s/wri (2/wri)*(1-c) 0;
        6*(s-wri*f) 1 0 -(2/wri)*(1-c) (4*s/wri)-3*f 0;
        0 0 c 0 0 s/wri;
        3*wri*s 0 0 c 2*s 0;
        -6*wri*(1-c) 0 0 -2*s 4*c-3 0;
        0 0 -wri*s 0 0 c];
    
    for i=1:3
        for j=1:3
            phirr(i,j) = PHI(i,j);
            phirrdot(i,j) = PHI(i,j+3);
            phirdotr(i,j) = PHI(i+3,j);
            phirdotrdot(i,j) = PHI(i+3,j+3);
        end
    end
    
    dr3(inc,1:3) = phirr*transpose(dri(3,1:3))+phirrdot*transpose(drdotp(3,1:3));
    
end
dr3;


time = t(5);
dt=t(5)/100;

dr4 = zeros(time/dt,3);
inc=0;
for f = 0:dt:time
    inc = inc+1;
    c = cos(wri*f);
    s = sin(wri*f);
    
    
    PHI = [4-3*c 0 0 s/wri (2/wri)*(1-c) 0;
        6*(s-wri*f) 1 0 -(2/wri)*(1-c) (4*s/wri)-3*f 0;
        0 0 c 0 0 s/wri;
        3*wri*s 0 0 c 2*s 0;
        -6*wri*(1-c) 0 0 -2*s 4*c-3 0;
        0 0 -wri*s 0 0 c];
    
    for i=1:3
        for j=1:3
            phirr(i,j) = PHI(i,j);
            phirrdot(i,j) = PHI(i,j+3);
            phirdotr(i,j) = PHI(i+3,j);
            phirdotrdot(i,j) = PHI(i+3,j+3);
        end
    end
    
    dr4(inc,1:3) = phirr*transpose(dri(4,1:3))+phirrdot*transpose(drdotp(4,1:3));
    
end
dr4;


time = t(6);
dt=t(6)/100;

dr5 = zeros(time/dt,3);
inc=0;
for f = 0:dt:time
    inc = inc+1;
    c = cos(wri*f);
    s = sin(wri*f);
    
    
    PHI = [4-3*c 0 0 s/wri (2/wri)*(1-c) 0;
        6*(s-wri*f) 1 0 -(2/wri)*(1-c) (4*s/wri)-3*f 0;
        0 0 c 0 0 s/wri;
        3*wri*s 0 0 c 2*s 0;
        -6*wri*(1-c) 0 0 -2*s 4*c-3 0;
        0 0 -wri*s 0 0 c];
    
    for i=1:3
        for j=1:3
            phirr(i,j) = PHI(i,j);
            phirrdot(i,j) = PHI(i,j+3);
            phirdotr(i,j) = PHI(i+3,j);
            phirdotrdot(i,j) = PHI(i+3,j+3);
        end
    end
    
    dr5(inc,1:3) = phirr*transpose(dri(5,1:3))+phirrdot*transpose(drdotp(5,1:3));
    
end
dr5;

time = t(7);
dt=t(7)/100;


dr6 = zeros(round(time/dt),3);
inc=0;
for f = 0:dt:time
    inc = inc+1;
    c = cos(wri*f);
    s = sin(wri*f);
    
    
    PHI = [4-3*c 0 0 s/wri (2/wri)*(1-c) 0;
        6*(s-wri*f) 1 0 -(2/wri)*(1-c) (4*s/wri)-3*f 0;
        0 0 c 0 0 s/wri;
        3*wri*s 0 0 c 2*s 0;
        -6*wri*(1-c) 0 0 -2*s 4*c-3 0;
        0 0 -wri*s 0 0 c];
    
    for i=1:3
        for j=1:3
            phirr(i,j) = PHI(i,j);
            phirrdot(i,j) = PHI(i,j+3);
            phirdotr(i,j) = PHI(i+3,j);
            phirdotrdot(i,j) = PHI(i+3,j+3);
        end
    end
    
    dr6(inc,1:3) = phirr*transpose(dri(6,1:3))+phirrdot*transpose(drdotp(6,1:3));
    
end
dr6;


plot(dr1(1:101,1),dr1(1:101,2),dr2(1:101,1),dr2(1:101,2),dr3(1:101,1),dr3(1:101,2),dr4(1:101,1),dr4(1:101,2),dr5(1:101,1),dr5(1:101,2),dr6(1:101,1),dr6(1:101,2))
xlabel('r (m)')
ylabel('i (m)')
title('Tug Path')
axis equal
hold on

r = 100;
theta = linspace(0,2*pi);
x = r*cos(theta);
y = r*sin(theta);
plot(x, y);
hold off

%% What is the total propellant mass required for this R-bar appoach?

% make sure to convert the dV back to km/s

g = 9.81; %m/s

m = zeros(7,1);
m(7) = 450; %kg
Isp = 300; %seconds
md = 1400; %kg
deltaVmag =zeros(8,1);

deltaVmag(1) = deltavinitialize/1000;
deltaVmag(2) = deltavcircularize/1000;


for i = 6:-1:1
    
    deltaVmag(i+2) = sqrt((deltaV(i,1)^2)+(deltaV(i,2)^2)+(deltaV(i,3)^2));
    
    %m(i) = m(i+1)*exp(deltaVmag(i)/(g*Isp));
end

m; %kg
deltaVmag %m/s



