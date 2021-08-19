clear all;  
clc;
%% Radius of Earth's SOI

a = 149597871; %Earth's semimajor axis (1AU) or in kms

muS = 1.327e11; %km3/s2

muE = 3.986e5; %km3/s2

SOI = a*((muE/muS)^(2/5))

%% Lambert targeting
%July 4 2021

%assume the earth is on a circular orbit, with radius a


rp = a;
ra = 1.3*a;

dVtrans = sqrt(muS/rp)*abs(sqrt((2*ra)/(ra+rp))-1);

C3 = dVtrans^2 %km^2/s^2


%first calculate the velocity of earth with respect to the sun

%the calculate the C3 based on the relative velocity and the apo and
%perihelion

%% Part c

rMD = 1.05*a; %AU

ai = (rp+ra)/2;

e = 1-rp/ai;

P = ai*(1-e*e);

v = sqrt(2*muE/rMD - muE/ai);

h = sqrt(P*muE);

phi = rad2deg(acos(h/(rMD*v)))


vellipse = sqrt(2*muS/rMD - muS/ai)
vcircle = sqrt(muS/rMD)


vt = vcircle-vellipse*cos(deg2rad(phi))

vr = vellipse*sin(deg2rad(phi))

VooIncoming = sqrt(vt*vt+vr*vr)

%% Propellant Mass
mf = 1000 %kg

ISP = 312 %s

%%need to define gravity in km/s^s

g = 9.81/1000;

%%need to make sure i have the right deltaV


m0 = mf*exp(VooIncoming/(g*ISP))

mp = m0-mf
%% Orion launch and entry

Vei = 12.4; %(km/s)

%assuming Hohmann transfer from Mars to Earth, what is the hperbolic excess
%velocity with which the capsule arrives at Earth's SOI

% Re = 149597871;

Re = 150e6;
Rm = 228e6; %km
muM = 4.283e4; %km3/s2

Ves = sqrt(muS/Re);
Vms = sqrt(muS/Rm);

Voodep = Vms*(sqrt((2*Re)/(Re+Rm))-1);

C3 = Voodep^2;

%rpo = 3396+200; %parking orbit for mars
rei = 6503;
%rei = 3522 %mars atmospheric radius

%deltaVdep = sqrt(C3+(2*muM)/rpo) - sqrt(muM/rpo);

Vooarr = abs(sqrt(muS/Re)*(1-sqrt(2*Rm/(Re+Rm))))


%Vei = sqrt(C3+(2*muE)/rei)


Voomax = sqrt(Vei*Vei - (2*muE)/rei)

deceleration = 0;

if Vooarr>Voomax
    deceleration = true
    deltaVRET = Vooarr-Voomax
else
    deceleration = false
end


