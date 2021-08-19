close all, clear all, clc ;

%define global constants:
tf = 86400; %length of the simulation
tday = 86400; %number of seconds in a day


re = 6378.1363; %mean equitorial radius of the Earth km
ee = 0.081819221456; %the eccentricity of the earth
we = 7.292115*(10^(-5)); %rad/second, Earth's rotational velocity

rho = 2.11164*(10^9); %kg/m^3
rinit = re+880; %Initial Rsat km


Cd = 2.2; %Coefficient of drag (dimensionless)
msat = 1.33; %Cubesat mass kg
A = 1*(10^(-8)); %Exposed cross sectional area in km^2 (100cm^2)
bc = msat/(Cd*A); %Ballistic coefficient

mu = 398600.4415; %km^3/s^2


thetainit = -76.937785;
thetainit = deg2rad(thetainit);
phiinit = 38.988542;
phiinit = deg2rad(phiinit);

xmdinit = re*sin(phiinit)*cos(thetainit);
ymdinit = re*sin(phiinit)*sin(thetainit);
zmdinit = re*cos(phiinit);

xinit = rinit*sin(phiinit)*cos(thetainit);
yinit = rinit*sin(phiinit)*sin(thetainit);
zinit = rinit*cos(phiinit);


%% Initialize system state variable with the initial conditions

R = zeros(6,1);

R(1) = xinit; %km
R(2) = yinit; %km
R(3) = zinit; %km


R(4) = 5.2; %km/s
R(5) = 5.2; %km/s
R(6) = 2.0; %km/s

SystemStateVectorLength = length(R);


%% Solve Differential Equation


%IN ORDER TO SIMULATE MULTIPLE SATELLITES AT THE SAME TIME, YOU MUST
%COMPUTE ALL OF THEIR STATE VECTORS IN THE SAME ODE45 NUMERICAL
%INTEGRATION. THAT WAY WE KNOW THE STATE AT ALL OF THE SAME TIME STEPS AND
%DO NOT HAVE TO ACCOUNT FOR DIFFERENT TIMESTEP SIZES

opts = odeset('RelTol',1e-6,'AbsTol',1e-10);
[t,r] = ode45(@Syssim2,[0 tf],R,opts,mu);


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


for i =1:length(t)
    rmag(i) = sqrt(sum(rvec(i,:).^2));
    vmag(i) = sqrt(sum(vvec(i,:).^2));
end
AverageAltitude = mean(rmag)-re

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

evec(1,:) = (1/mu)*cross(vvec(1,:),hvec) - rvec(1,:)./rmag(1);

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
w = acos(dot(n,evec(1,:))/(nmag*emag));

if evec(3) < 0
    w = 2*pi- w;
end

formatSpecw = 'The argument of perigee is %f radians\n';
fprintf(formatSpecw,w);

%% Energy Computations

% Specific Kinetic and Potential Energy

KE = transpose(0.5.*vmag.*vmag);  %want this get the velocity vector over time

% Specific Potential Energy

PE = -transpose(mu./rmag);

Eta = KE+PE;

figure(1)
A = max(abs(rx));
plot(t(:),KE(:),'--',t(:),PE(:),':',t(:),Eta(:))
xlabel('Time (s)')
ylabel('Specific Energy (J/kg)')
legend('T','U','E')
axis tight
grid on


%% Computing Latitude and Longitude
rdsat= zeros(length(t),1);
alpha = zeros(length(t),1);
delta = zeros(length(t),1);

C = zeros(length(t),1);
phigd = zeros(length(t),1);
phigc = zeros(length(t),1);
hellp = zeros(length(t),1);
elevation = zeros(length(t),1);

That = zeros(length(t),3);
Cross = zeros(length(t),3);
Nhat = zeros(length(t),3);
What = zeros(length(t),3);


for i = 1:length(t)
    %rdsat is the equitorial projection of the satellites position vector
    rdsat(i) = sqrt(rx(i)^2+ry(i)^2);
    
    %compute the right ascension (alpha)
    alpha(i) = rad2deg(asin(ry(i)/rdsat(i)));
    
    %determine delta as the starting value for the iteration, we can use
    %the position vector as a rough guess because the declination and
    %geocentric latitude are equal
    delta(i) = rad2deg(atan(rz(i)/rdsat(i)));
    
    %set phigd to delta and rd to rdsat, and rk to rksat 
    
    tolerance = 10;
    
    while tolerance > 0.0001
        
        dold = delta(i);
        C(i) = re/sqrt(1-(ee*ee*sin(delta(i))*sin(delta(i))));
        delta(i) = rad2deg(atan((rz(i)+C(i)*ee*ee*sin(delta(i)))/rdsat(i)));
        tolerance =abs(dold-delta(i));
    end
    
    %Produce the geodetic and geocentric latitude
    phigd(i) = deg2rad(delta(i));
    
    phigc(i) = atan((1-ee^2)*tan(phigd(i)));
    %height of the satellite from it's groundtack traced by the geodetic
    %latitude
    hellp(i) = (rdsat(i)/cos(phigd(i)))-C(i);
    
    
    
    %Compute the NTW coordinate frame
    That(i,:) = vvec(i,:)/vmag(i);
    Cross(i,:) = cross(rvec(i,:),vvec(i,:));
    What(i,:) = Cross(i,:)/sqrt(sum(Cross(i,:).^2));
    Nhat(i,:) = cross(That(i,:), What(i,:));
    
end


%% UMD rotates around the globe


theta = zeros(length(t),1);
mdx = zeros(length(t),1);
mdy = zeros(length(t),1);
mdz = zeros(length(t),1);
mdmag = zeros(length(t),1);

satx = zeros(length(t),1);
saty = zeros(length(t),1);
satz = zeros(length(t),1);
satmag = zeros(length(t),1);

long = zeros(length(t),1);

% perpx = zeros(length(t),1);
% perpy = zeros(length(t),1);
% perpz = zeros(length(t),1);
% perpmag = zeros(length(t),1);
% elevation= zeros(length(t),1);





for i = 1:length(t)
    
    %move the location of maryland based on time elapsed. The university
    %should rotate 2Pi every 24 hours or 86400 seconds
    
    deltat = t(i)/tday;
    theta(i) = thetainit + deltat*2*pi;
    
    
    %compute the longitude of the satellite
    long(i)=atan2(ry(i),rx(i))-deltat*2*pi;
    
  
    while long(i)<(-pi)
        long(i) = long(i)+2*pi;
    end
    

    % Now compute a vector from the origin to maryland at each timestep as
    % maryland rotates on the surface of the globe
    
    mdx(i,1) = re*sin(phiinit)*cos(theta(i));
    mdy(i,1) = re*sin(phiinit)*sin(theta(i));
    mdz(i,1) = re*cos(phiinit);
    mdmag(i) = sqrt(mdx(i)^2 + mdy(i)^2 + mdz(i)^2);
    
    % Now compute a vector from the location of maryland to each of
    % the satellites at every time step
 
    satx(i) = rx(i,1);
    saty(i) = ry(i,1);
    satz(i) = rz(i,1);
    satmag(i) = sqrt(satx(i)^2 + saty(i)^2 + satz(i)^2);
    
   
    
end


%% Compute line of sight (LOS)
 mdvec = cat(2,mdx,mdy);
 mdvec = cat(2, mdvec, mdz);
   
 sight = zeros(length(t),1);
 LOS = zeros(length(t),1);
 check = zeros(length(t),1);
 
for i=1:length(t)
    sight(i) = (mdmag(i)*mdmag(i) - dot(mdvec(i,:),rvec(i,:)))/(mdmag(i)*mdmag(i)+rmag(i)*rmag(i)-2*dot(mdvec(i,:),rvec(i,:)));
    if sight(i)<0 || sight(i)>1
        LOS(i) = 1;
    else
        check(i) = ((1-sight(i))*mdmag(i)*mdmag(i)+dot(mdvec(i,:),rvec(i,:))*sight(i))/(re*re);
        if check(i)>=1
            LOS(i) = 1;
        end
    end
end

%% Compute revisit and coverage time

totalRevisit = 0;
avgRevisit = 0;

totalCoverage = 0;
avgCoverage = 0;

current = 1;

for i = 1:(length(t)-1)
    change = LOS(i+1)-LOS(i);
    if change == -1
        last = i;
        
        if totalCoverage == 0
            totalCoverage = t(last)-t(current);
            avgCoverage = t(last)-t(current);
        else 
            coverageCheck = t(last) - t(current);
            totalCoverage = totalCoverage + (t(last)-t(current));
            avgCoverage = (avgCoverage + (t(last)-t(current)))/2;
        end
    end
    if change == 1
        current = i+1;
       
        if totalRevisit == 0
            totalRevisit = t(current)-t(last);
            avgRevisit = t(current)-t(last);
        else 
            revisitCheck = (t(current)-t(last));
            totalRevisit = totalRevisit + (t(current)-t(last));
            avgRevisit = (avgRevisit + (t(current)-t(last)))/2;
        end
    end
end

%This section accounts for the final segment that is not in the previous
%loop. This segment should not be included in the final code
% avgCoverage
% totalCoverage
% 
% avgRevisit
% totalRevisit
% 
% if last>current
%     totalRevisit = totalRevisit + t(length(t))-t(last)
%     avgRevisit = (avgRevisit + (t(length(t))-t(last)))/2
% else
%     totalCoverage = totalCoverage +t(length(t))-t(current)
%     avgCoverage = (avgCoverage +t(length(t))-t(current))/2
% end



%this way we only compute averages for the values we KNOW. if we use the
%commented loop from above, we assume that the end of the time array is a
%switch in the LOS, which isnt necessarily true.

TotalTime = totalRevisit+totalCoverage;
PercentCoverage = (totalCoverage/TotalTime)*100
PercentUncovered = 100-PercentCoverage



%% Plot the location of the satellite
loc = round(length(t)/4,0);

p0 = [0 0 0];

%define the line from the center of the ECI frame to the university of
%marylalnd at that particular timestep
pmd = [mdx(loc), mdy(loc), mdz(loc)];
mdline = [p0;pmd];

%define the line from the university of maryland to the satellite
% at that particular timestep
psat = [satx(loc), saty(loc), satz(loc)];
satline = [pmd;psat];


%define the line from the center of the ECI frame to the university of
%marylalnd at that particular timestep
% perp= [perpx(loc), perpy(loc), perpz(loc)];
% perpline = [pmd;perp];


%Add the NWT frame at the location of the satellite at that particular
%timestep
nwtmag = 100;
nvar=[Nhat(loc,1),Nhat(loc,2),Nhat(loc,3)];
nline = [psat;psat+nwtmag*nvar];

W=[What(loc,1),What(loc,2),What(loc,3)];
Wline = [psat;psat+nwtmag*W];

T=[That(loc,1),That(loc,2),That(loc,3)];
Tline = [psat;psat+nwtmag*T];


figure(2)
hold on
plot3(rx(:,1),ry(:,1),rz(:,1), 'k', xmdinit, ymdinit, zmdinit, '*', mdx(:,1), mdy(:,1),mdz(:,1), 'y',...
    mdline(:,1), mdline(:,2), mdline(:,3), 'r', satline(:,1), satline(:,2), satline(:,3), 'c',...
    nline(:,1), nline(:,2), nline(:,3), 'r',Wline(:,1), Wline(:,2), Wline(:,3), 'b',Tline(:,1), Tline(:,2), Tline(:,3), 'g',...
    satx(loc), saty(loc), satz(loc), 'o');
axis equal;
grid on;

xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

[i,j,k]=sphere(25);
earth=surf(i*re,j*re,k*re);
shading interp;
earth.FaceColor = 'b';
earth.EdgeColor = 'b';
earth.FaceAlpha = 0.1;
earth.EdgeAlpha = 0.1;



% Add axis for the IJK vectors
line([-re,re],[0,0],[0,0], 'LineWidth',1, 'Color', 'k');
line([0,0],[-re,re],[0,0], 'LineWidth',1, 'Color', 'k');
line([0,0],[0,0],[-re,re], 'LineWidth',1, 'Color', 'k');
text(0,0,re*1.1,'$\left| z \right>$', 'Interpreter', 'latex', ...
    'FontSize', 20, 'HorizontalAlignment', 'Center');
text(0,re*1.1,0,'$\left| y \right>$', 'Interpreter', 'latex', ...
    'FontSize', 20, 'HorizontalAlignment', 'Center');
text(re*1.1,0,0,'$\left| x \right>$', 'Interpreter', 'latex', ...
    'FontSize', 20, 'HorizontalAlignment', 'Center');

legend('s1','MdInit','MdCurrentLoc');

hold off;


%% Map the Ground track
figure(4)
geoshow('landareas.shp','FaceColor',[0.5 1 0.5]);
title('Satellite Ground Track')
hold on
plot(long*(180/pi),phigd*(180/pi),'.r')
hold off



%% Animate the Simulation if the Satellite is visible

% 
% for i =1:length(t)
%     
%     if LOS(i) ==1
%         
%         loc = i;
%         
%         p0 = [0 0 0];
%         
%         define the line from the center of the ECI frame to the university of
%         marylalnd at that particular timestep
%         pmd = [mdx(loc), mdy(loc), mdz(loc)];
%         mdline = [p0;pmd];
%         
%         define the line from the university of maryland to the satellite
%         at that particular timestep
%         psat = [satx(loc), saty(loc), satz(loc)];
%         satline = [pmd;psat];
%         
%         
%         define the line from the center of the ECI frame to the university of
%         marylalnd at that particular timestep
%         perp= [perpx(loc), perpy(loc), perpz(loc)];
%         perpline = [pmd;perp];
%         
%         
%         Add the NWT frame at the location of the satellite at that particular
%         timestep
%         nwtmag = 1000;
%         n=[Nhat(loc,1),Nhat(loc,2),Nhat(loc,3)];
%         nline = [psat;psat+nwtmag*n];
%         
%         W=[What(loc,1),What(loc,2),What(loc,3)];
%         Wline = [psat;psat+nwtmag*W];
%         
%         T=[That(loc,1),That(loc,2),That(loc,3)];
%         Tline = [psat;psat+nwtmag*T];
%         
%         currently the perp line is correctly in the same direction as the line to
%         the satellite but it is not perpendicular to the line to maryland.
%         figure(3)
%         
%         plot3(rx(:,1),ry(:,1),rz(:,1), 'k', xmdinit, ymdinit, zmdinit, '*', mdx(:,1), mdy(:,1),mdz(:,1), 'y',...
%             mdline(:,1), mdline(:,2), mdline(:,3), 'r', satline(:,1), satline(:,2), satline(:,3), 'c',...
%             nline(:,1), nline(:,2), nline(:,3), 'r',Wline(:,1), Wline(:,2), Wline(:,3), 'g',Tline(:,1), Tline(:,2), Tline(:,3), 'b',...
%             satx(loc), saty(loc), satz(loc), 'o');
%         axis equal;
%         grid on;
%         drawnow;
%         
%         
%         xlabel('x (km)')
%         ylabel('y (km)')
%         zlabel('z (km)')
%         pause(0.01);
%     end
%        
%     
% end


% 
% states = shaperead('usastatehi', 'UseGeoCoords', true);
% 
% symspec = makesymbolspec('Polygon', ...
%    {'Name', 'Alaska', 'FaceColor', 'red'}, ...
%    {'Name', 'Hawaii', 'FaceColor', 'red'});
% 
% figure(4)
% title('Satellite Ground Track')
% geoshow(states, 'SymbolSpec', symspec, ...
%    'DefaultFaceColor', 'blue', ...
%    'DefaultEdgeColor', 'black');
% hold on
% plot(long*(180/pi),phigd*(180/pi),'.r');
% 
% 
% 
% hold off
% 
% an = animatedline('Marker','*');
% for k = 1:length(t)
%     addpoints(an,long(k)*(180/pi),phigd(k)*(180/pi));
%     drawnow
%     pause(0.2);
%     clearpoints(an);
% end




%% Animate Ground Track movement
% an = animatedline('Marker','*');
% for k = 1:length(t)
%     addpoints(an,long(k)*(180/pi),phigd(k)*(180/pi));
%     drawnow
%     pause(0.1);
%     clearpoints(an);
% end



%% Map Ground track only when there is Line of Sight
% 
% 
% figure(5)
% 
% creates an array that contains longitudes, only if the satellite is within
% line of sight. 
% longLOS = long.*LOS;
% 
% geoshow('landareas.shp','FaceColor',[0.5 1 0.5]);
% title('Satellite Ground Track')
% hold on
% an = animatedline('Marker','*');
% for i = 1:length(t)
%     if longLOS(i) ~= 0
%         plot(longLOS(i)*(180/pi),phigd(i)*(180/pi),'.r')
%         
%     end
% end
% hold off
% 
% for i = 1:length(t)
%     if longLOS(i) ~= 0
%         addpoints(an,longLOS(i)*(180/pi),phigd(i)*(180/pi));
%         drawnow
%         pause(0.1);
%         clearpoints(an);
%     end
% end


