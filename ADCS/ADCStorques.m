% Author: Brady Beck
% AAE 450   
% Script Description: Used to solve for ADCS related torques and forces.
clear
clear all 

% Definitions from SolidWorks - OLD, DID NOT WORK
Rcom = [6.10; 2.50; 4.54]; % [m] Center of mass 
Rcross = vectorCross(Rcom); 
Jo = [63208.15 60307.92 33064.18; 
     60307.92 97535.03 24516.59; 
     33064.18 24516.59 130552.40]; % [kg*m^2] inertia tensor 
mDry = 1527.33; % [kg] estimated orbiter dry mass 
mfrac = mDry * 5574/2523; % [kg] estimated mass fraction from cassini scaling
Jsw = Jo - mfrac*(Rcross*Rcross');

% Define J based on dry mass fraction and cassini inertia tensor
J = 1527.33/2523 * diag([4711.1; 8138.9; 8.8399]); 
m = 1527.33*5574/2523; % [kg] estimated total mass based on cassini scaling

% Gravitational Torque 
% Setting constants of visited bodies
muEarth = 3.986e5; % [km^3/s^2] Earth gravitational parameter
muJupiter = 1.266e8; % [km^3/s^2] Jupiter gravitational parameter
muUranus = 5.793e6; % [km^3/s^2] Uranus gravitational parameter
G = 6.674e-20; % [km^3/kg*s^2] gravitational constant
muTitania = 3.4e21*G; % [km^3/s^2] Titania gravitational parameter
muAriel = 1.251e21*G; % [km^3/s^2] Ariel gravitational parameter
muOberon = 3.076e21*G; % [km^3/s^2] Oberon gravitational parameter
muUmbriel = 1.275e21*G; % [km^3/s^2] Umbriel gravitational parameter

% Setting radius values for closest vist to each body 
rEarth = 6578.14; % [km] Earth orbit radius 
rJupiter = 27.388*69911; % [km] Jupiter orbit radius 
rTitania = 2.214*788.4; % [km] Titania orbit radius 
rAriel = 2.914*578.9; % [km] Ariel orbit radius
rOberon = 2.708*761.4; % [km] Oberon orbit radius
rUmbriel = 1.1*584.7; % [km] Umbriel orbit radius
rUranus = 7.510*25362; % [km] Uranus orbit radius

% Calculating Lgg
n = [sqrt(3)/3; sqrt(3)/3; sqrt(3)/3]; % n vector to maximize Lgg
LggEarth = (3*muEarth/rEarth^3)*cross(n,J*n);
LggJupiter = (3*muJupiter/rJupiter^3)*cross(n,J*n);
LggUranus = (3*muUranus/rUranus^3)*cross(n,J*n);
LggTitania = (3*muTitania/rTitania^3)*cross(n,J*n);
LggAriel = (3*muAriel/rAriel^3)*cross(n,J*n);
LggOberon = (3*muOberon/rOberon^3)*cross(n,J*n);
LggUmbriel = (3*muUmbriel/rUmbriel^3)*cross(n,J*n);

% Printing Results
fprintf('\n-=- Gravity Gradient Torque Results -=-\n');
fprintf('LggEarth: %d [Nm]\n', norm(LggEarth));
fprintf('LggJupiter: %d [Nm]\n', norm(LggJupiter));
fprintf('LggUranus: %d [Nm]\n', norm(LggUranus));
fprintf('LggTitania: %d [Nm]\n', norm(LggTitania));
fprintf('LggArien: %d [Nm]\n', norm(LggAriel));
fprintf('LggOberon: %d [Nm]\n', norm(LggOberon));
fprintf('LggUmbriel: %d [Nm]\n\n', norm(LggUmbriel));
fprintf('Max Lgg: Earth, %d [Nm]\n\n', norm(LggEarth));
LggMAX = norm(LggEarth);

% Magnetic Torque 
% Note: using equation from SMAD 573, Tm = D(M*lambda/R^3)
D = m*4e-3; % [A*m^2] magnetic dipole moment, estimated with NASA paper and estimated mass
lambdaMax = 2; % unitless, from SMAD 573
Mearth = 7.8e15; % [Tesla*m^3] Earth magnetic moment, from SMAD 573
Mjupiter = 2.83e20; % [Tesla*m^3] Jupiter magnetic moment
LmagEarth = D*Mearth*lambdaMax/(rEarth*1000)^3;
LmagJupiter = D*Mjupiter*lambdaMax/(rJupiter*1000)^3;

% Print results
% fprintf('\n-=- Magnetic Moment Torque Results -=-\n');
% fprintf('LmagEarth: %d [Nm]\n', LmagEarth);
% fprintf('LmagJupiter: %d [Nm]\n', LmagJupiter);
% fprintf('LmagUranus: Not considered due to uncahracterized field\n\n');
% fprintf('Max Lmag: Earth, %d [Nm]\n\n', LmagEarth);
LmagMAX = LmagEarth;

% Atmospheric Torque 
r = [0; 0; 3]; % [m] estimated worst case CoM to Center Solar Radiation Pressure vector
A = 24; % [m^2] estimation for wetted area
rho = 4.1e-13; %[kg/m^3]
V = 7.79*1000; % [m/s] orbital veloctiy at 200km 
Cd = 2.5; % estimation from SMAD pg 573
Fd = 0.5*rho*V^2*Cd*A.*[0; 1; 0]; % vector to maximize cross product
Ldrag = norm(cross(r,Fd));
% fprintf('\n-=- Atmospheric Torque Results -=-\n');
% fprintf('LsrpEarth: %d [Nm]\n', Ldrag);
% fprintf('Other Atmospheric not considered\n because will be maximum condition at Earth\n\n');
% fprintf('Max Lsrp: Earth, %d [Nm]\n\n', Ldrag);
LdragMAX = Ldrag;

% Solar Radiation Torque 
r = [0; 0; 3]; % [m] estimated worst case CoM to Center Solar Radiation Pressure vector
phi = 1367; % [W/m^2] given from lecture 
c = 2.998e8; % [m/s] speed of light
q = 1; % unitless, reflectance factor, worst case 
A = 24; % [m^2] estimation for wetted area 
i = 0; % [rad] incidence angle of sun 
Fsrp = phi/c*A*(1+q)*cos(i).*[sqrt(2)/2; sqrt(2/2); 0]; % vector to maximize cross product
Lsrp = cross(r, Fsrp);
% fprintf('\n-=- Solar Radiation Pressure Torque Results -=-\n');
% fprintf('LsrpEarth: %d [Nm]\n', norm(Lsrp));
% fprintf('Other Solar Radiation Pressure Torques not considered\n because will be maximum when at Earth\n\n');
% fprintf('Max Lsrp: Earth, %d [Nm]\n\n', norm(Lsrp));
LsrpMAX = norm(Lsrp);

% Summing total environmental torque
Lenv = LggMAX + LmagMAX + LdragMAX + LsrpMAX;
fprintf('\n-=- Worst Case Scenario Environmental Torque -=-\n');
fprintf('Lenv: %d [Nm]\n', Lenv);

% Reaction Wheel Torque
% I = J(2,2);
I = 27430;
Lw = 4*pi*I/(3600)^2;
fprintf('\n-=- Required Reaction Wheel Torque for Slew Requirement -=-\n');
fprintf('Lw: %d [Nm]\n', Lw);

% Required total control torque
marg = 1.3; % added margin of error due to pre-phase A
Lreq = marg*(Lenv + Lw);
fprintf('\n-=- Total Required Control Torque -=-\n');
fprintf('Lreq: %d [Nm]\n', Lreq);

% Force of Thrusters Calculated 
Lvec = [0; 0; .1];
vRelPY = [0; 1/sqrt(2); 1/sqrt(2)];
vRelR = [1; 0; 0];
thrusterR = [4; -4; -6] - [-4; 4; -6]; % [m] distance between thrusters
Fexp = norm(Lvec)/norm(cross(thrusterR,vRelPY));
fprintf('\n-=- Required Force for Equivilant Thruster -=-\n');
fprintf('Fexp: %d [N]\n', Fexp);

