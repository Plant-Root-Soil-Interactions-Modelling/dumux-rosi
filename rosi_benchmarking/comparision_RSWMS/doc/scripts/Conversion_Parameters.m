% Conversion_Parameters.m
% converts dumux to rswms parameters
% writes rswms parameters in another file

%%
close all; clear all;
%cd 'C:\Users\k.huber\Dumux-rswms\';
fileDUMUX = 'inputDUMUX.txt';
fileRSWMS = 'inputRSWMS.txt';

%% conversion factors
dummy = dlmread(fileDUMUX,'',[2 0 2 2]);
rho = dummy(1); % density of water [kg/m^3]
mu = dummy(2); % dynamic viscosity, water 20� [Pa s]
g = dummy(3); % gravitational acceleration [m/s^2]
clear dummy;

%% DUMUX input; SI units
% time
dummy = dlmread(fileDUMUX,'',[6 0 6 2]);
TSd = dummy(1); % time step, [s]
TEd = dummy(2); % time end, [s]
TEpd = dummy(3); % time end, [s]
clear dummy;

% soil grid
LLd = dlmread(fileDUMUX,'',[12 0 12 2]); % lower left [m]
URd = dlmread(fileDUMUX,'',[14 0 14 2]); % upper right [m]
Cd = dlmread(fileDUMUX,'',[16 0 16 2]); % number of cells

% Van Genuchten parameters
dummy = dlmread(fileDUMUX,'',[20 0 20 2]);
SWRd = dummy(1);
ALPHAd = dummy(2); % [1/m]
Nd = dummy(3);
clear dummy;

dummy = dlmread(fileDUMUX,'',[23 0 23 1]);
PERMd = dummy(1); % permeabilty [m^2]
POROd = dummy(2); % porosity
clear dummy;

% Root hydraulic parameters
dummy = dlmread(fileDUMUX,'',[27 0 27 1]);
KXd = dummy(1); % [m^4/s/Pa]
KRd = dummy(2); % [m/s/Pa]
clear dummy;

% Boundary conditions
dummy = dlmread(fileDUMUX,'',[31 0 31 2]);
PHSINId = dummy(1); %[Pa]
PHRINId = dummy(2); % [Pa]
TPOTd = dummy(3); % [kg/s]
clear dummy;

%% rswms input, cm, day
% time
TSr = TSd/86400; % time step [d]
TEr = TEd/86400; % time end [d]
TEpr = TEpd/86400; % time end [d]

% soil grid
LLr = LLd*100; % lower left [cm]
URr = URd*100; % upper right [cm]
Cr = Cd; % number of cells (elements)
dx = (URr(1)+abs(LLr(1)))/Cr(1);
dy = (URr(2)+abs(LLr(2)))/Cr(2);
dz = (URr(3)+abs(LLr(3)))/Cr(3);

% Van Genuchten parameters
SWRr = SWRd; % residual soil water content
ALPHAr = ALPHAd*100; % [1/cm]
Nr = Nd;

% saturated hydraulic conductivity [cm/d], with the assumption that the 
% relative permeability kappa at saturation is 1
KSATr = PERMd*1*rho*g/mu*86400*100;
SWSr = POROd; % saturated soil water content

% Root hydraulic parameters
KXr = KXd*86400*100^4*rho*g/100; % axial hydraulic conductivity[cm^4/day/cm]
KRr = KRd*86400*rho*g; % radial hydraulic conductivity [cm/d/cm]


% Boundary conditions
PHSINIr = PHSINId/rho/g*100; % [cm]
htop = PHSINIr; % initial PH at the soil top [cm]
PHRINIr = PHRINId/rho/g*100; % [cm]
TPOTr = TPOTd*86400/rho*100^3; % [cm^3/d]

%% write rswms inputs
fid=fopen(fileRSWMS,'w'); 

fprintf(fid,'%s\n','RSWMS inputs');
fprintf(fid,'%s\n','TIME - control.in');
fprintf(fid,'%s\n','dt[d]     end time[d]');
fprintf(fid,'%4d\t %4d\t %4d\n\n',TSr,TEr,TEpr);

fprintf(fid,'%s\n','SOIL GRID - mesh.in [cm]');
fprintf(fid,'%s\n','dx   dy   dz   nex   ney   nez   xmin   ymin   zmax');
fprintf(fid,'%5d\t %5d\t %5d\t %5d\t %5d\t %5d\t  %5d\t %5d\t %5d\n\n',...
    dx,dy,dz,Cr(1),Cr(2),Cr(3),LLr(1),LLr(2),URr(3));
fprintf(fid,'%s\n','htop[cm]');
fprintf(fid,'%4.5d\n\n',PHSINIr);

fprintf(fid,'%s\n','SOIL HYDRAULIC PROPERTIES - soil.in');
fprintf(fid,'%s\n','thr[-]  ths[-]  a[1/cm]    n     Ks[L/T]');
fprintf(fid,'%4.5d\t %4.5d\t %4.5d\t %4.5d\t %4.5d\n\n',...
    SWRr,SWSr,ALPHAr,Nr,KSATr);

fprintf(fid,'%s\n','ROOT HYDRAULIC PROPERTIES - CondRoot.in');
fprintf(fid,'%s\n','Lr[cm/d]  KhRoot[cm^3/d]');
fprintf(fid,'%4.5d\t %4.5d\n\n',KRr,KXr);

fprintf(fid,'%s\n','ROOT BOUNDARY CONDITION - BCroot.in');
fprintf(fid,'%s\n','Tpot[cm^3/d]');
fprintf(fid,'%4.5d\n\n',TPOTr);

fclose(fid);






