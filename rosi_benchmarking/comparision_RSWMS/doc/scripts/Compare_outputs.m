% Compare_outputs.m
% compares dumux and rswms outputs

%% USER INPUT
close all; clear all;
folderML = '/home/h.mai/PROGRAM/rosi/rosi_master1/dumux-rosi/rosi_benchmark/vs_RSWMS/tools';
folderPLOT = '/home/h.mai/PROGRAM/rosi/rosi_master1/dumux-rosi/rosi_benchmark/vs_RSWMS/results_comparision';
folderDUMUX = '/home/h.mai/PROGRAM/rosi/rosi_master1/dumux-rosi/build-cmake/rosi_benchmark/vs_RSWMS/richardsrootsystem';
folderRSWMS = '/home/h.mai/PROGRAM/rosi/rosi_master1/dumux-rosi/rosi_benchmark/vs_RSWMS/RSWMS7/out20';
addpath(folderML,folderDUMUX,folderRSWMS);

% file prefix
soilDUMUX = 'rosi20-soil-';
soilRSWMS = 'veloci.';
rootDUMUX = 'rosi20-root-';
rootRSWMS = 'outRoot01.';

% number of outputs
nout = 11;
simtime = linspace(0,86400,11)/86400; % [d]

% grid information
% number of elements (voxels)
nex = 40; ney = 40; nez = 40;

xmin = -5; xmax = 5; % [cm]
ymin = -5; ymax = 5; % [cm]
zmin = -10; zmax = 0; % [cm]

% number of root segments (w/o tips)
nseg = 1739;
nsegR = 1796; % contains also tips!

% density of water
rho = 1000; % [kg/m^3]

%% soil grid
npx = nex+1; npy = ney+1; npz = nez+1;
nelm = nex*ney*nez;
nnodes = npx*npy*npz;

% node coordinates
Xg = linspace(xmin,xmax,npx);
Yg = linspace(ymin,ymax,npy);
Zg = linspace(zmin,zmax,npz);

% coordinates of element centres
dx = Xg(2)-Xg(1); dy = Yg(2)-Yg(1); dz = Zg(2)-Zg(1);
Xe = linspace(xmin+dx/2,xmax-dx/2,nex);
Ye = linspace(ymin+dy/2,ymax-dy/2,ney);
Ze = linspace(zmin+dz/2,zmax-dz/2,nez);

%% initiate matrices
PHsD_mean = zeros(nez,nout);
WCD_mean = zeros(nez,nout);
sourceD_sum = zeros(nez,nout);

PHsR_mean = zeros(npz,nout);
WCR_mean = zeros(npz,nout);
sinkR_sum = zeros(nez,nout);

PXD = zeros(nseg,nout);
PXR = zeros(nseg,nout);

%% SOIL
% DUMUX
% dumux coordiate system starts bottom left front
cd (folderDUMUX);
startPH = 7+11*(2+ceil(nelm/12)); % 11th data set
startWC = 7+12*(2+ceil(nelm/12)); 
startsource = 1+17*(2+ceil(nelm/12));
for ii=1:nout
    aa = sprintf('%05d',ii-1);
    filein = [soilDUMUX aa '.vtu'];
    
    PHs_dumux = read_dumuxsoil_single(filein,startPH, startPH+ceil(nelm/12)-1); %pressure head
%     PHs_dumux = read_dumuxsoil_single(filein,12044, 12710); %
%     matricPotential == pressure head
    WC_dumux = read_dumuxsoil_single(filein,startWC, startWC+ceil(nelm/12)-1);
    source_dumux = read_dumuxsoil_single(filein,startsource, startsource+ceil(nelm/12)-1); % sink [kg/s]
    
    % restructure matrices --> 1d
    for jj=1:length(PHs_dumux(:,2)) % rows in dumux output
        val1 = (jj-1)*12+1; % 12 columns in dumux output
        val2 = val1+11;
        PHsD(val1:val2) = PHs_dumux(jj,:);
        WCD(val1:val2) = WC_dumux(jj,:);
        sourceD(val1:val2) = -source_dumux(jj,:)*86400/rho*100^3; % [cm^3/d]
    end
    
    
    %% average variables over z
    for jj=1:nez
        val1 = (jj-1)*nex*ney+1;
        val2 = val1+nex*ney-1;
        PHsD_mean(jj,ii) = mean(PHsD(val1:val2));
        WCD_mean(jj,ii) = mean(WCD(val1:val2));
        sourceD_sum(jj,ii) = sum(sourceD(val1:val2));
    end
    
    clear PHs_dumux WC_dumux source_dumux PHsD WCD sourceD;
end


% RSMWS
cd (folderRSWMS);
startPH = 17+nnodes; %
startWC = 18+2*nnodes;
startS = 25+7*nnodes; 
for ii=1:nout
    filein = [soilRSWMS num2str(ii) '.vtk'];
    
    %     PHs, WC are node based, sinkElm is element based
    PHs_rsmws = read_rswmssoil(filein, startPH, startPH+nnodes-1);
    WC_rswms = read_rswmssoil(filein,startWC, startWC+nnodes-1);
    sink_rswms = read_rswmssoil(filein,startS, startS+nelm-1);
    
    %% average variables over z
    %     node based
    for jj=1:npz
        val1 = (jj-1)*npx*npy+1;
        val2 = val1+npx*npy-1;
        PHsR_mean(jj,ii) = mean(PHs_rsmws(val1:val2));
        WCR_mean(jj,ii) = mean(WC_rswms(val1:val2));
    end
    
    %     element based
    for jj=1:nez
        val1 = (jj-1)*nex*ney+1;
        val2 = val1+nex*ney-1;
        sinkR_sum(jj,ii) = sum(sink_rswms(val1:val2))*dx*dy*dz;
    end
    clear PHs_rswms WC_rswms sink_rswms;
end


%% ROOT
% Dumux
% PHX in dumux is a 'line' value (nseg-1), not a point value
cd (folderDUMUX);
startX = 7+15*2+14*ceil(nseg/12)+ceil(nseg*3/12); 
for ii=1:nout
    aa = sprintf('%05d',ii-1);
    filein = [rootDUMUX aa '.vtp'];
    PX_dumux = read_dumuxroot_out(filein,startX, startX+ceil(nseg/12)-1);
    
    % restructure matrices --> 1d
    for jj=1:length(PX_dumux(:,2)) % rows in dumux output
        val1 = (jj-1)*12+1; % 12 columns in dumux output
        val2 = val1+11;
        PXD(val1:val2,ii) = PX_dumux(jj,:);
    end
    
    clear PX_dumux;
end
PXD = PXD(1:nseg-1,:); % remove NAN values

% RSWMS
cd (folderRSWMS);
startX = 13+5*nsegR-1; 
for ii=1:nout
    filein = [rootRSWMS num2str(ii) '.vtk'];
    
    PXR(:,ii) = read_rswmsroot_out(filein,startX, startX+nseg-1);
end

%% PLOT
% !first value for dumux == bottom!
% !first value for rswms == top!
% !first value for grid == bottom!
cd (folderPLOT);
col_aut = colormap(autumn(nout)); col_sum = colormap(winter(nout));
close all;

% soil pressure head
figure; hold on; 
for ii=1:nout
    curveD(ii) = plot(PHsD_mean(:,ii),Ze,'LineWidth',1.8,'Color', col_aut(ii,:));
    curveR(ii) = plot(PHsR_mean(:,ii),linspace(zmax,zmin,npz),'LineWidth',1.8,...
        'LineStyle','--','Color', col_sum(ii,:));
end

hXLabel = xlabel('Soil water potential [cm]', 'FontSize', 14);
hYLabel = ylabel('Soil depth [cm]', 'FontSize', 14);

hLegend = legend( ...
    [ curveD(1),  curveR(1)], ...
    'Dumux',...
    'RSWMS',...
    'location', 'NorthEast');
set( [ hLegend, gca ], 'FontSize', 12,'Box','off');

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'LineWidth'   , 1         );

set(gcf, 'PaperPositionMode', 'auto');
saveas( gcf, 'PHs.png');


% soil water content
figure; hold on;
for ii=1:nout
    curveD(ii) = plot(WCD_mean(:,ii),Ze,'LineWidth',1.8,'Color', col_aut(ii,:));
    curveR(ii) = plot(WCR_mean(:,ii),linspace(zmax,zmin,npz),'LineWidth',1.8,...
        'LineStyle','--','Color', col_sum(ii,:));
end

hXLabel = xlabel('Soil water content [-]', 'FontSize', 14);
hYLabel = ylabel('Soil depth [cm]', 'FontSize', 14);

hLegend = legend( ...
    [ curveD(1),  curveR(1)], ...
    'Dumux',...
    'RSWMS',...
    'location', 'NorthEast');
set( [ hLegend, gca ], 'FontSize', 12,'Box','off');

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'LineWidth'   , 1         );

set(gcf, 'PaperPositionMode', 'auto');
saveas( gcf, 'WC.png');

% sink
figure; hold on;
for ii=1:nout
    curveD(ii) = plot(sourceD_sum(:,ii),Ze,'LineWidth',1.8,'Color', col_aut(ii,:));
    curveR(ii) = plot(sinkR_sum(:,ii),linspace(zmax-dz/2,zmin+dz/2,nez),'LineWidth',1.8,...
        'LineStyle','--','Color', col_sum(ii,:));
end

hXLabel = xlabel('Sink [cm^3 d^-^1]', 'FontSize', 14);
hYLabel = ylabel('Soil depth [cm]', 'FontSize', 14);

hLegend = legend( ...
    [ curveD(1),  curveR(1)], ...
    'Dumux',...
    'RSWMS',...
    'location', 'NorthEast');
set( [ hLegend, gca ], 'FontSize', 12,'Box','off');

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'LineWidth'   , 1         );

set(gcf, 'PaperPositionMode', 'auto');
saveas( gcf, 'sink.png');


% xylem pressure head in root collar
figure; hold on;

curveD = plot(simtime(2:nout),PXD(1,2:nout),'LineWidth',1.8,'Color', 'r');
curveR = plot(simtime(2:nout),PXR(1,2:nout),'LineWidth',1.8,...
    'LineStyle','--','Color', 'k');

hXLabel = xlabel('Time [d]', 'FontSize', 14);
hYLabel = ylabel('PHX_{collar} [cm]', 'FontSize', 14);

hLegend = legend( ...
    [ curveD,  curveR], ...
    'Dumux',...
    'RSWMS',...
    'location', 'NorthEast');
set( [ hLegend, gca ], 'FontSize', 12,'Box','off');

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'LineWidth'   , 1         );

set(gcf, 'PaperPositionMode', 'auto');
saveas( gcf, 'PHX.png');



