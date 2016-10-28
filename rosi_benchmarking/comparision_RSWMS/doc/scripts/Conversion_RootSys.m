% Conversion_RootSys.m
% reads dumux .dgf file
% writes RootSys file

%%
close all; clear all;
cd 'C:\Users\k.huber\Dumux-rswms\';
fileDUMUX = 'RootSysMRI_1times.dgf';
fileRSWMS = 'RootSysMRI_rswms';
n_seg = 1739; % number of root segments, !Dumux starts with node 0!

%% read data
% coordinates [m]
Coordinates = read_dumuxroot_coord(fileDUMUX,3,n_seg+2);

% node#, next node, root order, branch#, surface[m^2], mass[kg]???
% -->volume[m^3]
Parameters = read_dumuxroot_param(fileDUMUX,n_seg+6,n_seg*2+4);

%% Conversion to rswms
% convert next node to previous node
% Parameters is already sorted by next node
% rename to simplify writing
% ommit first node (ID=0), also has Z>0 in dumux
ID = Parameters(:,2);
prev = Parameters(:,1);
order = Parameters(:,3);
brID = Parameters(:,4);
segsur = Parameters(:,5)*100^2; %cylindrical surface [cm^2]
volume = Parameters(:,6)/100^3; %[cm^3]
Xr = Coordinates(2:n_seg,1)*100; % [cm]
Yr = Coordinates(2:n_seg,2)*100; % [cm]
Zr = Coordinates(2:n_seg,3)*100; % [cm]
n_seg = n_seg-1;

%
n_br = max(Parameters(:,4));

% Dumux starts with Node 0!
% NodeR = Parameters(:,1)+1;
% NxtR = Parameters(:,2)+1;
% RSWMS calculates length between segment and PREVIOUS segment
seglen = zeros(n_seg,1);
radius = zeros(n_seg,1);

for ii=1:n_seg
    seg = Parameters(ii,2); %rswms = seg, dumux = next
    prv = Parameters(ii,1); %prev seg
    
    if (prv == 0)   % connected to the collar
        Xr0 = Coordinates(1,1)*100;
        Yr0 = Coordinates(1,2)*100;
        Zr0 = Coordinates(1,3)*100;
        seglen(seg) = sqrt((Xr(seg)-Xr0)^2 +...
            (Yr(seg)-Yr0)^2 + (Zr(seg)-Zr0)^2);
        
    else        
        seglen(seg) = sqrt((Xr(seg)-Xr(prv))^2+...
            (Yr(seg)-Yr(prv))^2+(Zr(seg)-Zr(prv))^2);
        radius(seg) = segsur(seg)/2/pi/seglen(seg);
            
    end
end



% Tip information
brlen = zeros(n_br,1);
sgbhtip = zeros(n_br,1);
brorder = zeros(n_br,1);
brnumb = zeros(n_br,1);
Xg = zeros(n_br,1);
Yg = zeros(n_br,1);
Zg = zeros(n_br,1);

for ibr=1:n_br
    OK_br=find(brID==ibr); % segment IDs from branches in new RootSys
    
    tipID = OK_br(end);
    Xg(ibr) = Xr(tipID);
    Yg(ibr) = Yr(tipID);
    Zg(ibr) = Zr(tipID);
    sgbhtip(ibr) = tipID; %information of tip == information of last segment in branch
    
    brorder(ibr) = order(tipID);
    brnumb(ibr) = ibr;
    brlen(ibr) = sum(seglen(OK_br));
    clear OK_br tipID;
end


%% WRITE ROOTSYS
fid=fopen(fileRSWMS,'w');

fprintf(fid,'%s\n','Time: ');
fprintf(fid,'%4.5d\n\n',0);

fprintf(fid,'%s\n','Number of seeds ');
fprintf(fid,'%s\n\n','   1');

fprintf(fid,'%s\n','ID, X and Y coordinates of the seeds (one per line)');
fprintf(fid,'%s\n\n','   1       0.0000000       0.0000000    ');

fprintf(fid,'%s\n','IRoot DM, shoot DM, leaf area');
fprintf(fid,'%s\n\n','   0.0000000       0.0000000       0.0000000  ');

fprintf(fid,'%s\n','Average soil strength and solute concentration experienced by root system:');
fprintf(fid,'%s\n\n',' 0.0000000       0.0000000 ');

fprintf(fid,'%s\n','Total # of axes: ');
fprintf(fid,'%s\n\n','   1');

fprintf(fid,'%s\n','Total # of branches, including axis(es)');
fprintf(fid,'\t%d\n\n',n_br);

fprintf(fid,'%s\n','Total # of segment records:');
fprintf(fid,'\t%d\n\n',n_seg);

fprintf(fid,'%s\n','segID#        x               y              z            prev  or   br#   length    surface      volume');
fprintf(fid,'%s\n','origination time');

for i=1:n_seg;
    fprintf(fid,'  %5d\t %+6.4E\t %+6.4E\t %+6.4E\t %4d\t %d\t %d\t %6.4E %6.4E %6.4E\n',...
        ID(i), Xr(i),Yr(i), Zr(i),  prev(i), order(i), brID(i), seglen(i),segsur(i), volume(i));
    fprintf(fid,'%6.4E\n',0);
end

fprintf(fid,'\n%s\n','Total # of growing branch tips:') ;
fprintf(fid,'%d\n\n',n_br);

fprintf(fid,'%s\n','tipID#    xg          yg          zg      sg.bhd.tp. ord  br#  tot.br.lgth. axs#');
fprintf(fid,'%s\n','overlength  # of estblished points');
fprintf(fid,'%s\n','time of establishing (-->)');

for i=1:n_br
    fprintf(fid,'  %5d\t %+6.4E\t %+6.4E\t %+6.4E\t %4d\t %d\t %d\t %6.4E \t %d\n',...
        i, Xg(i), Yg(i), Zg(i), sgbhtip(i), brorder(i), brnumb(i),brlen(i), 1);
    fprintf(fid,'%6.4E %d\n',1,0);
end

fclose(fid);



