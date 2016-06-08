%% Load Nortek data files and store according to waves and currents data
% This function uses the both ways that the files are converted by the
%   program used to extract the data directly from the instruments. Vary
%   according to the update of the sofware.
% The program converts the Aquadopp data (the same as AWAC but does not
%   include STrk data).

%% Load different files
% a, b, c, d , and e correspond to BOEM folder and data naming scheme.
% Refer to 'main_data' script to get these for a specific deployment.
% To load the data, path in Matlab must be within 'Functions'. This script
% will go up one folder and then it will look for specific folders
% according to the BOEM scheme.

% Path and filename
% % % [pathname,filename]=filepath_data(a,b,c,d,e);
% function [curr,wave]=adcp01_LoadRaw(pathname,filename,coord);
pathname='/Users/pna/Dropbox/Data/Data02_ADCP/20090212_MtnzsDplmnt/AqDp/RawData/';
filename='mtz09a01';
coord=[29 42.840; -81 -13.298];
clc, close all
fprintf('======================================\n')
fprintf('LOAD RAW DATA FROM NORTEK AQUADOPP\n')
fprintf('Geomorphology Laboratory\n')
fprintf('Department of Geological Sciences\n')
fprintf('University of Florida\n')
fprintf('Gainesville, FL, USA\n')
fprintf('Summer of 2016\n')
fprintf('======================================\n')

tic;

fprintf('1. LOAD TEXT FILES\n')

% Load files for Nortek formatting
sen=load([pathname filename '.sen']);
fprintf('.sen file loaded...\n')

v1=load([pathname filename '.v1']);
v2=load([pathname filename '.v2']);
v3=load([pathname filename '.v3']);
fprintf('.v1, .v2, and .v3 files loaded...\n')

a1=load([pathname filename '.a1']);
a2=load([pathname filename '.a2']);
a3=load([pathname filename '.a3']);
fprintf('.a1, .a2, and .a3 files loaded...\n')

wad=load([pathname filename '.wad']);
fprintf('.wad file loaded...\n')

whd=load([pathname filename '.whd']);
fprintf('.whd file loaded...\n')


% Deployment information
% Retrieve data from .hdr file
fprintf('\nDeployment information from .hdr file\n')
fprintf('----------------------------------------------\n')

fid = fopen([pathname filename '.hdr']);

for k=1:23
    if k==11 % Number of cells
        tline=fgets(fid);
        disp(tline)
        Ncells=cell2mat(textscan(tline,'%*s %*s %*s %f %*s'));
        
    elseif k==12 % Cell size, m
        tline=fgets(fid);
        disp(tline)
        cellsize=cell2mat(textscan(tline,'%*s %*s %f %*s'))/100;
        
    elseif k==16 % Blanking distance, m
        tline=fgets(fid);
        disp(tline)
        blankdist=cell2mat(textscan(tline,'%*s %*s %f %*s'));
        
    elseif k==15 % Transmit pulse length, m
        tline=fgets(fid);
        disp(tline)
        Lp=cell2mat(textscan(tline,'%*s %*s %*s %f %*s'));
        
    elseif k==23 % Frequency of data collection (burst mode)
        tline=fgets(fid);
        disp(tline)
        freq=cell2mat(textscan(tline,'%*s %*s %*s %*s %f %*s'));
        
    elseif k==22 % Number of wave samples
        tline=fgets(fid);
        disp(tline)
        Nw=cell2mat(textscan(tline,'%*s %*s %*s %*s %*s %f %*s'));
        
    else
        % Read current line within text
        tline=fgets(fid);
        
    end
end

fclose(fid);

% Vector with the position of each cell center measured from the instrument
% Remember, first cell center is located at 'blankdist + 1 cellsize'
cells=(blankdist+(1:Ncells)*cellsize)';

fprintf('--------------------------------------------------------------\n')

fprintf('Deployment information retrieved...\n')

%% Currents data
fprintf('2. CURRENTS DATA\n')

% z_p: upward distance from the instrument) which are averaged to obtain
%   velocity values.
curr.z_p=cells;

% Velocities
curr.u=v1;
curr.v=v2;
curr.w=v3;

% Echo intensity
curr.RSSI1=a1;
curr.RSSI2=a2;
curr.RSSI3=a3;

% Time
yearsc=sen(:,3);
monthsc=sen(:,1);
daysc=sen(:,2);
hoursc=sen(:,4);
minsc=sen(:,5);
secsc=sen(:,6);
curr.t=datenum(yearsc,monthsc,daysc,hoursc,minsc,secsc);
curr.tJ=curr.t-datenum(yearsc(1),0,0);

% Pressure
curr.P=sen(:,14);

% Temperature
curr.T=sen(:,15);

% Heading
curr.heading=sen(:,11);

% Pitch
curr.pitch=sen(:,12);

% Roll
curr.roll=sen(:,13);

% Voltage
curr.volt=sen(:,9);

% Message
fprintf('Currents data retrieved...\n')

%% Wave data
% Message
fprintf('3. WAVES DATA\n')

% Time (complete wave data)
% Assuming that the .wad file includes the time data
yearsw=wad(:,3);
monthsw=wad(:,1);
daysw=wad(:,2);
hoursw=wad(:,4);
minsw=wad(:,5);
secsw=wad(:,6);
tw=datenum(yearsw,monthsw,daysw,hoursw,minsw,secsw); % Initial time for each burst

% Time (initial datum for each burst; .whd 'always' has the time in it, but
%   DOUBLE CHECK!)
yearswb=whd(:,3);
monthswb=whd(:,1);
dayswb=whd(:,2);
hourswb=whd(:,4);
minswb=whd(:,5);
secswb=whd(:,6);
twb=(datenum(yearswb,monthswb,dayswb,hourswb,minswb,secswb))'; % Initial time for each burst

% Data collection parameters
Nb=length(whd(:,1)); % Number of bursts

twvec=(0:1/freq:Nw/freq-1/freq)/(3600*24); % Time vector for each burst (in days)

% Preallocation
wave.P=nan(Nw,Nb); % Pressure recorded at bursts
wave.t=nan(Nw,Nb); % Time for each wave datum
wave.u=nan(Nw,Nb); % Velocity (East)
wave.v=nan(Nw,Nb); % Velocity (North)
wave.w=nan(Nw,Nb); % Velocity (Up)
wave.RSSI1=nan(Nw,Nb); % RSSI Amplitude (Beam 1)
wave.RSSI2=nan(Nw,Nb); % RSSI Amplitude (Beam 2)
wave.RSSI3=nan(Nw,Nb); % RSSI Amplitude (Beam 3)

% Heading
wave.heading=whd(:,12)';

% Pitch
wave.pitch=whd(:,13)';

% Roll
wave.roll=whd(:,14)';

% Voltage
wave.volt=whd(:,10)';


% 'Old' data arrangement: .wad file without time.
% 'New' data arrangement: .wad file including time at the beginning of each
%   row of data.
for i=1:Nb
    if round(i/Nb*100*100)/100 == 25 || round(i/Nb*100*100)/100 == 50 || round(i/Nb*100*100)/100 == 75 || round(i/Nb*100*100)/100 == 100
        fprintf(['Creating indexing for wave structure array: ' ...
            num2str(round(i/Nb*100*100)/100) '%% \n'])
    end
        
    % Generate an index and look for the times within each burst
    % For the last burst:
    if i==Nb
        indexburst{i}=find(tw>=twb(i));
    else
        % For every burst except the last one:
        indexburst{i}=find(tw>=twb(i)&tw<twb(i+1));
    end
end
   
if tw(1)==twb(1) % .wad and .whd files include the time data
    for i = 1:Nb
        if round(i/Nb*100*100)/100 == 25 || round(i/Nb*100*100)/100 == 50 || round(i/Nb*100*100)/100 == 75 || round(i/Nb*100*100)/100 == 100
            fprintf(['Creating wave structure array: ' ...
                num2str(round(i/Nb*100*100)/100) '%% \n'])
        end
        
        % Extract the variables associated with the index (concatenate NaNs if
        %   there is an incomplete burst)
        wave.t(:,i)=[tw(indexburst{i});nan(Nw-length(indexburst{i}),1)]; % Time in Matlab days
        wave.P(:,i)=[wad(indexburst{i},7);nan(Nw-length(indexburst{i}),1)]; % Pressure (dbar)
        wave.STrk(:,i) = ...
            ([wad(indexburst{i},8);nan(Nw-length(indexburst{i}),1)]+...
            [wad(indexburst{i},9);nan(Nw-length(indexburst{i}),1)])/2; % Acoustic Surface Tracking distance, m
        wave.u(:,i)=[wad(indexburst{i},12);nan(Nw-length(indexburst{i}),1)]; % u (East), m/s
        wave.v(:,i)=[wad(indexburst{i},13);nan(Nw-length(indexburst{i}),1)]; % v (North), m/s
        wave.w(:,i)=[wad(indexburst{i},14);nan(Nw-length(indexburst{i}),1)]; % w (Up), m/s
        wave.RSSI1(:,i)=[wad(indexburst{i},15);nan(Nw-length(indexburst{i}),1)]; % echo amplitude (East), in counts
        wave.RSSI2(:,i)=[wad(indexburst{i},16);nan(Nw-length(indexburst{i}),1)]; % echo amplitude (North), in counts
        wave.RSSI3(:,i)=[wad(indexburst{i},17);nan(Nw-length(indexburst{i}),1)]; % echo amplitude (Up), in counts
    end
    
else % .wad file DOES NOT include the time for each datum of wave data
    % This option has not been used for BOEM project deployments
    
    for i=1:Nb
        if round(i/Nb*100*100)/100 == 25 || round(i/Nb*100*100)/100 == 50 || round(i/Nb*100*100)/100 == 75 || round(i/Nb*100*100)/100 == 100
            fprintf(['Creating wave structure array: ' ...
                num2str(round(i/Nb*100*100)/100) '%% \n'])
        end
    
        wave.t(:,i) = twb(1,i)+twvec'; % Time in Matlab days
        wave.P(:,i) = [wad(indexburst{i},3);nan(Nw-length(indexburst{i}),1)]; % Pressure (dbar)
        wave.STrk(:,i) = ...
            ([wad(indexburst{i},8);nan(Nw-length(indexburst{i}),1)]+...
            [wad(indexburst{i},9);nan(Nw-length(indexburst{i}),1)])/2; % Acoustic Surface Tracking distance, m
        wave.u(:,i) = [wad(indexburst{i},6);nan(Nw-length(indexburst{i}),1)]; % u (East), m/s
        wave.v(:,i) = [wad(indexburst{i},7);nan(Nw-length(indexburst{i}),1)]; % v (North), m/s
        wave.w(:,i) = [wad(indexburst{i},8);nan(Nw-length(indexburst{i}),1)]; % w (Up), m/s
        wave.RSSI1(:,i) = [wad(indexburst{i},10);nan(Nw-length(indexburst{i}),1)]; % echo amplitude (Beam 1), in counts
        wave.RSSI2(:,i) = [wad(indexburst{i},11);nan(Nw-length(indexburst{i}),1)]; % echo amplitude (Beam 2), in counts
        wave.RSSI3(:,i) = [wad(indexburst{i},12);nan(Nw-length(indexburst{i}),1)]; %  echo amplitude (Beam 3), in counts
    end
end

wave.tJ=wave.t-datenum(yearswb(1),0,0);

fprintf('Waves data retrieved...\n')

%% Waves data
save([pathname filename '-rawXXX.mat'],'curr','wave','coord');

% Show messages about finished job and elapsed time
endt=toc;

fprintf(['Raw data for AWAC loaded and ' filename '-raw.mat file create\n'])
fprintf(['Elapsed time: ' num2str(endt/60) ' min\n'])
