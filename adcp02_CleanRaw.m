%% Clean Nortek instruments (AWAC and Aquadopp) raw data
close all,clear all
% After storing them in a .mat file using 'maindata' and
%   'loadrawdata_Nortek' scripts.

%% Load data
% [pathname,filename]=filepath_data(a,b,c,d,e);

% Define pathname and filename for non-BOEM data campaign (e.g. Matanzas)
filename='mtz09a01'; which_instr='AqDp'; % AqDp deployment Feb. 2009
% filename='mtz09B01'; which_instr='AWAC'; % AWAC deployment Feb. 2009
pathname=['/Users/pna/Dropbox/Data/Data02_ADCP/20090212_MtnzsDplmnt/' which_instr '/RawData/'];
data_out_dir=['/Users/pna/Dropbox/Publications/JournalManuscripts/MsXX_MtnzsCurrents/DataFiles/'];
fig_out_dir=['/Users/pna/Dropbox/Publications/JournalManuscripts/MsXX_MtnzsCurrents/Figs/' which_instr '/'];
adcp02a_ObtainHdrInfo;

load([pathname filename '-raw.mat']); % Raw data    % Data storage
wavep=load([pathname filename '.wap']); % Processed waves data

%% Rename un-cleaned data from currents
% Cell positions
z_p0=curr.z_p;
% Time
t0=curr.t; tJ0=curr.tJ;
% Pressure (dbar)
P0=curr.P; 
% Temperature (°C)
T0=curr.T; 
% Heading, pitch, roll (degrees)
heading0=curr.heading; pitch0=curr.pitch; roll0=curr.roll;
% Voltage
volt0=curr.volt;
% Velocities in (m/sec)
u0=curr.u; v0=curr.v; w0=curr.w;

%% Examine pressure data to select 'initial' and 'final' measurements
% Can be from pressure (exclude measurements in air) or at specific times. 
% If necessary to delimit both ends, change 'end' to whatever value in 
% 'good1(1:end)'.
figure(1)
set(gcf,'position',[821-800 800 800 600])
plot(P0,'-r'),hold on
plot(P0,'.b'),grid on
set(gca,'YLim',[0 1.2*max(P0)])
title('Current Measurements - Mean Pressure Time Series')
xlabel(['Current Measurement Interval (' num2str(ProfileInterval/60) ' min spacing)'])
ylabel('(dbars)')

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpdf ' fig_out_dir 'Fig_01.pdf'];
eval(prtstr)
%% Select data within 'initial' and 'final' times by inspecting pressure data
adcp02b_SelectTempCurr;

% Assign new variable names to the data that will be analyzed
t=t0(good1);
tJ=tJ0(good1);
P=P0(good1);
T=T0(good1);
heading=heading0(good1);
pitch=pitch0(good1);
roll=roll0(good1);
volt=volt0(good1);

u1=u0(good1,:);
v1=v0(good1,:);
w1=w0(good1,:);

%% Examine u-velocity data to select 'upper' and 'lower' boundaries
figure(2)
set(gcf,'position',[821   800   800   600])
pcolor(tJ,z_p0,u1'),shading flat,hold on
plot(tJ,P,'r','linewidth',1.5)
set(gca,'XLim',[min(tJ) max(tJ)],'YLim',[0 1.2*max(P)],'CLim',[-0.5 0.5])
title('Uncleaned u velocity (east-west)')
xlabel('Day of Year')
ylabel('Elev. above sensor (m)')
colorbar

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpng -r600 ' fig_out_dir 'Fig_02.png'];
eval(prtstr)

%% Select usable values from water column (R=cos(theta)*D)
adcp02c_SelectVertCurr;

z_p=z_p0(good2); % Distance of each cell from instrument
u=u1(:,good2);
v=v1(:,good2);
w=w1(:,good2);

figure(3)
set(gcf,'position',[821+800 800 800 600])
plot(tJ,P,'r','linewidth',1.5),hold on
pcolor(tJ,z_p,u'),shading flat
colormap('jet')
set(gca,'XLim',[min(tJ) max(tJ)],'YLim',[0 1.2*max(P)],'CLim',[-0.5 0.5])
title('Cleaned u velocity (east-west)')
xlabel('Day of Year')
ylabel('Elev. above sensor (m)')
colorbar

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpng -r600 ' fig_out_dir 'Fig_03.png'];
eval(prtstr)

%% Clean the RSSI matrices
% % % RSSI0=RSSI;
% % % RSSI=nan(length(good1),length(good2));
% % % % Received Signal Strength Indicator
% % % for j=1:length(good1)
% % %     RSSIinst=[RSSI0.a1(good1(j),good2);RSSI0.a2(good1(j),good2);RSSI0.a3(good1(j),good2)];
% % %     RSSI(j,:)=mean(RSSIinst);
% % % end

% Store time series of RSSI from each transducer
% Received Signal Strength Indicator (RSSI), in counts
currclean.RSSI1=curr.RSSI1(good1,good2);
currclean.RSSI2=curr.RSSI2(good1,good2);
currclean.RSSI3=curr.RSSI3(good1,good2);

%% Plot cleaned RSSI data for each of the 3 beams
figure(4)
set(gcf,'position',[821-800 800-750 800 600])
minRSSI=min([min(min(currclean.RSSI1)) min(min(currclean.RSSI2)) min(min(currclean.RSSI3))]);
maxRSSI=max([max(max(currclean.RSSI1)) max(max(currclean.RSSI2)) max(max(currclean.RSSI3))]);

subplot(3,1,1)
plot(tJ,P),hold on
pcolor(tJ,z_p,currclean.RSSI1'),shading interp
set(gca,'XLim',[min(tJ) max(tJ)],'YLim',[0 1.2*max(P)],'CLim',[minRSSI maxRSSI])
title('Echo Intensity - I1')
ylabel('Elev. above sensor (m)')

subplot(3,1,2)
plot(tJ,P),hold on
pcolor(tJ,z_p,currclean.RSSI2'),shading interp
set(gca,'XLim',[min(tJ) max(tJ)],'YLim',[0 1.2*max(P)],'CLim',[minRSSI maxRSSI])
title('Echo Intensity - I2')
ylabel('Elev. above sensor (m)')

subplot(3,1,3)
plot(tJ,P),hold on
pcolor(tJ,z_p,currclean.RSSI3'),shading interp
title('Echo Intensity - I3')
set(gca,'XLim',[min(tJ) max(tJ)],'YLim',[0 1.2*max(P)],'CLim',[minRSSI maxRSSI])
xlabel('Day of Year')
ylabel('Elev. above sensor (m)')

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpng -r600 ' fig_out_dir 'Fig_04.png'];
eval(prtstr)

%% Plot time series of instrument orientations - heading, pitch, and roll 
figure(5)
set(gcf,'position',[821 800-750 800 600])

subplot(3,1,1)
plot(tJ,heading,'.m'),grid on
set(gca,'XLim',[min(tJ) max(tJ)],'YLim',[-100 100])
title('Heading')
ylabel('Degrees')

subplot(3,1,2)
plot(tJ,pitch,'.g'),grid on
set(gca,'XLim',[min(tJ) max(tJ)],'YLim',[-3 3])
title('Pitch')
ylabel('Degrees')

subplot(3,1,3)
plot(tJ,roll,'.r'),grid on
set(gca,'XLim',[min(tJ) max(tJ)],'YLim',[-10 10])
title('Roll')
xlabel('Day of Year')
ylabel('Degrees')

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpdf ' fig_out_dir 'Fig_05.pdf'];
eval(prtstr)

%% Plot cleaned u and v velocities, depth-irrespective (not pcolor)
figure(6)
set(gcf,'position',[821+800 800-750 800 600])

subplot(2,1,1)
plot(tJ,u,'.k'),grid on
set(gca,'XLim',[min(tJ) max(tJ)],'YLim',[-1 1])
title('u velocity (east-west)')
ylabel('(m/s)')

subplot(2,1,2)
plot(tJ,v,'.b'),grid on
set(gca,'XLim',[min(tJ) max(tJ)],'YLim',[-1 1])
title('v velocity (north-south)')
xlabel('Day of Year')
ylabel('(m/s)')

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpdf ' fig_out_dir 'Fig_06.pdf'];
eval(prtstr)

%% Build strucuture array with cleaned data
currclean.t=t; % Time
currclean.tJ=tJ; % Time (Julian days)
currclean.u=u; % East-West velocity
currclean.v=v; % Nort-South velocity
currclean.w=w; % Up-down velocity
currclean.P=P; % Water depth (pressure)
currclean.T=T; % Temperature (°C)
currclean.heading=heading; % Temperature (°C)
currclean.pitch=pitch; % Temperature (°C)
currclean.roll=roll; % Temperature (°C)
currclean.volt=volt; % Voltage (V)
currclean.z_p=z_p; % Distance from instrument (m)

%% WAVES DATA =======================================================
% % % wave_clean=wave.P.*(wave.P>5); % Look for pressure less than 5 m
% % % wave_clean(wave_clean==0)=NaN;
% % % indexNaN=isnan(wave_clean);
% % % colsNaN=sum(indexNaN,1);
% % % 
% % % % Create an index to search for at least one NaN within each burst
% % % indexwaveclean=find(colsNaN==0);
% % % 
% % % pcolor(wave.P),shading interp
% % % 
figure(7)
set(gcf,'position',[821-800 800 800 600])
plot(wave.P(end,:),'-b'),hold on
plot(wave.P(end,:),'.r'),grid on
title('Wave Burst Summary - Time Series of final pressure obs. in each burst ')
xlabel('Wave Burst Number (hourly spacing)')
ylabel('dbars')

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpdf ' fig_out_dir 'Fig_07.pdf'];
eval(prtstr)

%% Select wave data within 'initial' and 'final' times by inspecting pressure data
adcp02d_SelectTempWave;

% Assign new variable names to the cleaned wave data
waveclean.P=wave.P(:,indexwaveclean);
waveclean.heading=wave.heading(:,indexwaveclean);
waveclean.pitch=wave.pitch(:,indexwaveclean);
waveclean.roll=wave.roll(:,indexwaveclean);
waveclean.volt=wave.volt(:,indexwaveclean);

waveclean.u=wave.u(:,indexwaveclean); waveclean.u(abs(waveclean.u)>1)=NaN;
waveclean.v=wave.v(:,indexwaveclean); waveclean.v(abs(waveclean.v)>1)=NaN;
waveclean.w=wave.w(:,indexwaveclean); waveclean.w(abs(waveclean.w)>1)=NaN;

waveclean.RSSI1=wave.RSSI1(:,indexwaveclean); waveclean.w(abs(waveclean.w)>1)=NaN;
waveclean.RSSI2=wave.RSSI2(:,indexwaveclean); waveclean.w(abs(waveclean.w)>1)=NaN;
waveclean.RSSI3=wave.RSSI3(:,indexwaveclean); waveclean.w(abs(waveclean.w)>1)=NaN;

%% Surface Tracking Work (omitted here)
% % If AWAC
% waveclean.STrk=wave.STrk(:,indexwaveclean);
% waveclean.STrk(detrend(waveclean.STrk)>3)=NaN;
% waveclean.STrk(detrend(waveclean.STrk)<-3)=NaN;
% 
% figure
% plot(detrend(waveclean.STrk),'.')
% 
% % Spring 2014 - Canaveral swale east (AWAC)
% waveclean.STrk=waveclean.STrk(detrend(waveclean.STrk)>-0.5 & detrend(waveclean.STrk)<0.7);
% 
% % Summer A 2015 - Chester Swale West (AWAC)
% waveclean.STrk=waveclean.STrk(detrend(waveclean.STrk)>-0.5 & detrend(waveclean.STrk)<0.7);
% 
% % Fall B 2015 - Canaveral Swale West (AWAC)
% waveclean.STrk=waveclean.STrk(detrend(waveclean.STrk)>-1 & detrend(waveclean.STrk)<1);
% 
% figure
% plot(waveclean.STrk,'.')

%% Store Average RSSI for each burst of waves (omitted here)
% % Average received signal strength indicator per datum within each burst
% RSSIwav=nan(length(wave.a1(:,1)),length(indexwaveclean));
% % Received Signal Strength Indicator
% for k=1:length(indexwaveclean)
%     rows=[wave.RSSI1(:,indexwaveclean(k)),wave.RSSI2(:,indexwaveclean(k)),wave.RSSI3(:,indexwaveclean(k))];
%     RSSIwav(:,k)=mean(rows,2);
% end

figure(8)
set(gcf,'position',[821 800 800 600])
pcolor(waveclean.P),shading interp
title('Cleaned Wave Burst Individual Pressure Measurements')
xlabel('Wave Burst Number (hourly spacing)')
ylabel('Observation Number Within Burst (2Hz sampling rate)')

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpng -r300 ' fig_out_dir 'Fig_08.png'];
eval(prtstr)

figure(9)
set(gcf,'position',[821+800 800 800 600])
pcolor(waveclean.u),shading interp
title('Wave Burst Individual U Velocity Measurements')
xlabel('Wave Burst Number (hourly spacing)')
ylabel('Observation Number Within Burst (2Hz sampling rate)')

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpng -r300 ' fig_out_dir 'Fig_09.png'];
eval(prtstr)
% figure
% pcolor(RSSIwav),shading interp
% % If AWAC
% figure
% pcolor(waveclean.STrk),shading interp

Pmean=mean(waveclean.P); % Mean pressure for each burst

figure(10)
set(gcf,'position',[821-800 800-750 800 600])
plot(Pmean,'-b'),hold on
plot(Pmean,'.r'),grid on % Plot mean pressure for each burst
title('Cleaned Wave Burst Mean Pressure Measurements')
xlabel('Wave Burst Number (hourly spacing)')
ylabel('(dbars)')

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpdf ' fig_out_dir 'Fig_10.pdf'];
eval(prtstr)
% Select final range of values to store in 'goodP'
% % % goodP=1:length(indexwaveclean);
n1=1;
nend=length(indexwaveclean);
goodP=n1:1:nend;

% Pmean=mean(waveclean.P); % Mean pressure for each burst

% figure('position',[821 800-750 800 600])
% plot(Pmean) % Plot mean pressure for each burst
% 
% figure('position',[821+800 800-750 800 600])
% pcolor(waveclean.P),shading interp
% figure
% pcolor(RSSIwav),shading interp

% Select bursts with good pressure data (make sure that the pressure data
% is good, zoom at the previous plot)
waveclean.P=waveclean.P(:,goodP);
waveclean.t=waveclean.t(:,goodP);
waveclean.u=waveclean.u(:,goodP);
waveclean.v=waveclean.v(:,goodP);
waveclean.w=waveclean.w(:,goodP);
waveclean.RSSI1=wave.RSSI1(:,goodP);
% waveclean.RSSIm=RSSIwav(:,goodP);
% waveclean.STrk=waveclean.STrk(:,goodP); % If AWAC
disp('processed waves')
%% PROCESSED WAVES (after using Quick Waves software) DATA ================
% Select only valid data (look at pressure plot)
if length(wavep(1,:))==23 % AWAC
    disp('AWAC')
    figure
    subplot(2,1,1)
    plot(wavep(:,18),'-b'),hold on
    plot(wavep(:,18),'.r')
elseif length(wavep(1,:))==17 % Aquadopp
    disp('Aquadopp')
    figure
    subplot(2,1,1)
    plot(wavep(:,14),'-b'),hold on
    plot(wavep(:,14),'.r')
end

%% Select QUICKwave data within 'initial' and 'final' times by linking to previous
adcp02e_SelectTempQuickWave;

if length(wavep(1,:))==23 % AWAC
    % Time vector
    YYYYmin=min(wavep(:,3));
    YYYY=wavep(:,3);
    MM=wavep(:,1);
    DD=wavep(:,2);
    hh=wavep(:,4)+h_corr;
    mm=wavep(:,5);
    ss=wavep(:,6);
    wp_t=datenum(YYYY,MM,DD,hh,mm,ss);
    wp_tJ=wp_t-datenum(YYYYmin,0,0);
    wavepclean.Hs=wavep(indexgoodwp,7); % Significant height (m)
        wavepclean.Hs(wavepclean.Hs<0)=NaN;
    wavepclean.Tp=wavep(indexgoodwp,12); % Peak period (sec)
        wavepclean.Tp(wavepclean.Tp<0)=NaN;
    wavepclean.Dp=wavep(indexgoodwp,14); % Peak direction (°)
        wavepclean.Dp(wavepclean.Dp<0)=NaN;
    wavepclean.P=wavep(indexgoodwp,18); % Mean pressure (dbar)
    wavepclean.t=wp_t(indexgoodwp); % Time (MATLAB days)
    wavepclean.tJ=wp_tJ(indexgoodwp); % Time (Julian days)
elseif length(wavep(1,:))==17 % Aquadopp
    % Time vector
    YYYYmin=min(wavep(:,3));
    YYYY=wavep(:,3);
    MM=wavep(:,1);
    DD=wavep(:,2);
    hh=wavep(:,4)+h_corr;
    mm=wavep(:,5);
    ss=wavep(:,6);
    wp_t=datenum(YYYY,MM,DD,hh,mm,ss);
    wp_tJ=wp_t-datenum(YYYYmin,0,0);
    wavepclean.Hs=wavep(indexgoodwp,7); % Significant height (m)
        wavepclean.Hs(wavepclean.Hs<0)=NaN;
    wavepclean.Tp=wavep(indexgoodwp,9); % Peak period (sec)
        wavepclean.Tp(wavepclean.Tp<0)=NaN;
    wavepclean.Dp=wavep(indexgoodwp,10); % Peak direction (°)
        wavepclean.Dp(wavepclean.Dp<0)=NaN;
    wavepclean.P=wavep(indexgoodwp,14); % Mean pressure (dbar)    
    wavepclean.t=wp_t(indexgoodwp); % Time (MATLAB days)
    wavepclean.tJ=wp_tJ(indexgoodwp); % Time (Julian days)
end

figure(11)
set(gcf,'position',[821 800-750 800 600])
subplot(2,1,2)
plot(wavepclean.tJ,wavepclean.P,'b')
xlabel('Day of Year')

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpdf ' fig_out_dir 'Fig_11.pdf'];
eval(prtstr)

figure(12)
set(gcf,'position',[821+800 800-750 800 600])
subplot(3,1,1)
plot(wavepclean.tJ,wavepclean.Hs,'r'),grid on
title('Significant Wave Height')
ylabel('m')

subplot(3,1,2)
plot(wavepclean.tJ,wavepclean.Tp,'.b'),grid on
title('Dominant Wave Period')
ylabel('seconds')

subplot(3,1,3)
plot(wavepclean.tJ,wavepclean.Dp,'.g'),grid on
title('Peak Wave Direction')
xlabel('Day of Year')
ylabel('Degrees')

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpdf ' fig_out_dir 'Fig_12.pdf'];
eval(prtstr)

%% Coordinates
coord=load([pathname filename '_coord.txt']); % Coordinates

%% Save fixed file
save([data_out_dir filename '-clean'],'currclean','waveclean','wavepclean','coord');
