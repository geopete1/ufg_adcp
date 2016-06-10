close all,clear all
% Sensor='a'; InstrName='Aquadopp';
Sensor='B'; InstrName='AWAC';
DataDir='/Users/pna/Dropbox/Publications/JournalManuscripts/MsXX_MtnzsCurrents/DataFiles/';
% DataDir='DataFiles/';

load([DataDir 'mtz09' Sensor '01-clean.mat']);

load('/Users/pna/Dropbox/Data/Data02_ADCP/20090212_MtnzsDplmnt/AqDp/AqDp_0907_WvHtPress.mat')

if strcmp(InstrName,'Aquadopp')
    fig_out_dir='/Users/pna/Dropbox/Publications/JournalManuscripts/MsXX_MtnzsCurrents/Figs/AqDp/';
elseif strcmp(InstrName,'AWAC')
    fig_out_dir='/Users/pna/Dropbox/Publications/JournalManuscripts/MsXX_MtnzsCurrents/Figs/AWAC/';
end

timewindow=[datenum('Feb-12-2009') datenum('Apr-24-2009')];

fs=2;
findNaNs=isnan(wlevel);
a=diff(findNaNs);

starts=(find(a<0))+1;
ends=find(a>0);

figure
set(gcf,'units','inches','position',[1 1 8 14])

for i=1:length(starts)
    intrvl_edge_inds(i,:)=[starts(i) ends(i)];
    w=detrend(wlevel(starts(i):ends(i)));
    [pxx(:,i),f]=pwelch(w,[],[],[],fs);
   
%     subplot(3,1,1)
%     plot(dateWaves,wlevel)
%     set(gca,'YLim',[-1.2 1.2],'XLim',[min(dateWaves) max(dateWaves)])
%     grid on
%     ylabel('Water level (m)')
%     datetick('x')
%     title(['Water Levels - ' InstrName])

    subplot(2,1,1)
    plot(dateWaves(starts(i):ends(i)),wlevel(starts(i):ends(i)))
    set(gca,'YLim',[-1.2 1.2])
    grid on
    ylabel('Water level (m)')
    datetick('x')
    title(datestr(dateWaves(starts(i)),'dd-mmm-yyyy'))
    
    subplot(2,1,2)
    semilogx(f(2:end),pxx(2:end,i),'b','linewidth',2),grid on
    set(gca,'YLim',[0 1])
    xlabel('Frequency (Hz)')
    ylabel('Spectral Power')
    drawnow
    pause(0.01)
end

%% Draw Spectrogram
timewindow2=[datenum('Mar-18-2009') datenum('Mar-31-2009')];

figure

h1a=axes;
pcolor(dateWaves(starts),f(3:end),pxx(3:end,:))
pos = [timewindow2(1) 0 timewindow2(2)-timewindow2(1) 0.05];
rectangle('Position',pos,'EdgeColor','r','linewidth',2)
shading interp
set(gca,'CLim',[0 1.3],'XLim',timewindow,'YLim',[0 0.3],'fontSize',14)
ylabel('Frequency (Hz)')
datetick('x','keepticks','keeplimits')

h1b=axes;
PdYaxTicks=get(h1a,'YTick');
PdYaxLims=get(h1a,'YLim');
PdYaxLabels=num2str(round(10*(1./str2num(char(get(h1a,'YTickLabel')))))/10);
set(h1b,'Color','None','YAxisLocation','Right',...
    'XLim',[min(dateWaves(starts)) max(dateWaves(starts))],'YLim',PdYaxLims,...
    'XTick',[],'YTick',PdYaxTicks,'YTickLabel',PdYaxLabels,'fontSize',14)
ylabel('Wave Period (s)')

set(gcf,'units','inches','papersize',[8 3],'paperposition',[0 0 8 2.5])
prtstr=['print -dpng -r300 ' fig_out_dir 'Fig_20_Spectrogram' InstrName '.png'];
eval(prtstr)

figure

h2=axes;
pcolor(dateWaves(starts),f(2:end),pxx(2:end,:))
shading interp
set(gca,'CLim',[0 0.05],'XLim',timewindow2,'YLim',[0 0.05],'fontSize',14)
PdYaxLabels=num2str(round(10*(1./str2num(char(get(h2,'YTickLabel')))))/10);
set(gca,'YTickLabel',[PdYaxLabels])
ylabel('Wave Period (s)')
datetick('x','keepticks','keeplimits')

set(gcf,'units','inches','papersize',[8 3],'paperposition',[0 0 8 2.5])
prtstr=['print -dpng -r300 ' fig_out_dir 'Fig_21_InfragravSpectrogram' InstrName '.png'];
eval(prtstr)
