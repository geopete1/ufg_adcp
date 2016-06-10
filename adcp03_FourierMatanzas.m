close all,clear all
% Sensor='a'; InstrName='Aquadopp';
Sensor='B'; InstrName='AWAC';
DataDir='/Users/pna/Dropbox/Publications/JournalManuscripts/MsXX_MtnzsCurrents/DataFiles/';
% DataDir='DataFiles/';

if strcmp(InstrName,'Aquadopp')
    fig_out_dir='/Users/pna/Dropbox/Publications/JournalManuscripts/MsXX_MtnzsCurrents/Figs/AqDp/';
elseif strcmp(InstrName,'AWAC')
    fig_out_dir='/Users/pna/Dropbox/Publications/JournalManuscripts/MsXX_MtnzsCurrents/Figs/AWAC/';
end

load([DataDir 'mtz09' Sensor '01-clean.mat']);

P=currclean.P;
% Demean the Pressure Data
P_dm=P-mean(P);

% Select which treatment for the data: demeaned, detrended, or both
P_dmdt=detrend(P_dm); titlestr=['FFT of Detrended and Demeaned'];treatmentStr='DeMeDeTr';
% P_dmdt=P_dm; titlestr=['FFT of Demeaned Only (Not Detrended)']; treatmentStr='DeMeUnDeTr';
% P_dmdt=P; titlestr='FFT of Original Data (Not Detrended and Not Demeaned)'; treatmentStr='UnDeMeUnDeTr';

mb=polyfit(currclean.tJ,P,1);
P_tr=polyval(mb,currclean.tJ);

%% Make plot of pressure time series
figure
subplot(3,1,1),plot(currclean.tJ,P),hold on
plot(currclean.tJ,ones(length(P),1)*mean(P),'r','linewidth',1.5)
plot(currclean.tJ,P_tr,'g--','linewidth',1.5),grid on
title('Original Pressure Time Series')
ylabel('dbar')
subplot(3,1,2),plot(currclean.tJ,P_dm),grid on
ylabel('dbar')
title('Demeaned Pressure Time Series')
subplot(3,1,3),plot(currclean.tJ,P_dmdt),grid on
title('Demeaned, Detrended Pressure Time Series')
ylabel('dbar')
xlabel('Day of Year')
suptitle([InstrName])

set(gcf,'units','inches','papersize',[5 5],'paperposition',[0 0 5 5])
prtstr=['print -dpdf ' fig_out_dir 'Fig_13_MtzPressureTS' InstrName '.pdf'];
eval(prtstr)
%% Insert Task Here: Code the Discrete Fourier Transform "by hand"

%% Calculate and plot FFTs
L=length(P_dmdt);
NFFT=2^nextpow2(L);
Y=fft(P_dmdt,NFFT)/L;
if strcmp(InstrName,'Aquadopp')
    Fs=48; % 48 samples per day, 30 min. sampling frequency
elseif strcmp(InstrName,'AWAC')
    Fs=144; % 144 samples per day, 10 min. sampling frequency
end
f_nyq=Fs/2; % Nyquist Frequency
f=f_nyq*linspace(0,1,NFFT/2+1);

figure
subplot(2,1,1)
plot(f,2*abs(Y(1:NFFT/2+1)))
set(gca,'XLim',[0 25])
grid on
xlabel('frequency (cpd)')
title(titlestr)

subplot(2,1,2)
plot(f,2*abs(Y(1:NFFT/2+1)))
set(gca,'XLim',[0 4])
grid on
xlabel('frequency (cpd)')
title('Same as above, but zoomed in on x axis limits')

suptitle([InstrName])

set(gcf,'units','inches','papersize',[6 6],'paperposition',[0 0 6 6])
prtstr=['print -dpdf ' fig_out_dir 'Fig_14_MtzPressureFFT' InstrName treatmentStr '.pdf'];
eval(prtstr)

%% Calculate the Fourier Coefficients
N=length(P_dmdt);
tic
for n=1:N/2
    for j=1:N
        a(n,j)=cos(2*pi*(n-1)*j/N);
        b(n,j)=sin(2*pi*(n-1)*j/N);
%         disp(['n=' num2str(n) ', j=' num2str(j) ', a(n,j)=' num2str(a(n,j))])
    end
    A(n)=(2/N)*sum(P_dmdt.*a(n,:)');
    B(n)=(2/N)*sum(P_dmdt.*b(n,:)');
%     disp(['Freq. loop ' num2str(n) ' of ' num2str(N/2)]);
end
toc
harm=[0:1:n-1];
C=sqrt((A.^2)+(B.^2));
Pd_hrs=24./(harm./(length(P_dmdt)/Fs));
Phase=atan(B./A);

% Plot the Fourier Coefficients
figure
subplot(2,1,1),hold on
plot(harm,A,'b')
plot(harm,B,'r')
plot(harm,C,'g')
grid on
legend('A','B','C')
xlabel('harmonic number')
ylabel('dbar')

subplot(2,1,2)
semilogx(Pd_hrs,A,'b'),hold on
semilogx(Pd_hrs,B,'r')
semilogx(Pd_hrs,C,'g')
grid on
xlabel('Periodicity (hrs)')
ylabel('dbar')
suptitle(['Fourier Coefficients - ' InstrName])

set(gcf,'units','inches','papersize',[6 6],'paperposition',[0 0 6 6])
prtstr=['print -dpdf ' fig_out_dir 'Fig_15_MtzPressureFourierCoeffs' InstrName treatmentStr '.pdf'];
eval(prtstr)

% Plot the Phase
figure
subplot(2,1,1),hold on
plot(harm,Phase)
grid on
ylabel('Radians')
xlabel('harmonic number')

subplot(2,1,2)
semilogx(Pd_hrs,Phase)
grid on
xlabel('Periodicity (hrs)')
ylabel('Radians')

suptitle(['Phase Angle - ' InstrName])

set(gcf,'units','inches','papersize',[6 6],'paperposition',[0 0 6 6])
prtstr=['print -dpdf ' fig_out_dir 'Fig_16_MtzPressurePhaseAngle' InstrName treatmentStr '.pdf'];
eval(prtstr)

%% Reconstruct the Pressure Time series
ncomp=100;
[Csorted,C_ind]=sort(C,'descend');

for j=1:N
    for n=1:ncomp
        Comp(j,n)=C(C_ind(n))*cos((2*pi*(C_ind(n)-1)*j/N)-Phase(C_ind(n)));
        CompAlt(j,n)=A(C_ind(n))*cos((2*pi*(C_ind(n)-1)*j/N))+B(C_ind(n))*sin((2*pi*(C_ind(n)-1)*j/N));
    end
    P_recon(j)=(0.5*C(1))+sum(CompAlt(j,:));
end

figure('position',[20 300 800 600])
subplot(2,1,1)
plot(currclean.tJ,P_dmdt,'b','linewidth',2),hold on
plot(currclean.tJ,P_recon,'r--','linewidth',2),grid on
ylabel('Pressure deviation from mean (dbar)')
subplot(2,1,2)
plot(currclean.tJ,P_dmdt,'b','linewidth',2),hold on
plot(currclean.tJ,P_recon,'r--','linewidth',2),grid on
set(gca,'XLim',[70 80])
xlabel('Day of Year')
ylabel('Pressure deviation from mean (dbar)')

suptitle(['Reconstruction using ' num2str(ncomp) ' Harmonic Components - ' InstrName])

set(gcf,'units','inches','papersize',[6 6],'paperposition',[0 0 6 6])
prtstr=['print -dpdf ' fig_out_dir 'Fig_17_MtzPressureReconst' InstrName treatmentStr num2str(ncomp) 'comp.pdf'];
eval(prtstr)

figure
plot(P_dmdt,P_recon,'b*'),grid on
set(gca,'XLim',[-1.5 1.5],'YLim',[-1.5 1.5])
xlabel('Measured Pressure')
ylabel('Reconstructed Pressure')
title(['Reconstruction Performance - ' num2str(ncomp) ' Harmonic Components - ' InstrName])

set(gcf,'units','inches','papersize',[6 6],'paperposition',[0 0 6 6])
prtstr=['print -dpdf ' fig_out_dir 'Fig_18_MtzPressureReconstFitTest' InstrName treatmentStr num2str(ncomp) 'comp.pdf'];
eval(prtstr)

