%%
clc;clear all;
X_00=xlsread('C:\Users\abhij\Desktop\WaveModel\StiffPlate\UT2DLA0S16-1T0003_M00025.xlsx');
X_05=xlsread('C:\Users\abhij\Desktop\WaveModel\Dam5mm\UT2DLA0S1-16T0003_M00025.xlsx');
X_10=xlsread('C:\Users\abhij\Desktop\WaveModel\Dam10mm\UT2DLA0S16-1T0003_M00025.xlsx');
X_15=xlsread('C:\Users\abhij\Desktop\WaveModel\Dam15mm\UT2DLA0S1-16T0003_M00025.xlsx');
X_20=xlsread('C:\Users\abhij\Desktop\WaveModel\Dam20mm\UT2DLA0S1-16T0003_M00025.xlsx');
%X_25=xlsread('C:\Users\abhij\Desktop\WaveModel\Dam25mm\UT2DLA0S1-16T0003_M00025.xlsx');
%%
sn=9;  %Sensor Number started from 2nd row 
X_sen=[X_00(:,(19-sn)) X_05(:,sn) X_10(:,(19-sn)) X_15(:,sn) X_20(:,sn)];
%% Normalization time series sensor data in range -1 to 1
for i=1:5
    %X_scaled=rescale(X_sen(:,i),-1,1);
    X_scaled=X_sen(:,i);
    X_norm(:,i)=X_scaled;
end
X=X_norm;
fs = 1e8;
t = 0:1/fs:0.0003; 
%%
figure
plot(t,X(:,2))
%% Filtering of sensor data through HighPass Fill and then BandPass filter
load('BandPassFil_1.mat')
load('BandPassFil450kHz.mat')
load('BandPassFil450kHz_2.mat')
load('BandPassFil450kHz_3.mat')
load('HighPassFil_1.mat')
load('HighPassFil400kHz.mat')
load('HighPassFil400kHz_2.mat')

fs = 1e8;
t = 0:1/fs:0.0003; 
L = length(X);
for j=1:5
    X_1=filter(HighPassFil_1,X(:,j));      %HIGHPASS Filtering
    X_fil=filter(BandPassFil_1,X_1);       %BANDPASS Filtering
    x_fil(:,j)=X_fil;
end
%% Ploting Filtered Data
figure
for k=2:4
    hold on
    plot(t,x_fil(:,k))
end
xlabel('time')
ylabel('Normalized Amplitude')
legend('undamaged','5mm','10mm','15mm','20mm','Location','northwest')
%%
figure
plot(t,X(:,1))
hold on
plot(t,x_fil(:,1))
%% FFT of filtered sensor data

for l=1:5
    Y=fft(x_fil(:,l));
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;
    P_1(:,l)=P1;
    F(:,l)=f;
end
%% Ploting FFT
figure
for m=1:5
    hold on
    plot(F(:,m),P_1(:,m),'marker','.') 
end
%title('FFT for sensor at 330mm')
xlabel('f(Hz)')
ylabel('Amplitute')
legend('undamaged','5mm','10mm','15mm','20mm','Location','northeast')
%% Hilbert Transform
for n=1:5
    y= hilbert(x_fil(:,n));
    Y(:,n)=abs(y);
    SUM=sum(abs(y));
    DamIndex(:,n)=SUM;
end
figure
for r=1:5
    hold on 
    plot(Y(:,r),'.-')
end
xlabel('time')
ylabel('Hilbert transform coefficient')
legend('undamaged','5mm','10mm','15mm','20mm','Location','northwest')
%% Damage Indexing 
figure
dam=[0 5 10 15 20]; %Size of damage in mm "0" means undamage
plot(dam(:,1:4),DamIndex(:,1:4),'rs-')
ylabel('Damage Idex')
xlabel('Debonding Size (mm)')
figure
s=2; t=3;
dam_2=[dam(:,s) dam(:,t:4)];
DI=[DamIndex(:,s) DamIndex(:,t:4)];
plot(dam(:,2:5),DamIndex(:,2:5),'rs-')
ylabel('Damage Idex')
xlabel('Debonding Size (mm)')
%% Ploting Filtered Data
% figure
% plot(t,x_fil(:,1));
% hold on
% plot(t,x_fil(:,5));
% xlabel('time')
% ylabel('Normalized Amplitude')
% legend('undamaged','10 mm','Location','northwest')
