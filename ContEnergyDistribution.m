%%
clear; clear all;
x_SCP=xlsread('C:\Users\abhij\Desktop\WaveModel\Dam_Size_Anl\U2DLA0SCP42-01Dam12mm150khzStiffNM.xlsx');
x_SCS=xlsread('C:\Users\abhij\Desktop\WaveModel\Dam_Size_Anl\U2DLA0SCS42-01Dam12mm150khzStiffNM.xlsx');
x_CFN=xlsread('C:\Users\abhij\Desktop\WaveModel\Dam_Size_Anl\U2DLA0CFMN01-42Dam12mm150khzStiffNM.xlsx');
%%
% Calculation of (u2_plat-u2_stiff) #Displacament difference for each time step
for i=1:42
    x_diff=x_SCP(2:end,i+1)-x_SCS(2:end,i+1);
    X_Diff(:,i)=x_diff;
end
%%
% Energy calculation for each time step
for j=1:42
    energy_t = x_CFN(2:end,j+1).*(abs(X_Diff(:,43-j)));
    ENERGY_t(:,j)=energy_t;
end
%%
figure
for j=1:42
    hold on
    plot(ENERGY_t(:,j))
end
%%
% Summation energy
for j=1:42
    enr_sum =sum(ENERGY_t(:,j));
    ENR_SUM(:,j)=enr_sum;
end
%%
ENR_Matrix = zeros(7,6);
for i=1:7
    for j=1:6
        ENR_Matrix(i,j) = ENR_SUM(:,(i-1)*6+j);
    end
end   
ENR_Matrix(1,1) = 0.75*ENR_Matrix(1,1);
%%
n = 5;
beta = 20;
ENR_MAT_A=zeros(7,12);
for i =1:7
    ENR_MAT_A(i,:)=resample(ENR_Matrix(i,:),12,6,n,beta);
end
ENR_MAT_B=zeros(13,12);
for i =1:12
    ENR_MAT_B(:,i)=resample(ENR_MAT_A(:,i),13,7,n,beta);
end
%%
ENR_MAT_AA =zeros(13,12);
for i=1:13
    ENR_MAT_AA(14-i,:)=ENR_MAT_B(i,:);
end
ENR_MAT=zeros(26,12);
ENR_MAT =[ENR_MAT_B; ENR_MAT_AA];
%%
figure
X=1:1:12;
Y=1:1:26;
Z = bar3(ENR_MAT,1.1);
% xlabel('Lenght mm')
% ylabel('Width mm')
zlabel('Energy')
for k = 1:length(Z)
    zdata = Z(k).ZData;
    Z(k).CData = zdata;
    Z(k).FaceColor = 'interp';
end
ENR_MAT_Modified =[ENR_MAT(1:12,:); ENR_MAT(15:26,:)];
figure
surf(ENR_MAT_Modified)
shading interp
%%
% Contact Force distribution
for j=1:42
    cfn_sum =sum(x_CFN(2:end,j+1));
    CFN_SUM(:,j)=cfn_sum;
end
CFN_Matrix = zeros(7,6);
for i=1:7
    for j=1:6
        CFN_Matrix(i,j) = CFN_SUM(:,(i-1)*6+j);
    end
end 
CFN_A = zeros(6,6);
for i=1:6
    CFN_A(i,:)=CFN_Matrix(7-i,:);
end
CFN_MAT = [CFN_Matrix; CFN_A];
figure 
F = bar3(CFN_MAT,1.05);
for k = 1:length(F)
    zdata = F(k).ZData;
    F(k).CData = zdata;
    F(k).FaceColor = 'interp';
end
figure
surf(CFN_MAT)
shading interp
%%
% Phase Difference
figure
c_p = 22; % contact position
plot(x_SCP(2:end,1),x_SCP(2:end,44-c_p),'b.')
hold on
plot(x_SCP(2:end,1),x_SCS(2:end,44-c_p),'r.')
%%
%// Define the x values
t=x_SCP(2:end,1).*1000000;
tMat = repmat(t,1,6); %// For plot3
%// Define y values
y = 1:6:31;
yMat = repmat(y, numel(t), 1); %//For plot3
%// Define z values
c_p1 = 1;c_p2 = 7;c_p3 = 13;c_p4 = 19;c_p5 = 25;c_p6 = 31; % contact position
zSCP1 = x_SCP(2:end,44-c_p1);
zSCP2 = x_SCP(2:end,44-c_p2);
zSCP3 = x_SCP(2:end,44-c_p3);
zSCP4 = x_SCP(2:end,44-c_p4);
zSCP5 = x_SCP(2:end,44-c_p5);
zSCP6 = x_SCP(2:end,44-c_p6);
zMat1 =[zSCP1 zSCP2 zSCP3 zSCP4 zSCP5 zSCP6];
zSCS1 = x_SCS(2:end,44-c_p1);
zSCS2 = x_SCS(2:end,44-c_p2);
zSCS3 = x_SCS(2:end,44-c_p3);
zSCS4 = x_SCS(2:end,44-c_p4);
zSCS5 = x_SCS(2:end,44-c_p5);
zSCS6 = x_SCS(2:end,44-c_p6);
zMat2 =[zSCS1 zSCS2 zSCS3 zSCS4 zSCS5 zSCS6];   %// For plot3
figure
plot3(tMat, yMat, zMat1, 'b','LineWidth',2); %// Make all traces blue
hold on
plot3(tMat, yMat, zMat2, 'g','LineWidth',2)
grid;
xlabel('Time (µ sec)'); ylabel('Location'); zlabel('Transverse Displacement');
view(40,40); %// Adjust viewing angle so you can clearly see data
%%
% Contact presure energy distribution
%// Define the x values
x_e = 1:1:12; 
x_e = x_e';
X_Mat  = repmat(x_e,1,26);%// For plot3
%// Define y values
y_e = 1:1:26;
y_e = y_e;
Y_Mat = repmat(y_e, numel(x_e), 1); %//For plot3
%// Define z values
Z_Mat = ENR_MAT_B;   %// For plot3
figure
plot3(X_Mat, Y_Mat, Z_Mat, 'b.-','LineWidth',2); %// Make all traces blue
hold on
plot3(X_Mat', Y_Mat', Z_Mat', 'g.-','LineWidth',2);
plot3(X_Mat', Y_Mat', Z_Mat', 'r*','LineWidth',2);
grid;
%xlabel('Time (µ sec)'); ylabel('Location'); zlabel('Transverse Displacement');
view(40,40); %// Adjust viewing angle so you can clearly see data
%%
% Contact Force distribution
%// Define the x values
x_e = 1:1:6; 
x_e = x_e';
X_Mat  = repmat(x_e,1,13);%// For plot3
%// Define y values
y_e = 1:1:13;
y_e = y_e;
Y_Mat = repmat(y_e, numel(x_e), 1); %//For plot3
%// Define z values
Z_Mat = CFN_MAT;   %// For plot3
figure
plot3(X_Mat, Y_Mat, Z_Mat, 'b.-','LineWidth',2); %// Make all traces blue
hold on
plot3(X_Mat', Y_Mat', Z_Mat', 'g.-','LineWidth',2);
plot3(X_Mat', Y_Mat', Z_Mat', 'r*','LineWidth',2);
grid;
%xlabel('Time (µ sec)'); ylabel('Location'); zlabel('Transverse Displacement');
view(40,40); %// Adjust viewing angle so you can clearly see data
%%
n = 5;
beta = 20;
ENR_MAT_A=zeros(13,12);
for i =1:13
    ENR_MAT_A(i,:)=resample(ENR_MAT(i,:),12,6,n,beta);
end
ENR_MAT_B=zeros(26,12);
for i =1:12
    ENR_MAT_B(:,i)=resample(ENR_MAT_A(:,i),26,13,n,beta);
end
figure
surf(abs(ENR_MAT_B(1:26,:)))%;shading interp
%%