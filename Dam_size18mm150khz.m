%%
clear; clear all;
x_all_node_cpress=xlsread('C:\Users\abhij\Desktop\WaveModel\Dam_size18mm\CPRESS_all_node_dam18mm150kHz.xlsx');
x_cpp=xlsread('C:\Users\abhij\Desktop\WaveModel\Dam_size18mm\U2CPP_+112_-99_dam18mm150kHz.xlsx');
x_csp=xlsread('C:\Users\abhij\Desktop\WaveModel\Dam_size18mm\U2CSP_-112_+99_dam18mm150kHz.xlsx');
x_sen=xlsread('C:\Users\abhij\Desktop\WaveModel\Dam_size18mm\U2SEN_1-3_-1_dam18mm150kHz.xlsx');
%% Extraction of CPRESS_-112_+99_dam10mm180kHz i.e cpress
x_cpress =zeros(1001,17*102);
for i=1:17
    for j=1:102
        x_cpress(:,(i-1)*102+j)=x_all_node_cpress(:,107+(i-1)*208+j);
    end
end
%% CPRESS extraction and arrangment
X_CPRESS=zeros(1001,238);
X_CPRESS(:,1)=x_cpress(:,1);
X_CPRESS(:,15)=x_cpress(:,103);
X_CPRESS(:,29)=x_cpress(:,205);
X_CPRESS(:,43)=x_cpress(:,307);
X_CPRESS(:,57)=x_cpress(:,409);
X_CPRESS(:,71)=x_cpress(:,511);
X_CPRESS(:,85)=x_cpress(:,613);
X_CPRESS(:,99)=x_cpress(:,715);
X_CPRESS(:,113)=x_cpress(:,817);
X_CPRESS(:,127)=x_cpress(:,919);
X_CPRESS(:,141)=x_cpress(:,1021);
X_CPRESS(:,155)=x_cpress(:,1123);
X_CPRESS(:,169)=x_cpress(:,1225);
X_CPRESS(:,183)=x_cpress(:,1327);
X_CPRESS(:,197)=x_cpress(:,1429);
X_CPRESS(:,211)=x_cpress(:,1531);
X_CPRESS(:,225)=x_cpress(:,1633);
for i=1:17
    for j=1:13
       X_CPRESS(:,(i-1)*14+j+1)=x_cpress(:,(i-1)*102+(j-1)*8+6);
    end
end
%% x_cpp column arrangment
 X_CPP =zeros(1001,239);
 X_CPP(:,1) = x_cpp(:,1);
for i=1:17
    for j= 1:14
        X_CPP(:,(i-1)*14+j+1)=x_cpp(:,(17-i)*14+1+j);
    end
end
%% Calculation of (u2_plat-u2_stiff) #Displacament difference for each time step
for i=1:238
    x_diff=X_CPP(1:end,i+1)-x_csp(1:end,i+1);
    X_Diff(:,i)=x_diff;
end
%% Energy calculation for each time step
for j=1:238
    energy_t = X_CPRESS(1:end,j).*(abs(X_Diff(:,j)));
    ENERGY_t(:,j)=energy_t;
end
%%
figure
t =x_csp(:,1)*1000000;
for j=1:238
    hold on
    plot(t,abs(ENERGY_t(:,j)))
end
xlabel('time (\mu sec)')
ylabel('Cont Energy Intensity')
%% Summation energy
for j=1:238
    enr_sum =sum(ENERGY_t(:,j));
    ENR_SUM(:,j)=enr_sum;
end
%%
figure
plot(ENR_SUM,'b.')
%%
ENR_Matrix = zeros(17,14);
for i=1:17
    for j=1:14
        ENR_Matrix(i,j) = ENR_SUM(:,(i-1)*14+j);
    end
end   
%figure;surf(ENR_Matrix);shading interp
%% Energy distribution over contact area
for i=1:17
    for j=1:14
        enr_mat_1(i,j)=ENR_Matrix(i,j);
    end
end
for i=1:17
    for j=1:13
        enr_mat_2(i,j)=ENR_Matrix(i,14-j);
    end
end
ENR_MAT=[enr_mat_1 enr_mat_2];
figure;surf(ENR_MAT);shading interp
%% Contact Pressure Distribution
% CPress Summation
for j=1:238
    cpress_sum =sum(X_CPRESS(:,j));
    CPRESS_SUM(:,j)=cpress_sum;
end

CPRESS_Matrix =zeros(17,14);
for i=1:17
    for j=1:14
        CPRESS_Matrix(i,j) = CPRESS_SUM(:,(i-1)*14+j);
    end
end   

for i=1:17
    for j=1:14
        cpress_mat_1(i,j)=CPRESS_Matrix(i,j);
    end
end
for i=1:17
    for j=1:13
        cpress_mat_2(i,j)=CPRESS_Matrix(i,14-j);
    end
end
CPRESS_MAT=[cpress_mat_1 cpress_mat_2];
figure;surf(CPRESS_MAT);shading interp
%% Plotting transverse displacement
%// Define the x values
t=x_csp(2:end,1).*1000000;
tMat = repmat(t,1,9); %// For plot3
%// Define y values
y = 1:2:17;
yMat = repmat(y, numel(t), 1); %//For plot3
%// Define z values
p=14;
c_p1 = p+1;c_p2 = p+14*2;c_p3 =p+28*2;c_p4 = p+42*2;c_p5 = p+56*2; c_p6 =p+70*2; c_p7 =p+ 84*2; c_p8 = p+98*2; % contact position
c_p9 = p+112*2; c_p10 = p+126*2;c_p11 = p+140;
zSCP1 = X_CPP(2:end,c_p1);
zSCP2 = X_CPP(2:end,c_p2);
zSCP3 = X_CPP(2:end,c_p3);
zSCP4 = X_CPP(2:end,c_p4);
zSCP5 = X_CPP(2:end,c_p5);
zSCP6 = X_CPP(2:end,c_p6);
zSCP7 = X_CPP(2:end,c_p7);
zSCP8 = X_CPP(2:end,c_p8);
zSCP9 = X_CPP(2:end,c_p9);
%zSCP10 = X_CPP(2:end,c_p10);
%zSCP11 = X_CPP(2:end,c_p11);
zMat1 =[zSCP1 zSCP2 zSCP3 zSCP4 zSCP5 zSCP6 zSCP7 zSCP8 zSCP9];% zSCP10];% zSCP11];
zSCS1 = x_csp(2:end,c_p1);
zSCS2 = x_csp(2:end,c_p2);
zSCS3 = x_csp(2:end,c_p3);
zSCS4 = x_csp(2:end,c_p4);
zSCS5 = x_csp(2:end,c_p5);
zSCS6 = x_csp(2:end,c_p6);
zSCS7 = x_csp(2:end,c_p7);
zSCS8 = x_csp(2:end,c_p8);
zSCS9 = x_csp(2:end,c_p9);
%zSCS10 = x_csp(2:end,c_p10);
%zSCS11 = x_csp(2:end,c_p11);
zMat2 =[zSCS1 zSCS2 zSCS3 zSCS4 zSCS5 zSCS6 zSCS7 zSCS8 zSCS9];% zSCS10];% zSCS11];   %// For plot3
figure
plot3(tMat, yMat, zMat1, 'b','LineWidth',1); %// Make all traces blue
hold on
plot3(tMat, yMat, zMat2, 'g','LineWidth',1)
grid;
xlabel('Time (µ sec)'); %ylabel('Location'); 
zlabel('Transverse Displacement');
view(40,40); %// Adjust viewing angle so you can clearly see data