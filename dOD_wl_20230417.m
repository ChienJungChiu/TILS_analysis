%{
input:
PreProcessed2 output:
DeltaOD_All_trail(.mat) or excel output
(x-axis: time; y-axis: wavelength)

output:

Chien-Jung Chiu
Last Update: 2023/4/21
%}

clc; clear all; close all;

%% check the input data
wavelength_selection = [660:5:980];  %9~49
subject = 'DMSBefore_Subject_5_Day1_Round1';
channel = 'Ch2';
DeltaOD_All = load(fullfile('input','TILS-810nm','DeltaOD_MA1',subject,channel,'DMSBefore_Subject_5_Day1_Round1_ALL_DeltaOD.mat'));
DeltaOD_ma_long = DeltaOD_All.DeltaOD_All{1,1};  %load from PreProcessed2 .mat output 
DeltaOD_ma_short = DeltaOD_All.DeltaOD_All{1,2};
%the 2-row (Mark) of long_excel sign symbol difference determinate the change of different stages 
long_excel = xlsread(fullfile('Excel Result','TILS-810nm','DeltaOD_MA1',subject,channel,'DeltaOD.xlsx'),1);  %load from PreProcessed2 output excel file
short_excel = xlsread(fullfile('Excel Result','TILS-810nm','DeltaOD_MA1',subject,channel,'DeltaOD.xlsx'),2);
% compare two different sources data. check whether they are the same.
a = long_excel(4:68,2:750);  
b = short_excel(3:67,2:750);
if (isequal(a,DeltaOD_ma_long) == 0)  %0 is false; 1 is true
    disp('Please check your source data (long)!')
end
if (isequal(b,DeltaOD_ma_short) == 0)
    disp('Please check your source data (short)!')
end
%plot(a(:,130));

%% split DeltaOD_All into 30 trails
k = 1;
trail_first = 1;
for i = 1:size(DeltaOD_ma_long,2)
   if (long_excel(2,i+1)*long_excel(2,i+2)) < 0
       DeltaOD_EachTrail(1,k) = {DeltaOD_ma_long(:,trail_first:i)};
       DeltaOD_EachTrail(2,k) = {DeltaOD_ma_short(:,trail_first:i)};
       k = k +1; 
       trail_first = i + 1;
   end
   if i == size(DeltaOD_ma_long,2)-2
       DeltaOD_EachTrail(1,k) = {DeltaOD_ma_long(:,trail_first:size(DeltaOD_ma_long,2))};
       DeltaOD_EachTrail(2,k) = {DeltaOD_ma_short(:,trail_first:size(DeltaOD_ma_long,2))};
       break;
   end
   
end


%% Delta_OD_EachTrail block average
for trial = 1:30
    %firt shift the delta_OD for each time gate
%     shift_long =  DeltaOD_EachTrail{1,trial+1} - DeltaOD_EachTrail{1,trial+1}(1,:);
%     shift_short = DeltaOD_EachTrail{2,trial+1} - DeltaOD_EachTrail{2,trial+1}(1,:);
    if trial == 1
        DeltaOD_long = [];
        DeltaOD_short = [];
        DeltaOD_long = DeltaOD_EachTrail{1,trial+1};
        DeltaOD_short = DeltaOD_EachTrail{2,trial+1};
    else
        DeltaOD_long = DeltaOD_long + DeltaOD_EachTrail{1,trial+1};
        DeltaOD_short = DeltaOD_short + DeltaOD_EachTrail{2,trial+1};
    end
%       plot(DeltaOD_EachTrail{1,trial+1}(:,10));
%       hold on;
%       plot(DeltaOD_EachTrail{2,trial+1}(:,10));
%       hold on;
%       plot(shift_long(:,10));
%       hold on;
%       plot(shift_short(:,10));
%       legend('long ori','short ori','long shift','short shift')
end
DeltaOD_short_block_average = DeltaOD_short./30;
DeltaOD_long_block_average = DeltaOD_long./30;

%plot delta_OD after block average

% figure;
% plot(DeltaOD_short_block_average(9:49,10));
% %xlabel('Wavelength(nm)');ylabel('Delta OD');
% hold on;
% plot(DeltaOD_long_block_average(9:49,10));
% title('After Block Average');
% xlabel('Wavelength(nm)');ylabel('Delta OD');%xticklabels({'660', num2str(660+wl_ticks_interval), num2str(660+2*wl_ticks_interval), num2str(660+3*wl_ticks_interval), num2str(660+4*wl_ticks_interval), num2str(660+5*wl_ticks_interval), num2str(660+6*wl_ticks_interval), num2str(660+7*wl_ticks_interval)})
% legend('short','long')
% %xticklabels(num2str(xticks*5+660){1})

% %shift
% for wl = 1:65
%     for t = 1:21
%         temp_short(wl,t) = DeltaOD_short_block_average(wl,t) - DeltaOD_short_block_average(9,t);
%         temp_long(wl,t) = DeltaOD_long_block_average(wl,t) - DeltaOD_long_block_average(wl,1);
%     end
% end

DeltaOD_block_average = [DeltaOD_short_block_average; DeltaOD_long_block_average];
%DeltaOD_block_average = DeltaOD_block_average - DeltaOD_block_average(1,:);


%%
wavelength_selection = wavelength_selection';
mean_pathlength = load(fullfile('mean_path','TCThesis','Pathlength_TCThesis_Subject_5.mat'));
mean_pathlength = mean_pathlength.B;
mean_pathlength = interp1(mean_pathlength(:,1), mean_pathlength(:,2:end),wavelength_selection);
mean_pathlength = cat(2,wavelength_selection,mean_pathlength);

molar_extinction_coefficient = load(fullfile('molar_extinction_coefficient','MolarExtinctionCoefficient.mat'));
molar_extinction_coefficient = molar_extinction_coefficient.molar_extinction_coefficient;
HbO2_molar_coefficient =((interp1(molar_extinction_coefficient(:,1), molar_extinction_coefficient(:,2),wavelength_selection)).*2.303); %M
Hb_molar_coefficient =(( interp1(molar_extinction_coefficient(:,1), molar_extinction_coefficient(:,3),wavelength_selection)).*2.303);  %M
cytoxidase_molar_coefficient =(( interp1(molar_extinction_coefficient(:,1), molar_extinction_coefficient(:,4),wavelength_selection)).*2.303);  %M
new_molar_extinction_coefficient=[wavelength_selection HbO2_molar_coefficient Hb_molar_coefficient cytoxidase_molar_coefficient];

approach1_ratio = 0;
sensitivity_matrix_3=[(mean_pathlength(:,2)+approach1_ratio.*mean_pathlength(:,3)).* new_molar_extinction_coefficient(:,2:3)  mean_pathlength(:,5).*new_molar_extinction_coefficient(:,2:4) ; 
                           (mean_pathlength(:,7)+approach1_ratio.*mean_pathlength(:,8)).* new_molar_extinction_coefficient(:,2:3)  mean_pathlength(:,10).*new_molar_extinction_coefficient(:,2:4)]; 
sensitivity_matrix_2=[(mean_pathlength(:,2)+approach1_ratio.*mean_pathlength(:,3)).* new_molar_extinction_coefficient(:,2:3)  mean_pathlength(:,5).*new_molar_extinction_coefficient(:,2:3) ; 
                           (mean_pathlength(:,7)+approach1_ratio.*mean_pathlength(:,8)).* new_molar_extinction_coefficient(:,2:3)  mean_pathlength(:,10).*new_molar_extinction_coefficient(:,2:3)]; 

%%                       
delta_concentration_2=sensitivity_matrix_2\DeltaOD_block_average; %unit: molar
delta_concentration_3=sensitivity_matrix_3\DeltaOD_block_average;

% plot(delta_concentration_2);
% plot(delta_concentration_3);
calculate_deltaOD_2 = sensitivity_matrix_2*delta_concentration_2;
calculate_deltaOD_3 = sensitivity_matrix_3*delta_concentration_3;

t = 10;
% delta_OD_2_3_short = calculate_deltaOD_2(1:65,t)-calculate_deltaOD_3(1:65,t);
% delta_OD_2_3_long = calculate_deltaOD_2(66:130,t)-calculate_deltaOD_3(66:130,t);
% plot(calculate_deltaOD_3(66:130,t))
% plot(gca,calculate_deltaOD_2(1:65,10));
% a = (980-660)/70*xticks+660;
% %set(gca,'xticks',a);
% %xticks(660:5:980)
% hold on;
% 
% %DeltaOD_block_average = DeltaOD_block_average*-1;
% plot(DeltaOD_block_average(1:65,t));

%%
figure;
for t = 1:21
    figure;
    plot(DeltaOD_short_block_average(9:49,t),'r','LineWidth',3);
    hold on;
    plot(calculate_deltaOD_2(9:49,t),'g','LineWidth',3)
    hold on;
    plot(calculate_deltaOD_3(9:49,t),'b','LineWidth',3)
    xlim([1 41]);  %700, 705, 710, 715, ..., 900
    %title('Short channel');
    xlabel('Wavelength(nm)');ylabel('\DeltaOD');xticks([1:10:90]);xticklabels({'700','750','800','850','900'});
    legend('measured','fit__2 chromophore','fit__3 chromophore','Location','northwest','FontSize',16)
    legend('boxoff')
    ax = gca;  
    ax.FontSize = 16;

    figure;
    plot(DeltaOD_long_block_average(9:49,t),'r','LineWidth',3);
    hold on;
    plot(calculate_deltaOD_2(74:114,t),'g','LineWidth',3)
    hold on;
    plot(calculate_deltaOD_3(74:114,t),'b','LineWidth',3)
    hold on;
    xlim([1 41]);  %700, 705, 710, 715, ..., 900
    %title('long channel');
    xlabel('Wavelength(nm)');ylabel('\DeltaOD');xticks([1:10:90]);xticklabels({'700','750','800','850','900'});
    legend('measured','fit__2 chromophore','fit__3 chromophore','Location','northwest','FontSize',16)
    legend('boxoff')
    ax = gca;  
    ax.FontSize = 16;
end

%% residual & residual difference
for wl = 1:65
    for time = 1:21
        residual_long_lo2(wl,time) = DeltaOD_block_average(wl+65,time) - calculate_deltaOD_2(wl+65,time);
        residual_long_lo3(wl,time) = DeltaOD_block_average(wl+65,time) - calculate_deltaOD_3(wl+65,time);
        residual_long_diff(wl,time) = residual_long_lo2(wl,time) - residual_long_lo3(wl,time);
        residual_short_lo2(wl,time) = DeltaOD_block_average(wl,time) - calculate_deltaOD_2(wl,time);
        residual_short_lo3(wl,time) = DeltaOD_block_average(wl,time) - calculate_deltaOD_3(wl,time);
        residual_short_diff(wl,time) = residual_short_lo2(wl,time) - residual_short_lo3(wl,time);
%         plot(residual_lo2(wl,time));
%         hold on;
%         plot(residual_lo3(1:40,time));
%         hold on;
%         plot(residual_diff(1:40,time))
%         legend('2','3','residual difference')
    end
end
 
%plot
for time = 1 :21
    figure;
    %subplot(1,2,1)
    plot(residual_long_lo2(9:49,time),'g','LineWidth',3);
    hold on;
    plot(residual_long_lo3(9:49,time),'b','LineWidth',3);
    %title('Residual');
    xlabel('Wavelength(nm)');ylabel('\DeltaOD');xticks([1:10:90]);xticklabels({'700','750','800','850','900'});
    legend('fit__2 chromophore','fit__3 chromophore','Location','northwest','FontSize',16)
    ax = gca;  
    ax.FontSize = 16;
    %subplot(1,2,1)
    figure;
    %subplot(1,2,2)
    yyaxis left
    plot(residual_short_diff(9:49,time),'LineWidth',3)
    %subplot(1,2,2)
    %title('Residual Difference');
    xlabel('Wavelength(nm)');ylabel('Residual Difference (\DeltaOD)');
    xticks([1:10:90]);xticklabels({'700','750','800','850','900'});
    hold on;
    %legend('residual difference')
    %figure
    yyaxis right
    plot(cytoxidase_molar_coefficient(9:49,:),'LineWidth',3)
    ylabel('molar absorption coefficient (1/M/cm)) ');
    ax = gca;  
    ax.FontSize = 16;
end

%% residual difference to oxCCO concebtraion
% short_sm = sensitivity_matrix_3(9:49,:);
% long_sm = sensitivity_matrix_3(74:114,:);
% new_sensitivity_matrix_3 = [short_sm; long_sm];

GM_HbO_con = delta_concentration_3(3,:) - delta_concentration_3(3,1) ;
GM_Hb_con = delta_concentration_3(4,:) - delta_concentration_3(4,1) ;
GM_oxCCO_con = delta_concentration_3(5,:) - delta_concentration_3(5,1) ;
figure;
plot(GM_HbO_con);
hold on;
plot(GM_Hb_con);
hold on;
plot(GM_oxCCO_con);
title('DMS before ch2');
xlabel('time');ylabel('concentration change');xlim([1 21]);%xticklabels({'1','5','10','850','900'});
legend('HbO', 'Hb', 'oxCCO')

molar absorption coefficient (1/M/cm)) 
residual_diff = [residual_short_diff; residual_long_diff];
oxCCO_con = sensitivity_matrix_3\residual_diff;
plot(DeltaOD_long_block_average(9:49,t));
title('After Block Average');
xlabel('Wavelength(nm)');ylabel('Delta OD');%xticklabels({'660', num2str(660+wl_ticks_interval), num2str(660+2*wl_ticks_interval), num2str(660+3*wl_ticks_interval), num2str(660+4*wl_ticks_interval), num2str(660+5*wl_ticks_interval), num2str(660+6*wl_ticks_interval), num2str(660+7*wl_ticks_interval)})
legend('short','long')