clear; clf;clc; close all; fclose('all');
addpath([pwd,'/utils']); addpath([pwd,'/output']);
closevar; 

time_index_max=49;

filename='input.xlsx';
    input_solid=readtable(filename).solid';
    input_liquid=readtable(filename).liquid';
    N=readtable(filename).N'; % normalization
    Kd_cpx=readtable(filename).Kd_cpx';
    data=xlsread(filename,'data');

TE_list = split('La Ce Pr Nd Sm Eu Gd Tb Dy Ho Er Tm Yb Lu'); nTE = length(TE_list);
TE_full = split('La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu');

%normalization
solid_N_aux=input_solid./N;
liquid_N_aux=input_liquid./N;
eq_solid_N_aux=Kd_cpx.*input_liquid./N;
data_N_aux=data./N;

%Pm interpolation
Pm_pos=find(contains(TE_full,"Pm"));
solid_N=[solid_N_aux(1:Pm_pos-1) sqrt(solid_N_aux(Pm_pos-1)*solid_N_aux(Pm_pos)) solid_N_aux(Pm_pos:end)] ; 
liquid_N=[liquid_N_aux(1:Pm_pos-1) sqrt(liquid_N_aux(Pm_pos-1)*liquid_N_aux(Pm_pos)) liquid_N_aux(Pm_pos:end)] ;
eq_solid_N=[eq_solid_N_aux(1:Pm_pos-1) sqrt(eq_solid_N_aux(Pm_pos-1)*eq_solid_N_aux(Pm_pos)) eq_solid_N_aux(Pm_pos:end)] ;
data_N=[data_N_aux(:,1:Pm_pos-1) sqrt(data_N_aux(:,Pm_pos-1).*data_N_aux(:,Pm_pos)) data_N_aux(:,Pm_pos:end)] ; 

%% REE patterns (Fig. 1)
figure(1); clf(1)
ymin = 0.1; ymax = 1000;

subplot(1,3,1); % models at 20 kyr & depth 1
depth_index = 1; %out of 100 (from 1-bottom to 100-top)
time_index=time_index_max;
baseFileName = sprintf('savetime_%d', time_index);
load(baseFileName);       

semilogy(0:nTE,solid_N(:),'k'); hold on
semilogy(0:nTE,liquid_N(:),'--r');
semilogy(0:nTE,eq_solid_N(:),'r');

for zoning_index = 1:length(xmesh_Cpx)
    TE_pattern_aux(zoning_index,:) = ME_Cpx_save(depth_index*nTE-nTE+1:depth_index*nTE,zoning_index)./N';
    Pm(zoning_index) = sqrt(TE_pattern_aux(zoning_index,4)*TE_pattern_aux(zoning_index,5)); %Pm interpolation
    TE_pattern(zoning_index,:) = [TE_pattern_aux(zoning_index,1:4) Pm(zoning_index) TE_pattern_aux(zoning_index,5:end)];

    if rem(zoning_index,round(length(xmesh_Cpx)/10))==0
        semilogy(0:nTE,TE_pattern(zoning_index,:),'b');
    end
    hold on
end

xticks([0:nTE]); xticklabels(TE_full)
ylabel(['Normalized REE concentrations'])  
ylim([ymin ymax])

subplot(1,3,2); % models at 5 kyr & depth 3
depth_index = 3; %out of 100 (from 1-bottom to 100-top)
time_index=10;
baseFileName = sprintf('savetime_%d', time_index);
load(baseFileName);

semilogy(0:nTE,solid_N(:),'k'); hold on
semilogy(0:nTE,liquid_N(:),'--r');
semilogy(0:nTE,eq_solid_N(:),'r');

for zoning_index = 1:length(xmesh_Cpx)
    TE_pattern_aux(zoning_index,:) = ME_Cpx_save(depth_index*nTE-nTE+1:depth_index*nTE,zoning_index)./N';
    Pm(zoning_index) = sqrt(TE_pattern_aux(zoning_index,4)*TE_pattern_aux(zoning_index,5)); %Pm interpolation
    TE_pattern(zoning_index,:) = [TE_pattern_aux(zoning_index,1:4) Pm(zoning_index) TE_pattern_aux(zoning_index,5:end)];

    if rem(zoning_index,round(length(xmesh_Cpx)/10))==0
        semilogy(0:nTE,TE_pattern(zoning_index,:),'b');
    end
    hold on
end

xticks([0:nTE]); xticklabels(TE_full)
ylabel(['Normalized REE concentrations'])  
ylim([ymin ymax])

subplot(1,3,3); %data

for i = 1:length(data_N(:,1))
    semilogy(0:nTE,data_N(i,:),'k');
    hold on
end

xticks([0:nTE]); xticklabels(TE_full)
ylabel(['Normalized REE concentrations']) 
ylim([ymin ymax])

%% Covariation (Fig. 2)
figure(2);

data_Eu_anomaly=data_N(:,7)./(sqrt(data_N(:,6).*data_N(:,8)));      
scatter(data(:,6),data_Eu_anomaly(:),10,'filled','k');                                               
hold on

for time_index=1:time_index_max

    baseFileName = sprintf('savetime_%d', time_index);
    load(baseFileName);       

    for zoning_index = 1:length(xmesh_Cpx)
        depth_index = 1;
        Sm_N(time_index,zoning_index) = ME_Cpx_save(depth_index*nTE-nTE+5,zoning_index)./N(5);
        Eu_N(time_index,zoning_index) = ME_Cpx_save(depth_index*nTE-nTE+6,zoning_index)./N(6);
        Gd_N(time_index,zoning_index) = ME_Cpx_save(depth_index*nTE-nTE+7,zoning_index)./N(7);
        Eu_anomaly_N(time_index,zoning_index) = Eu_N(time_index,zoning_index)/(sqrt(Sm_N(time_index,zoning_index)*Gd_N(time_index,zoning_index)));
        color(time_index,zoning_index) = timeMy_save/timeMy_mat_save(end);
    end

    if rem(time_index,round(time_index_max/10,0))==0
        plot(ME_Cpx_save(depth_index*nTE-nTE+6,:), Eu_anomaly_N(time_index,:));
    end

end

xmin = 0.1; xmax = 10;
ymin = 0.4; ymax = 2;
xlim([xmin xmax]);
ylim([ymin ymax]);  
set(gca,'xscale','log')