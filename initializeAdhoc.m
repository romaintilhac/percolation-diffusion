function [benchmark_2Cpx, benchmark_Eu, r_Eu,input_solid,input_liquid, diff_nTE,ME_solid,NR_solid,NR_fluid,NGamma,fix_radi,diff_nTP,Kd_coeff,D_coeff, E_coeff,V_coeff,NRho0,NRho,NRho_tp,NRho_tp0,NV_mat0,NV_mat,NP_mat,NDiv0,NDiv,NPhi0,NPhi,NTp_wt0,NTp_wt,NTp_vol0,NTp_vol,NT,timestep,time_end,save_list,input_TP_v,input_TP_w] = ...
    initializeAdhoc(nTP, nTE,ME_solid,NR_solid,NR_fluid,NGamma,NRho0,NRho,NRho_tp,NRho_tp0,NV_mat0,NV_mat,NP_mat,NDiv0,NDiv,NPhi0,NPhi,NTp_wt0,NTp_wt,NTp_vol, NTp_vol0,NT,ystpt,ysize)

%% Load inputs & initialize
filename='input.xlsx';
    input_solid=readtable(filename).solid';
    input_liquid=readtable(filename).liquid';
    normalization=readtable(filename).N';
    Kd_cpx=readtable(filename).Kd_cpx';
    D_cpx=readtable(filename).D_cpx';
    E_cpx=readtable(filename).E_cpx';
    V_cpx=readtable(filename).V_cpx';

    diff_nTE = ones(1,nTE);
    Kd_coeff=zeros(nTP, nTE);
    D_coeff=zeros(nTP, nTE);
    E_coeff=zeros(nTP, nTE);
    V_coeff=zeros(nTP, nTE);

%% Main parameters

    r_Eu = 0.1; % Eu2+/Eu3+ proportion
    benchmark_Eu=1; % 1 to use synthetic diopside (conservative), 0 to use natural diopside (both Sr diffusivities from Sneeringer et al., 1984)

    fix_poros = 0.01;  % constant melt porosity (vol)
    fix_P = 1.2;  % constant pressure [GPa]
    fix_T = 1200;  % constant temperature [ºC]
    fix_Vy = 5; % constant melt velocity [cm/y]
    fix_Vy = fix_Vy/100/365/24/60/60; % constant melt velocity [m/s]
                      
%     time_critical = ysize/fix_Vy;
%     n_time_critical = 2; % number of critical times
%     time_end=n_time_critical*time_critical; 

    time_end_My = 0.02; % duration [My]
    time_end = time_end_My*1e6*365*24*60*60; % [s]
    
    timestep = ystpt/fix_Vy/5; %check with number of nodes and particle spacing
    %timestep = time_end/10;
    ntimesteps = round(time_end/timestep);

    save_interval = 1; % saving at every timestep
    save_list = 1:save_interval:ntimesteps;
    
%% Mineral compositions

    benchmark_2Cpx =1; % 1 to allocate 2 populations of cpx  % (cpx1: porph & cpx2: matrix)

                                                %cpx2   cpx1        
    nTP_vol =                           [0.5     0.5     0    0    0    0]; % modal proportion [vol]
    diff_nTP  =  double(nTP_vol>0);
    rho_nTP   =  diff_nTP.*     [  1        1       1    1    1    1];
    ini_radi =                            [ 0.1      1       0    0    0    0]*1E-3; % grain size in [mm] (to [m])
    diffT_nTP =  diff_nTP.*     [  1        1       1    1    1    1];  % P dependency of the diffusivity
    diffP_nTP = diffT_nTP.*    [  1        1       0    1    0    0];  % T dependency of the diffusivity

    Kd_coeff(1,:) = Kd_cpx;
    Kd_coeff(2,:) =Kd_cpx;

    D_coeff(1,:) = D_cpx;
    E_coeff(1,:) = E_cpx*diffT_nTP(1);
    V_coeff(1,:) = V_cpx*diffP_nTP(1);
    D_coeff(2,:) = D_cpx;
    E_coeff(2,:) = E_cpx*diffT_nTP(2);
    V_coeff(2,:) = V_cpx*diffP_nTP(2);
    
%% Assign variables

    rho_nPH   = [nTP_vol*rho_nTP' 1];
    Rho_ave   = [(1-fix_poros) fix_poros]*rho_nPH';
    input_TP_v = [nTP_vol/sum(nTP_vol)*(1-fix_poros) fix_poros];
    input_TP_w = input_TP_v.*[rho_nTP rho_nPH(2)]./Rho_ave;

    n_Tp = [5E5 5E5 5E5 5E5 5E5 5E5];
    n_Tp(logical(diff_nTP)) = [nTP_vol(logical(diff_nTP))/sum(nTP_vol)*(1-fix_poros)]./((4/3)*pi.*ini_radi(logical(diff_nTP)).^3);
    ME_solid(:,3+nTE+1:end) = repmat(n_Tp,size(ME_solid,1),1);

    aux_radi=0;
    fix_radi  = aux_radi*diff_nTP.*ini_radi; 
    fix_poros = [1-fix_poros fix_poros];    
    fix_press = [fix_P fix_P];
    fix_velo  = [0 0 0 fix_Vy];

% from initializeDynamics
    NR_solid = 0*NR_solid;
    NR_fluid = 0*NR_fluid;
    NGamma = 0*NGamma;
    NRho0 = ones(size(NRho0));
    NRho =  ones(size(NRho));
    NRho_tp = repmat(rho_nTP,size(NRho_tp,1),1);
    NRho_tp0 = repmat(rho_nTP,size(NRho_tp0,1),1);
    NV_mat0 = repmat(fix_velo,size(NV_mat0,1),1);
    NV_mat = repmat(fix_velo,size(NV_mat,1),1);
    NP_mat = repmat(fix_press,size(NP_mat,1),1);
    NDiv0 = 0*NDiv0;
    NDiv = 0*NDiv;
    NPhi0 = repmat(fix_poros,size(NPhi0,1),1);
    NPhi = repmat(fix_poros,size(NPhi,1),1);
    NTp_wt0 = repmat(input_TP_w,size(NTp_wt0,1),1);
    NTp_wt = repmat(input_TP_w,size(NTp_wt,1),1);
    NTp_vol0 = repmat(input_TP_v,size(NTp_vol0,1),1);
    NTp_vol = repmat(input_TP_v,size(NTp_vol,1),1);
    NT = repmat(fix_T,size(NT,1),1);
