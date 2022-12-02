%% MAIN PERCOLATION-DIFFUSION
    % Computes the diffusional re-equilibration of REE in a solid matrix
    % of spherical mantle minerals percolated by a melt in a 1D column

    % Modified from the original MPMCRT code developped by Beñat Oliveira Bravo
    % By Romain Tilhac, Decembre 2022
    % Contact: romain.tilhac@csic.es

    clear; clf;clc; close all;
    addpath([pwd,'/utils']);
    closevar; 
    
%% PROBLEM SETUP

    % PROBLEM STATEMENT 
        x1 = 0;                     % Domain X first vertice 
        x2 = 300;                 % Domain X second vertice        
        xsize = x2 - x1;
        
        y1 = 0;                     % Domain Y first vertice
        y2 = 1000;              % Domain Y second vertice
        ysize = y2 - y1;

    % INFO EULERIAN MESH
        % Temperature
            elemTypeT = 1;               % Quads

            nent = 4;                    % Number of nodes per element (velocity)
            nxt = 10;                    % Number of nodes -1 in x direction
            nyt = 10;                    % Number of nodes -1 in y direction

            nnt = (nxt+1)*(nyt+1);       % Number of total temperature nodes

            xstpt = xsize/nxt;           % Horizontal grid step
            ystpt = ysize/nyt;           % Vertical grid step
            
        % Create the meshes for temperature
            [XT,TT]=createVelocityMesh(elemTypeT,nent,x1,x2,y1,y2,nxt,nyt);

        % Making vectors for nodal points positions (basic nodes)
            xsize_aux = 2000;
            ysize_aux = 5000;

                gridyt   =   (0:ystpt:ysize);   % Vertical
                gridxt   =   (0:xstpt:xsize);   % Horizontal

    % INFO LAGRANGIAN PARTICLES
        % Markers
            mxelem  = 5;                    % Markers in x direction per element
            myelem  = 5;                    % Markers in y direction per element

        % Moving Markers: 
            markmove=4; % 0 = not moving at all, 1 = simple 1-st order advection, 4 = 4-th order in space Runge-Kutta

        % Create Particles
            % Preliminary info
                TE_list = split("La Ce Pr Nd Sm Eu Gd Tb Dy Ho Er Tm Yb Lu");
                Eu_pos=find(contains(TE_list,"Eu")); %position of Eu in the list of trace elements
                nTE = length(TE_list); % trace elements
                nTP = 6; % thermodynamic phases

                [ME_solid, ME_fluid]  = createParticles_Elements(mxelem,myelem,nxt,nyt,xsize,ysize,nTE,nTP);
                
        nc = 51; % number of nodes for the diffusion profiles

        tol_NTE = 1e-8;
        tol_TE = 1E-9;
        tol_TE_abs = 5E-9;

        R = 8.3144598;    % [J/(mol K)]
        TK=273.15; PGPa=1E9; 

%% INITIALIZATION

    [NR_solid,NR_fluid,NGamma,NRho0,NRho,NRho_tp,NRho_tp0,NV_mat0,NV_mat,NP_mat,NDiv0,NDiv,NPhi0,NPhi,NTp_wt0,NTp_wt,NTp_vol0,NTp_vol,NT] = ...
    initializeDynamics('default.mat',nxt,nyt); 

    [benchmark_2Cpx, benchmark_Eu, r_Eu,TE_input_solid,TE_input_liquid, diff_nTE, ME_solid,NR_solid,NR_fluid,NGamma,fix_radi,diff_nTP,Kd_coeff,D_coeff, E_coeff,V_coeff,NRho0,NRho,NRho_tp,NRho_tp0,NV_mat0,NV_mat,NP_mat,NDiv0,NDiv,NPhi0,NPhi,NTp_wt0,NTp_wt,NTp_vol0,NTp_vol,NT,timestep,time_end,save_list,input_TP_v,input_TP_w] = ...
    initializeAdhoc(nTP, nTE ,ME_solid,NR_solid,NR_fluid,NGamma,NRho0,NRho,NRho_tp,NRho_tp0,NV_mat0,NV_mat,NP_mat,NDiv0,NDiv,NPhi0,NPhi,NTp_wt0,NTp_vol,NTp_vol0,NTp_wt,NT,ystpt,ysize);

    if size(Kd_coeff,2) ~= nTE
        error("Kd_coeff does not match nTE")
    end
    
    disp(['Initialization completed - ' num2str(time_end/timestep) ' timesteps'])
    keyboard

%% MAIN TIME LOOP
        
        iter = 1;  save_iter = 0;
        time = 0; timeMy_mat = [];

        while time<time_end
            
            disp(['=== TIMESTEP ' num2str(iter) ' ==='])
            
            if iter==1     
                %% FIRST TIMESTEP - CREATE VARIABLES
    
                    % TRACE ELEMENTS
                        % Compute initial profiles inside
                            ME_Oli    = zeros(size(ME_solid,1),nTE);
                            ME_Cpx   = zeros(size(ME_solid,1),nTE);
                            ME_Opx   = zeros(size(ME_solid,1),nTE);
                            ME_Grt    = zeros(size(ME_solid,1),nTE);
                            ME_Spl    = zeros(size(ME_solid,1),nTE);
                            ME_Plg    = zeros(size(ME_solid,1),nTE);
                            ME_melt  = zeros(size(ME_solid,1),nTE);
    
                            ME_solid_Rho_TP = zeros(size(ME_solid,1),nTP+1);
                            ME_solid_Rho_TP(:,1)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,1),ME_solid(:,2));
                            ME_solid_Rho_TP(:,2)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,2),ME_solid(:,2));
                            ME_solid_Rho_TP(:,3)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,3),ME_solid(:,2));
                            ME_solid_Rho_TP(:,4)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,4),ME_solid(:,2));
                            ME_solid_Rho_TP(:,5)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,5),ME_solid(:,2));
                            ME_solid_Rho_TP(:,6)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,6),ME_solid(:,2));
                            NRho0(NRho0(:,2)==0,2)  = NRho0(NRho0(:,2)==0,1);
                            ME_solid_Rho_TP(:,7)  = interp1(XT(1:nxt+1:end,2),NRho0(1:nxt+1:end,2),ME_solid(:,2));
    
                            ME_fluid_Rho_TP = zeros(size(ME_fluid,1),nTP+1);
                            ME_fluid_Rho_TP(:,1)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,1),ME_fluid(:,2));
                            ME_fluid_Rho_TP(:,2)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,2),ME_fluid(:,2));
                            ME_fluid_Rho_TP(:,3)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,3),ME_fluid(:,2));
                            ME_fluid_Rho_TP(:,4)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,4),ME_fluid(:,2));
                            ME_fluid_Rho_TP(:,5)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,5),ME_fluid(:,2));
                            ME_fluid_Rho_TP(:,6)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,6),ME_fluid(:,2));
                            NRho0(NRho0(:,2)==0,2)  = NRho0(NRho0(:,2)==0,1);
                            ME_fluid_Rho_TP(:,7)  = interp1(XT(1:nxt+1:end,2),NRho0(1:nxt+1:end,2),ME_fluid(:,2));
    
                            ME_solid_TP_w0 = zeros(size(ME_solid,1),nTP+1);
    
                            ME_solid_TP_w0(:,1)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,1),ME_solid(:,2));
                            ME_solid_TP_w0(:,2)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,2),ME_solid(:,2));
                            ME_solid_TP_w0(:,3)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,3),ME_solid(:,2));
                            ME_solid_TP_w0(:,4)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,4),ME_solid(:,2));
                            ME_solid_TP_w0(:,5)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,5),ME_solid(:,2));
                            ME_solid_TP_w0(:,6)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,6),ME_solid(:,2));
                            ME_solid_TP_w0(:,7)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,7),ME_solid(:,2));
                            ME_solid_TP_w0       = ME_solid_TP_w0./repmat(sum(ME_solid_TP_w0,2),1,size(ME_solid_TP_w0,2));
    
                            ME_fluid_TP_w0 = zeros(size(ME_fluid,1),nTP+1);
                            ME_fluid_TP_w0(:,1)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,1),ME_fluid(:,2));
                            ME_fluid_TP_w0(:,2)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,2),ME_fluid(:,2));
                            ME_fluid_TP_w0(:,3)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,3),ME_fluid(:,2));
                            ME_fluid_TP_w0(:,4)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,4),ME_fluid(:,2));
                            ME_fluid_TP_w0(:,5)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,5),ME_fluid(:,2));
                            ME_fluid_TP_w0(:,6)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,6),ME_fluid(:,2));
                            ME_fluid_TP_w0(:,7)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,7),ME_fluid(:,2));
                            ME_fluid_TP_w0        = ME_fluid_TP_w0./repmat(sum(ME_fluid_TP_w0,2),1,size(ME_fluid_TP_w0,2));
    
                            ME_solid_TP_v0 = zeros(size(ME_solid,1),nTP+1);
    
                            ME_solid_TP_v0(:,1)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,1),ME_solid(:,2));
                            ME_solid_TP_v0(:,2)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,2),ME_solid(:,2));
                            ME_solid_TP_v0(:,3)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,3),ME_solid(:,2));
                            ME_solid_TP_v0(:,4)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,4),ME_solid(:,2));
                            ME_solid_TP_v0(:,5)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,5),ME_solid(:,2));
                            ME_solid_TP_v0(:,6)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,6),ME_solid(:,2));
                            ME_solid_TP_v0(:,7)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,7),ME_solid(:,2));
                            ME_solid_TP_v0        = ME_solid_TP_v0./repmat(sum(ME_solid_TP_v0,2),1,size(ME_solid_TP_v0,2));
    
                            ME_fluid_TP_v0 = zeros(size(ME_fluid,1),nTP+1);
                            ME_fluid_TP_v0(:,1)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,1),ME_fluid(:,2));
                            ME_fluid_TP_v0(:,2)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,2),ME_fluid(:,2));
                            ME_fluid_TP_v0(:,3)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,3),ME_fluid(:,2));
                            ME_fluid_TP_v0(:,4)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,4),ME_fluid(:,2));
                            ME_fluid_TP_v0(:,5)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,5),ME_fluid(:,2));
                            ME_fluid_TP_v0(:,6)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,6),ME_fluid(:,2));
                            ME_fluid_TP_v0(:,7)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,7),ME_fluid(:,2));
                            ME_fluid_TP_v0        = ME_fluid_TP_v0./repmat(sum(ME_fluid_TP_v0,2),1,size(ME_fluid_TP_v0,2));
    
                        % Partition coefficients - Rows # particles, columns TE
                            Kd_Oli = repmat(Kd_coeff(1,:),size(ME_Oli,1),1);
                            Kd_Cpx = repmat(Kd_coeff(2,:),size(ME_Cpx,1),1);
                            Kd_Opx = repmat(Kd_coeff(3,:),size(ME_Opx,1),1);
                            Kd_Grt = repmat(Kd_coeff(4,:),size(ME_Grt,1),1);
                            Kd_Spl = repmat(Kd_coeff(5,:),size(ME_Spl,1),1);
                            Kd_Plg = repmat(Kd_coeff(6,:),size(ME_Plg,1),1);                

                       if benchmark_2Cpx ==1; Kd_Oli=Kd_Cpx; end %WARNING allocating porphyroclastic matrix cpx to opx

                        % Initial TE (ppm) concentration in TP and Melt 
                            ME_TE    = repmat(TE_input_solid,size(ME_Spl,1),1);
                            auxE_Oli   = repmat(ME_solid_TP_v0(:,1)>0,1,nTE);
                            auxE_Cpx  = repmat(ME_solid_TP_v0(:,2)>0,1,nTE);
                            auxE_Opx  = repmat(ME_solid_TP_v0(:,3)>0,1,nTE);
                            auxE_Grt   = repmat(ME_solid_TP_v0(:,4)>0,1,nTE);
                            auxE_Spl   = repmat(ME_solid_TP_v0(:,5)>0,1,nTE);
                            auxE_Plg   = repmat(ME_solid_TP_v0(:,6)>0,1,nTE);
                            auxE_Melt = repmat(ME_solid_TP_v0(:,7)>0,1,nTE);
    
                            ME_solid_w0_Oli   = repmat(ME_solid_TP_w0(:,1),1,nTE);
                            ME_solid_w0_Cpx  = repmat(ME_solid_TP_w0(:,2),1,nTE);
                            ME_solid_w0_Opx  = repmat(ME_solid_TP_w0(:,3),1,nTE);
                            ME_solid_w0_Grt   = repmat(ME_solid_TP_w0(:,4),1,nTE);
                            ME_solid_w0_Spl   = repmat(ME_solid_TP_w0(:,5),1,nTE);
                            ME_solid_w0_Plg   = repmat(ME_solid_TP_w0(:,6),1,nTE);
                            ME_solid_w0_Melt = repmat(ME_solid_TP_w0(:,7),1,nTE);
    
                        % Input TE is the solid
                            ME_solid(:,4:3+nTE) = ME_TE;
    
                            ME_solid_w0_solid = ME_solid_w0_Oli + ME_solid_w0_Cpx + ME_solid_w0_Opx + ME_solid_w0_Grt + ME_solid_w0_Spl + ME_solid_w0_Plg;  % normalize to 1 every TP contribution
                            ME_solid_w0_total = ME_solid_w0_Oli + ME_solid_w0_Cpx + ME_solid_w0_Opx + ME_solid_w0_Grt + ME_solid_w0_Spl + ME_solid_w0_Plg + ME_solid_w0_Melt;  % normalize to 1 every TP contribution
                            ME_melt_aux = ME_solid(:,4:3+nTE)./(auxE_Oli.*Kd_Oli.*ME_solid_w0_Oli./ME_solid_w0_solid    ...
                                                + auxE_Cpx.*Kd_Cpx.*ME_solid_w0_Cpx./ME_solid_w0_solid ...
                                                + auxE_Opx.*Kd_Opx.*ME_solid_w0_Opx./ME_solid_w0_solid ...
                                                + auxE_Grt.*Kd_Grt.*ME_solid_w0_Grt./ME_solid_w0_solid    ...
                                                + auxE_Spl.*Kd_Spl.*ME_solid_w0_Spl./ME_solid_w0_solid    ...
                                                + auxE_Plg.*Kd_Plg.*ME_solid_w0_Plg./ME_solid_w0_solid);    % [melt] = [solid]/Kd_bulk    (normalized to solid)
    
                            ME_melt(ME_solid_TP_v0(:,7)>0,:) = ME_melt_aux(ME_solid_TP_v0(:,7)>0,:);  % TE mass / fluid mass 
                            ME_melt_eq = ME_melt_aux;
    
                            for TE_pos = 1:nTE
                                ME_fluid(:,3+TE_pos) = interp1(ME_solid(:,2),ME_melt(:,TE_pos),ME_fluid(:,2),'nearest','extrap');
                            end
                           
                            % TE mass / TP mass
                                ME_Oli_eq   =  ME_melt_aux.*auxE_Oli.*Kd_Oli;           ME_Oli(ME_solid_TP_v0(:,1)>0,:)  = ME_Oli_eq(ME_solid_TP_v0(:,1)>0,:);
                                ME_Cpx_eq  =  ME_melt_aux.*auxE_Cpx.*Kd_Cpx ;        ME_Cpx(ME_solid_TP_v0(:,2)>0,:) = ME_Cpx_eq(ME_solid_TP_v0(:,2)>0,:);
                                ME_Opx_eq  =  ME_melt_aux.*auxE_Opx.*Kd_Opx;         ME_Opx(ME_solid_TP_v0(:,3)>0,:) = ME_Opx_eq(ME_solid_TP_v0(:,3)>0,:);
                                ME_Grt_eq   =  ME_melt_aux.*auxE_Grt.*Kd_Grt;           ME_Grt(ME_solid_TP_v0(:,4)>0,:)  = ME_Grt_eq(ME_solid_TP_v0(:,4)>0,:);
                                ME_Spl_eq   =  ME_melt_aux.*auxE_Spl.*Kd_Spl;           ME_Spl(ME_solid_TP_v0(:,5)>0,:)  = ME_Spl_eq(ME_solid_TP_v0(:,5)>0,:);
                                ME_Plg_eq   =  ME_melt_aux.*auxE_Plg.*Kd_Plg;           ME_Plg(ME_solid_TP_v0(:,6)>0,:)  = ME_Plg_eq(ME_solid_TP_v0(:,6)>0,:);
    
                        % Prepare BC for next iteration
                            MBC_Oli_0  = ME_Oli;
                            MBC_Cpx_0 = ME_Cpx;
                            MBC_Opx_0 = ME_Opx;
                            MBC_Grt_0  = ME_Grt;
                            MBC_Spl_0  = ME_Spl;
                            MBC_Plg_0  = ME_Plg;
    
                            MBC_Oli_0(ME_solid_TP_v0(:,7)>0,:)  = ME_melt(ME_solid_TP_v0(:,7)>0,:).*Kd_Oli(ME_solid_TP_v0(:,7)>0,:);
                            MBC_Cpx_0(ME_solid_TP_v0(:,7)>0,:) = ME_melt(ME_solid_TP_v0(:,7)>0,:).*Kd_Cpx(ME_solid_TP_v0(:,7)>0,:);
                            MBC_Opx_0(ME_solid_TP_v0(:,7)>0,:) = ME_melt(ME_solid_TP_v0(:,7)>0,:).*Kd_Opx(ME_solid_TP_v0(:,7)>0,:);
                            MBC_Grt_0(ME_solid_TP_v0(:,7)>0,:)  = ME_melt(ME_solid_TP_v0(:,7)>0,:).*Kd_Grt(ME_solid_TP_v0(:,7)>0,:);
                            MBC_Spl_0(ME_solid_TP_v0(:,7)>0,:)  = ME_melt(ME_solid_TP_v0(:,7)>0,:).*Kd_Spl(ME_solid_TP_v0(:,7)>0,:);
                            MBC_Plg_0(ME_solid_TP_v0(:,7)>0,:)  = ME_melt(ME_solid_TP_v0(:,7)>0,:).*Kd_Plg(ME_solid_TP_v0(:,7)>0,:);
    
                            MBC_Oli  = MBC_Oli_0;
                            MBC_Cpx = MBC_Cpx_0;
                            MBC_Opx = MBC_Opx_0;
                            MBC_Grt  = MBC_Grt_0;
                            MBC_Spl  = MBC_Spl_0;
                            MBC_Plg  = MBC_Plg_0;
    
                        % TE compositions per mineral 
                                for TE_pos = 1:nTE
                                    TE = TE_list(TE_pos);
                                    TEprofile_Oli.(TE) = repmat(ME_Oli(:,TE_pos)',nc,1)';
                                    TEprofile_Cpx.(TE) = repmat(ME_Cpx(:,TE_pos)',nc,1)';
                                    TEprofile_Opx.(TE) = repmat(ME_Opx(:,TE_pos)',nc,1)';
                                    TEprofile_Grt.(TE) = repmat(ME_Grt(:,TE_pos)',nc,1)';
                                    TEprofile_Spl.(TE) = repmat(ME_Spl(:,TE_pos)',nc,1)';
                                    TEprofile_Plg.(TE) = repmat(ME_Plg(:,TE_pos)',nc,1)';
                                    if TE_pos == 1
                                        Mat_Oli_aux = TEprofile_Oli.(TE);
                                        Mat_Cpx_aux = TEprofile_Cpx.(TE);
                                        Mat_Opx_aux = TEprofile_Opx.(TE);
                                        Mat_Grt_aux = TEprofile_Grt.(TE);
                                        Mat_Spl_aux = TEprofile_Spl.(TE);
                                        Mat_Plg_aux = TEprofile_Plg.(TE);
                                    else
                                        Mat_Oli_aux = cat(2,Mat_Oli_aux,TEprofile_Oli.(TE));
                                        Mat_Cpx_aux = cat(2,Mat_Cpx_aux,TEprofile_Cpx.(TE));
                                        Mat_Opx_aux = cat(2,Mat_Opx_aux,TEprofile_Opx.(TE));
                                        Mat_Grt_aux = cat(2,Mat_Grt_aux,TEprofile_Grt.(TE));
                                        Mat_Spl_aux = cat(2,Mat_Spl_aux,TEprofile_Spl.(TE));
                                        Mat_Plg_aux = cat(2,Mat_Plg_aux,TEprofile_Plg.(TE));
                                    end
                                end
    
                            Mat_Oli_order  = reshape(Mat_Oli_aux',nc,[])';   ME_Oli  = mat2cell(Mat_Oli_order,nTE*ones(1,size(ME_solid,1)),nc);                
                            Mat_Cpx_order = reshape(Mat_Cpx_aux',nc,[])';  ME_Cpx = mat2cell(Mat_Cpx_order,nTE*ones(1,size(ME_solid,1)),nc);
                            Mat_Opx_order = reshape(Mat_Opx_aux',nc,[])';  ME_Opx = mat2cell(Mat_Opx_order,nTE*ones(1,size(ME_solid,1)),nc);
                            Mat_Grt_order  = reshape(Mat_Grt_aux',nc,[])';   ME_Grt  = mat2cell(Mat_Grt_order,nTE*ones(1,size(ME_solid,1)),nc); 
                            Mat_Spl_order  = reshape(Mat_Spl_aux',nc,[])';   ME_Spl  = mat2cell(Mat_Spl_order,nTE*ones(1,size(ME_solid,1)),nc); 
                            Mat_Plg_order  = reshape(Mat_Plg_aux',nc,[])';   ME_Plg  = mat2cell(Mat_Plg_order,nTE*ones(1,size(ME_solid,1)),nc);
    
                        % Diffusion
    
                            % P-T in particles;
                                ME_T = interp1(XT(1:nxt+1:end,2),NT(1:nxt+1:end,1),ME_solid(:,2))+TK;
                                ME_P = interp1(XT(1:nxt+1:end,2),NP_mat(1:nxt+1:end,1),ME_solid(:,2))*PGPa;
    
                            % Coefficients
                                D_Oli  = repmat(D_coeff(1,:),size(ME_T,1),1);   E_Oli  = repmat(E_coeff(1,:),size(ME_T,1),1);   V_Oli = repmat(V_coeff(1,:),size(ME_T,1),1);
                                D_Cpx = repmat(D_coeff(2,:),size(ME_T,1),1);  E_Cpx = repmat(E_coeff(2,:),size(ME_T,1),1);  V_Cpx = repmat(V_coeff(2,:),size(ME_T,1),1);
                                D_Opx = repmat(D_coeff(3,:),size(ME_T,1),1);  E_Opx = repmat(E_coeff(3,:),size(ME_T,1),1);  V_Opx = repmat(V_coeff(3,:),size(ME_T,1),1);
                                D_Grt  = repmat(D_coeff(4,:),size(ME_T,1),1);   E_Grt = repmat(E_coeff(4,:),size(ME_T,1),1);    V_Grt = repmat(V_coeff(4,:),size(ME_T,1),1);
                                D_Spl  = repmat(D_coeff(5,:),size(ME_T,1),1);   E_Spl  = repmat(E_coeff(5,:),size(ME_T,1),1);   V_Spl = repmat(V_coeff(5,:),size(ME_T,1),1);
                                D_Plg  = repmat(D_coeff(6,:),size(ME_T,1),1);   E_Plg  = repmat(E_coeff(6,:),size(ME_T,1),1);   V_Plg = repmat(V_coeff(6,:),size(ME_T,1),1);
    
                            % Diffusivites
                                Mdiff_Oli    = double(D_Oli.*exp((-E_Oli + repmat(ME_P,1,nTE).*V_Oli)./(R.*repmat(ME_T,1,nTE))));
                                Mdiff_Cpx   = double(D_Cpx.*exp((-E_Cpx + repmat(ME_P,1,nTE).*V_Cpx)./(R.*repmat(ME_T,1,nTE))));
                                Mdiff_Opx   = double(D_Opx.*exp((-E_Opx + repmat(ME_P,1,nTE).*V_Opx)./(R.*repmat(ME_T,1,nTE))));
                                Mdiff_Grt    = double(D_Grt.*exp((-E_Grt + repmat(ME_P,1,nTE).*V_Grt)./(R.*repmat(ME_T,1,nTE))));
                                Mdiff_Spl    = double(D_Spl.*exp((-E_Spl + repmat(ME_P,1,nTE).*V_Spl)./(R.*repmat(ME_T,1,nTE))));
                                Mdiff_Plg    = double(D_Plg.*exp((-E_Plg + repmat(ME_P,1,nTE).*V_Plg)./(R.*repmat(ME_T,1,nTE))));
    
                            % Eu diffusivities in Cpx based on Sr diffusivities from Sneeringer et al. (1984)
                                if benchmark_Eu==1 %synthetic diopside
                                    D_Cpx_Eu2_ref = 1200*1e-4; % [m^2/s]
                                    E_Cpx_Eu2_ref = 122*4.184*1e3; % [J/mol]
                                elseif benchmark_Eu==0 %natural diopside
                                    D_Cpx_Eu2_ref = 54*1e-4; % [m^2/s]
                                    E_Cpx_Eu2_ref = 97*4.184*1e3; % [J/mol]
                                end
                                D_Cpx_Eu2 = repmat(D_Cpx_Eu2_ref,size(ME_T,1),1); 
                                E_Cpx_Eu2 = repmat(E_Cpx_Eu2_ref,size(ME_T,1),1);
                                Sn84a=repmat(2.94,size(ME_T,1),1); 
                                Sn84b=repmat(4640,size(ME_T,1),1);
                                Mdiff_Cpx_Eu2   = double(D_Cpx_Eu2.*exp(-E_Cpx_Eu2./(R.*ME_T)).*exp(-(ME_P.*1e-9)./((Sn84a.*ME_T-Sn84b).*R))); % [m^2/s]   
                                Mdiff_Cpx_Eu = Mdiff_Cpx_Eu2.*r_Eu + Mdiff_Cpx(:,Eu_pos).*(1-r_Eu); 
                                
                                Mdiff_Cpx = [Mdiff_Cpx(:,1:Eu_pos-1) Mdiff_Cpx_Eu2 Mdiff_Cpx(:,Eu_pos+1:nTE)];
                                if benchmark_2Cpx ==1; Mdiff_Oli=Mdiff_Cpx; end %WARNING allocating porphyroclastic matrix cpx to opx      
    
                        % Compute radii - check and correct n_tp if the radii is fixed
                            Rradi_Oli_0  = ((ME_solid_TP_v0(:,1))./((4/3)*pi.*ME_solid(:,3+nTE+1))).^(1/3);  if diff_nTP(1)==0;  Rradi_Oli_0  = 0*Rradi_Oli_0; end;   if fix_radi(1)~=0; Rradi_Oli_0   = fix_radi(1)*ones(size(Rradi_Oli_0));  ME_solid(:,3+nTE+1) = ME_solid_TP_v0(:,1)./((4/3)*pi.*Rradi_Oli_0.^3); end
                            Rradi_Cpx_0 = ((ME_solid_TP_v0(:,2))./((4/3)*pi.*ME_solid(:,3+nTE+2))).^(1/3);  if diff_nTP(2)==0;  Rradi_Cpx_0 = 0*Rradi_Cpx_0; end;  if fix_radi(2)~=0; Rradi_Cpx_0  = fix_radi(2)*ones(size(Rradi_Cpx_0)); ME_solid(:,3+nTE+2) = ME_solid_TP_v0(:,2)./((4/3)*pi.*Rradi_Cpx_0.^3); end
                            Rradi_Opx_0 = ((ME_solid_TP_v0(:,3))./((4/3)*pi.*ME_solid(:,3+nTE+3))).^(1/3);  if diff_nTP(3)==0;  Rradi_Opx_0 = 0*Rradi_Opx_0; end;  if fix_radi(3)~=0; Rradi_Opx_0  = fix_radi(3)*ones(size(Rradi_Opx_0)); ME_solid(:,3+nTE+3) = ME_solid_TP_v0(:,3)./((4/3)*pi.*Rradi_Opx_0.^3); end
                            Rradi_Grt_0  = ((ME_solid_TP_v0(:,4))./((4/3)*pi.*ME_solid(:,3+nTE+4))).^(1/3);  if diff_nTP(4)==0;  Rradi_Grt_0  = 0*Rradi_Grt_0; end;   if fix_radi(4)~=0; Rradi_Grt_0   = fix_radi(4)*ones(size(Rradi_Grt_0));  ME_solid(:,3+nTE+4) = ME_solid_TP_v0(:,4)./((4/3)*pi.*Rradi_Grt_0.^3); end
                            Rradi_Spl_0  = ((ME_solid_TP_v0(:,5))./((4/3)*pi.*ME_solid(:,3+nTE+5))).^(1/3);  if diff_nTP(5)==0;  Rradi_Spl_0  = 0*Rradi_Spl_0; end;   if fix_radi(5)~=0; Rradi_Spl_0   = fix_radi(5)*ones(size(Rradi_Spl_0));  ME_solid(:,3+nTE+5) = ME_solid_TP_v0(:,5)./((4/3)*pi.*Rradi_Spl_0.^3); end
                            Rradi_Plg_0  = ((ME_solid_TP_v0(:,6))./((4/3)*pi.*ME_solid(:,3+nTE+6))).^(1/3);  if diff_nTP(6)==0;  Rradi_Plg_0  = 0*Rradi_Plg_0; end;   if fix_radi(6)~=0; Rradi_Plg_0   = fix_radi(6)*ones(size(Rradi_Plg_0));  ME_solid(:,3+nTE+6) = ME_solid_TP_v0(:,6)./((4/3)*pi.*Rradi_Plg_0.^3); end
       
                            ME_solid_Rho_TP0 = ME_solid_Rho_TP;
    
                        % ME_TP initial (per mineral)
                            ME_Oli_0  = ME_Oli;
                            ME_Cpx_0 = ME_Cpx;
                            ME_Opx_0 = ME_Opx;
                            ME_Grt_0  = ME_Grt;
                            ME_Spl_0  = ME_Spl;
                            ME_Plg_0  = ME_Plg;
                            
                            tic
                                % This function will not change anything at this stage
                                correct_TE_ind = ones(size(MBC_Oli))>0;
                                correct_TE_met = -ones(size(MBC_Oli))>0;
                                [Gamma_Oli,Gamma_Cpx,Gamma_Opx,Gamma_Grt,Gamma_Spl,Gamma_Plg,ME_Oli,ME_Cpx,ME_Opx,ME_Grt,ME_Spl,ME_Plg,TE_Oli,TE_Cpx,TE_Opx,TE_Grt,TE_Spl,TE_Plg,MRadi,MRadi_0,ME_solid] = ...
                                    diffusion1D(ME_Oli_eq,ME_Cpx_eq,ME_Opx_eq,ME_Grt_eq,ME_Spl_eq,ME_Plg_eq,ME_melt_eq,MBC_Oli_0,MBC_Cpx_0,MBC_Opx_0,MBC_Grt_0,MBC_Spl_0,MBC_Plg_0,MBC_Oli,MBC_Cpx,MBC_Opx,MBC_Grt,MBC_Spl,MBC_Plg,ME_solid,ME_solid,ME_solid_TP_v0,ME_solid_TP_v0,ME_solid_TP_w0,ME_Oli_0,ME_Cpx_0,ME_Opx_0,ME_Grt_0,ME_Spl_0,ME_Plg_0,ME_solid_Rho_TP,ME_solid_Rho_TP0,Mdiff_Oli,Mdiff_Cpx,Mdiff_Opx,Mdiff_Grt,Mdiff_Spl,Mdiff_Plg,timestep,nc,nTE,diff_nTP,diff_nTE,fix_radi,correct_TE_ind,correct_TE_met);
                            toc
    
                            MGamma_e_solid = Gamma_Oli+Gamma_Cpx+Gamma_Opx+Gamma_Grt+Gamma_Spl+Gamma_Plg;
                            MGamma_e_fluid = zeros(size(ME_fluid,1),nTE);
    
                            ME_fluid_solid = zeros(size(ME_solid,1),3+nTE);
                            ME_fluid_solid(:,1:3) = ME_solid(:,1:3);
    
                            for TE_pos = 4:3+nTE
                                    ME_fluid_solid(:,TE_pos) = interp1(ME_fluid(:,2),ME_fluid(:,TE_pos),ME_solid(:,2),'nearest','extrap');
                            end 
    
                            ME_solid_TP_v = ME_solid_TP_v0;
                            ME_solid_TP_w = ME_solid_TP_w0;
                            ME_fluid_TP_v = ME_fluid_TP_v0;
                            ME_fluid_TP_w = ME_fluid_TP_w0;
    
                            Error_Oli = 0*MBC_Oli_0;
                            Error_Cpx = 0*MBC_Cpx_0;
                            Error_Opx = 0*MBC_Opx_0;
                            Error_Grt = 0*MBC_Grt_0;
                            Error_Spl = 0*MBC_Spl_0;
                            Error_Plg = 0*MBC_Plg_0;
    
                            iter_TE = 0; 
    
                else
                %% SECOND TIMESTEP ONWARDS
    
                    % OPENING
    
                        % Recover
                            ME_solid0 = ME_solid;
                            ME_fluid0 = ME_fluid;
                            ME_fluid_solid0 = ME_fluid_solid;
    
                            ME_solid_Rho_TP0 = ME_solid_Rho_TP;
                            ME_Rho_solid0 = ME_Rho_solid;
                            ME_Rho_fluid0 = ME_Rho_fluid;
                            ME_Oli_0  = ME_Oli;
                            ME_Cpx_0 = ME_Cpx;
                            ME_Opx_0 = ME_Opx;
                            ME_Grt_0  = ME_Grt;
                            ME_Spl_0  = ME_Spl;
                            ME_Plg_0  = ME_Plg;
    
                            TE_Oli_0  = TE_Oli;
                            TE_Cpx_0 = TE_Cpx;
                            TE_Opx_0 = TE_Opx;
                            TE_Grt_0  = TE_Grt;
                            TE_Spl_0  = TE_Spl;
                            TE_Plg_0  = TE_Plg;
    
                            % Radi
                            MRadi_0 = MRadi;
                            
                            ME_solid_TP_v0 = ME_solid_TP_v; % Volume fraction thermo phase
                            ME_solid_TP_w0 = ME_solid_TP_w; % Weight fraction thermo phase
                            ME_fluid_TP_v0 = ME_fluid_TP_v; % Volume fraction thermo phase
                            ME_fluid_TP_w0 = ME_fluid_TP_w; % Weight fraction thermo phase
    
                            MEGamma_fluid = interp1(XT(1:nxt+1:end,2),NGamma(1:nxt+1:end,1),ME_fluid0(:,2),'linear','extrap'); MEGamma_fluid = repmat(MEGamma_fluid,1,nTE);
                            MEGamma_solid = interp1(XT(1:nxt+1:end,2),NGamma(1:nxt+1:end,1),ME_solid0(:,2),'linear','extrap'); MEGamma_solid = repmat(MEGamma_solid,1,nTE);
    
                            MERho_fluid = interp1(XT(1:nxt+1:end,2),NRho0(1:nxt+1:end,2),ME_solid0(:,2),'linear','extrap');
                            MERho_solid = interp1(XT(1:nxt+1:end,2),NRho0(1:nxt+1:end,1),ME_solid0(:,2),'linear','extrap');
    
                        % Move particles
                            NV_mat_TE = NV_mat;
    
                            NV_mat_TE(NV_mat_TE(:,2)>NV_mat_TE(:,4),4) = NV_mat_TE(NV_mat_TE(:,2)>NV_mat_TE(:,4),2);
                            [ME_solid,ME_fluid] = reactivetransportREE(ME_solid0,ME_fluid0,NV_mat_TE,NDiv,gridxt,gridyt,nxt,0,nyt,0,TT,markmove,timestep,nTP);
                            [ME_fluid0_solid,ME_fluid0_fluid] = reactivetransportREE(ME_solid,ME_fluid,[-NV_mat_TE(:,3) -NV_mat_TE(:,4) -NV_mat_TE(:,3) -NV_mat_TE(:,4)],NDiv,gridxt,gridyt,nxt,0,nyt,0,TT,markmove,timestep,nTP);
                                
                            % Compute TP abundances (weight %)
                                ME_solid_TP_w = zeros(size(ME_solid,1),7);
    
                                ME_solid_TP_w(:,1)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,1),ME_solid0(:,2));
                                ME_solid_TP_w(:,2)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,2),ME_solid0(:,2));
                                ME_solid_TP_w(:,3)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,3),ME_solid0(:,2));
                                ME_solid_TP_w(:,4)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,4),ME_solid0(:,2));
                                ME_solid_TP_w(:,5)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,5),ME_solid0(:,2));
                                ME_solid_TP_w(:,6)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,6),ME_solid0(:,2));
                                ME_solid_TP_w(:,7)  = interp1(XT(1:nxt+1:end,2),NTp_wt0(1:nxt+1:end,7),ME_solid0(:,2));
                                ME_solid_TP_w       = ME_solid_TP_w./repmat(sum(ME_solid_TP_w,2),1,size(ME_solid_TP_w,2));
    
                                ME_solid_TP_w(:,1)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,1),ME_solid(:,2));
                                ME_solid_TP_w(:,2)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,2),ME_solid(:,2));
                                ME_solid_TP_w(:,3)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,3),ME_solid(:,2));
                                ME_solid_TP_w(:,4)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,4),ME_solid(:,2));
                                ME_solid_TP_w(:,5)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,5),ME_solid(:,2));
                                ME_solid_TP_w(:,6)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,6),ME_solid(:,2));
                                ME_solid_TP_w(:,7)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,7),ME_solid(:,2));
                                ME_solid_TP_w(ME_solid_TP_w(:,7)<tol_NTE,7)=0;
                                ME_solid_TP_w       = ME_solid_TP_w./repmat(sum(ME_solid_TP_w,2),1,size(ME_solid_TP_w,2));
    
                                ME_fluid_TP_w(:,1)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,1),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_w(:,2)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,2),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_w(:,3)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,3),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_w(:,4)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,4),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_w(:,5)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,5),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_w(:,6)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,6),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_w(:,7)  = interp1(XT(1:nxt+1:end,2),NTp_wt(1:nxt+1:end,7),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_w(ME_fluid_TP_w(:,7)<tol_NTE,7)=0;
                                ME_fluid_TP_w       = ME_fluid_TP_w./repmat(sum(ME_fluid_TP_w,2),1,size(ME_fluid_TP_w,2));
    
                            % Compute TP abundances (vol %)
                                ME_solid_TP_v = zeros(size(ME_solid,1),7);
        
                                ME_solid_TP_v(:,1)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,1),ME_solid0(:,2));
                                ME_solid_TP_v(:,2)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,2),ME_solid0(:,2));
                                ME_solid_TP_v(:,3)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,3),ME_solid0(:,2));
                                ME_solid_TP_v(:,4)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,4),ME_solid0(:,2));
                                ME_solid_TP_v(:,5)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,5),ME_solid0(:,2));
                                ME_solid_TP_v(:,6)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,6),ME_solid0(:,2));
                                ME_solid_TP_v(:,7)  = interp1(XT(1:nxt+1:end,2),NTp_vol0(1:nxt+1:end,7),ME_solid0(:,2));
                                ME_solid_TP_v        = ME_solid_TP_v./repmat(sum(ME_solid_TP_v,2),1,size(ME_solid_TP_v,2));
        
                                ME_solid_TP_v(:,1)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,1),ME_solid(:,2));
                                ME_solid_TP_v(:,2)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,2),ME_solid(:,2));
                                ME_solid_TP_v(:,3)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,3),ME_solid(:,2));
                                ME_solid_TP_v(:,4)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,4),ME_solid(:,2));
                                ME_solid_TP_v(:,5)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,5),ME_solid(:,2));
                                ME_solid_TP_v(:,6)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,6),ME_solid(:,2));
                                ME_solid_TP_v(:,7)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,7),ME_solid(:,2));
                                ME_solid_TP_v(ME_solid_TP_v(:,7)<tol_NTE,7)=0;
                                ME_solid_TP_v        = ME_solid_TP_v./repmat(sum(ME_solid_TP_v,2),1,size(ME_solid_TP_v,2));
        
                                ME_solid_TP_v0 = ME_solid_TP_v;
        
                                ME_fluid_TP_v(:,1)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,1),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_v(:,2)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,2),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_v(:,3)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,3),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_v(:,4)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,4),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_v(:,5)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,5),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_v(:,6)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,6),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_v(:,7)  = interp1(XT(1:nxt+1:end,2),NTp_vol(1:nxt+1:end,7),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_TP_v(ME_fluid_TP_v(:,7)<tol_NTE,7)=0;
                                ME_fluid_TP_v        = ME_fluid_TP_v./repmat(sum(ME_fluid_TP_v,2),1,size(ME_fluid_TP_v,2));
    
                        % Update phi
                                %Interpolate differences from nodes
                                FPhi_solid = scatteredInterpolant(XT(:,1),XT(:,2),NPhi(:,1)-NPhi0(:,1),'nearest','nearest');
                                ME_solid(:,3) = ME_solid0(:,3) + FPhi_solid(ME_solid0(:,1),ME_solid0(:,2));
                                ME_fluid(:,3) = ME_fluid0(:,3) - FPhi_solid(ME_fluid0(:,1),ME_fluid0(:,2));
    
                            ME_fluid(ME_fluid(:,3)<tol_NTE,3)=0;
                        
                            % ME_solid_rho & ME_fluid_rho
                                ME_solid_Rho_TP = zeros(size(ME_solid0,1),nTP+1);
                                ME_solid_Rho_TP(:,1)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,1),ME_solid0(:,2));
                                ME_solid_Rho_TP(:,2)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,2),ME_solid0(:,2));
                                ME_solid_Rho_TP(:,3)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,3),ME_solid0(:,2));
                                ME_solid_Rho_TP(:,4)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,4),ME_solid0(:,2));
                                ME_solid_Rho_TP(:,5)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,5),ME_solid0(:,2));
                                ME_solid_Rho_TP(:,6)  = interp1(XT(1:nxt+1:end,2),NRho_tp0(1:nxt+1:end,6),ME_solid0(:,2));

                                NRho0(NRho0(:,2)==0,2)  = NRho0(NRho0(:,2)==0,1);
                                ME_solid_Rho_TP(:,7)  = interp1(XT(1:nxt+1:end,2),NRho0(1:nxt+1:end,2),ME_solid0(:,2));
    
                                ME_solid_Rho_TP = zeros(size(ME_solid,1),nTP+1);
                                ME_solid_Rho_TP(:,1)  = interp1(XT(1:nxt+1:end,2),NRho_tp(1:nxt+1:end,1),ME_solid(:,2));
                                ME_solid_Rho_TP(:,2)  = interp1(XT(1:nxt+1:end,2),NRho_tp(1:nxt+1:end,2),ME_solid(:,2));
                                ME_solid_Rho_TP(:,3)  = interp1(XT(1:nxt+1:end,2),NRho_tp(1:nxt+1:end,3),ME_solid(:,2));
                                ME_solid_Rho_TP(:,4)  = interp1(XT(1:nxt+1:end,2),NRho_tp(1:nxt+1:end,4),ME_solid(:,2));
                                ME_solid_Rho_TP(:,5)  = interp1(XT(1:nxt+1:end,2),NRho_tp(1:nxt+1:end,5),ME_solid(:,2));
                                ME_solid_Rho_TP(:,6)  = interp1(XT(1:nxt+1:end,2),NRho_tp(1:nxt+1:end,6),ME_solid(:,2));

                                NRho(NRho(:,2)==0,2)  = NRho(NRho(:,2)==0,1);
                                ME_solid_Rho_TP(:,7)  = interp1(XT(1:nxt+1:end,2),NRho(1:nxt+1:end,2),ME_solid(:,2));
    
                                ME_fluid_Rho_TP = zeros(size(ME_fluid,1),nTP+1);
                                ME_fluid_Rho_TP(:,1)  = interp1(XT(1:nxt+1:end,2),NRho_tp(1:nxt+1:end,1),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_Rho_TP(:,2)  = interp1(XT(1:nxt+1:end,2),NRho_tp(1:nxt+1:end,2),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_Rho_TP(:,3)  = interp1(XT(1:nxt+1:end,2),NRho_tp(1:nxt+1:end,3),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_Rho_TP(:,4)  = interp1(XT(1:nxt+1:end,2),NRho_tp(1:nxt+1:end,4),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_Rho_TP(:,5)  = interp1(XT(1:nxt+1:end,2),NRho_tp(1:nxt+1:end,5),ME_fluid(:,2),'linear','extrap');
                                ME_fluid_Rho_TP(:,6)  = interp1(XT(1:nxt+1:end,2),NRho_tp(1:nxt+1:end,6),ME_fluid(:,2),'linear','extrap');
                                NRho(NRho(:,2)==0,2)  = NRho(NRho(:,2)==0,1);
                                ME_fluid_Rho_TP(:,7)  = interp1(XT(1:nxt+1:end,2),NRho(1:nxt+1:end,2),ME_fluid(:,2),'linear','extrap');
    
                                MERho_fluid = ME_solid_Rho_TP(:,7);
    
                        % Diffusion
    
                            % P-T in particles;
                                ME_T = interp1(XT(1:nxt+1:end,2),NT(1:nxt+1:end,1),ME_solid(:,2))+TK;
                                ME_P = interp1(XT(1:nxt+1:end,2),NP_mat(1:nxt+1:end,1),ME_solid(:,2))*PGPa;
    
                            % Coefficients
                                D_Oli  = repmat(D_coeff(1,:),size(ME_T,1),1);   E_Oli  = repmat(E_coeff(1,:),size(ME_T,1),1);   V_Oli = repmat(V_coeff(1,:),size(ME_T,1),1);
                                D_Cpx = repmat(D_coeff(2,:),size(ME_T,1),1);  E_Cpx = repmat(E_coeff(2,:),size(ME_T,1),1);  V_Cpx = repmat(V_coeff(2,:),size(ME_T,1),1);
                                D_Opx = repmat(D_coeff(3,:),size(ME_T,1),1);  E_Opx = repmat(E_coeff(3,:),size(ME_T,1),1);  V_Opx = repmat(V_coeff(3,:),size(ME_T,1),1);
                                D_Grt  = repmat(D_coeff(4,:),size(ME_T,1),1);   E_Grt = repmat(E_coeff(4,:),size(ME_T,1),1);    V_Grt = repmat(V_coeff(4,:),size(ME_T,1),1);
                                D_Spl  = repmat(D_coeff(5,:),size(ME_T,1),1);   E_Spl  = repmat(E_coeff(5,:),size(ME_T,1),1);   V_Spl = repmat(V_coeff(5,:),size(ME_T,1),1);
                                D_Plg  = repmat(D_coeff(6,:),size(ME_T,1),1);   E_Plg  = repmat(E_coeff(6,:),size(ME_T,1),1);   V_Plg = repmat(V_coeff(6,:),size(ME_T,1),1);
    
                            % Diffusivites
                                Mdiff_Oli    = double(D_Oli.*exp((-E_Oli + repmat(ME_P,1,nTE).*V_Oli)./(R.*repmat(ME_T,1,nTE))));
                                Mdiff_Cpx   = double(D_Cpx.*exp((-E_Cpx + repmat(ME_P,1,nTE).*V_Cpx)./(R.*repmat(ME_T,1,nTE))));
                                Mdiff_Opx   = double(D_Opx.*exp((-E_Opx + repmat(ME_P,1,nTE).*V_Opx)./(R.*repmat(ME_T,1,nTE))));
                                Mdiff_Grt    = double(D_Grt.*exp((-E_Grt + repmat(ME_P,1,nTE).*V_Grt)./(R.*repmat(ME_T,1,nTE))));
                                Mdiff_Spl    = double(D_Spl.*exp((-E_Spl + repmat(ME_P,1,nTE).*V_Spl)./(R.*repmat(ME_T,1,nTE))));
                                Mdiff_Plg    = double(D_Plg.*exp((-E_Plg + repmat(ME_P,1,nTE).*V_Plg)./(R.*repmat(ME_T,1,nTE))));
    
                            % Eu diffusivities in Cpx based on Sr diffusivities from Sneeringer et al. (1984)
                                if benchmark_Eu==1 %synthetic diopside
                                    D_Cpx_Eu2_ref = 1200*1e-4; % [m^2/s]
                                    E_Cpx_Eu2_ref = 122*4.184*1e3; % [J/mol]
                                elseif benchmark_Eu==0 %natural diopside
                                    D_Cpx_Eu2_ref = 54*1e-4; % [m^2/s]
                                    E_Cpx_Eu2_ref = 97*4.184*1e3; % [J/mol]
                                end
                                D_Cpx_Eu2 = repmat(D_Cpx_Eu2_ref,size(ME_T,1),1); 
                                E_Cpx_Eu2 = repmat(E_Cpx_Eu2_ref,size(ME_T,1),1);
                                Sn84a=repmat(2.94,size(ME_T,1),1); 
                                Sn84b=repmat(4640,size(ME_T,1),1);
                                Mdiff_Cpx_Eu2   = double(D_Cpx_Eu2.*exp(-E_Cpx_Eu2./(R.*ME_T)).*exp(-(ME_P.*1e-9)./((Sn84a.*ME_T-Sn84b).*R))); % [m^2/s]   
                                Mdiff_Cpx_Eu = Mdiff_Cpx_Eu2.*r_Eu + Mdiff_Cpx(:,Eu_pos).*(1-r_Eu); 
                                
                                Mdiff_Cpx = [Mdiff_Cpx(:,1:Eu_pos-1) Mdiff_Cpx_Eu2 Mdiff_Cpx(:,Eu_pos+1:nTE)];
                                if benchmark_2Cpx ==1; Mdiff_Oli=Mdiff_Cpx; end %WARNING allocating porphyroclastic matrix cpx to opx

                            % No diffusion if there is not fluid phase nearby
                                zeros_aux = zeros(size(Mdiff_Oli));
                                Mdiff_Oli(ME_solid(:,3)>1-tol_NTE,:)  = zeros_aux(ME_solid(:,3)>1-tol_NTE,:);
                                Mdiff_Cpx(ME_solid(:,3)>1-tol_NTE,:) = zeros_aux(ME_solid(:,3)>1-tol_NTE,:);
                                Mdiff_Opx(ME_solid(:,3)>1-tol_NTE,:) = zeros_aux(ME_solid(:,3)>1-tol_NTE,:);
                                Mdiff_Grt(ME_solid(:,3)>1-tol_NTE,:)  = zeros_aux(ME_solid(:,3)>1-tol_NTE,:);
                                Mdiff_Spl(ME_solid(:,3)>1-tol_NTE,:)  = zeros_aux(ME_solid(:,3)>1-tol_NTE,:);
                                Mdiff_Plg(ME_solid(:,3)>1-tol_NTE,:)  = zeros_aux(ME_solid(:,3)>1-tol_NTE,:);
    
                    % TRACE ELEMENTS
    
                        % Boundary condition
                            ME_fluid_solid = zeros(size(ME_solid,1),3+nTE);
                            ME_fluid_solid(:,1:2) = ME_solid(:,1:2);
                            ME_fluid_solid(:,3)   = 1-ME_solid(:,3);
                            ME_fluid_solid(ME_fluid_solid(:,3)<tol_NTE,7)=0;
    
                            for TE_pos = 4:3+nTE
                                [~,ia,~] = unique(ME_fluid(:,2));
                                ME_fluid_solid(:,TE_pos) = interp1(ME_fluid(ia,2),ME_fluid(ia,TE_pos),ME_solid(:,2),'nearest','extrap');
                            end
    
                            ME_RhoTP_solid0_adv = interp1(ME_solid0(:,2),ME_Rho_solid0.*ME_solid_TP_w0(:,7),ME_fluid0_solid(:,2),'linear','extrap');
                            ME_fluid_solid_adv = ME_fluid_solid;
                            ind_fluid = any(ME_fluid_solid(:,4:3+nTE),2);
    
                        % Prepare BC for next iteration
                            MBC_Oli  = MBC_Oli_0;
                            MBC_Cpx = MBC_Cpx_0;
                            MBC_Opx = MBC_Opx_0;
                            MBC_Grt  = MBC_Grt_0;
                            MBC_Spl  = MBC_Spl_0;
                            MBC_Plg  = MBC_Plg_0;
    
                            Kd_Oli = repmat(Kd_coeff(1,:),size(ME_solid,1),1);
                            Kd_Cpx = repmat(Kd_coeff(2,:),size(ME_solid,1),1);
                            Kd_Opx = repmat(Kd_coeff(3,:),size(ME_solid,1),1);
                            Kd_Grt = repmat(Kd_coeff(4,:),size(ME_solid,1),1);
                            Kd_Spl = repmat(Kd_coeff(5,:),size(ME_solid,1),1);
                            Kd_Plg = repmat(Kd_coeff(6,:),size(ME_solid,1),1);

                            if benchmark_2Cpx ==1; Kd_Oli=Kd_Cpx;end %WARNING allocating porphyroclastic matrix cpx to opx

                            MBC_Oli(ind_fluid,:)  = ME_fluid_solid(ind_fluid,4:3+nTE).*Kd_Oli(ind_fluid,:); 
                            MBC_Cpx(ind_fluid,:) = ME_fluid_solid(ind_fluid,4:3+nTE).*Kd_Cpx(ind_fluid,:);
                            MBC_Opx(ind_fluid,:) = ME_fluid_solid(ind_fluid,4:3+nTE).*Kd_Opx(ind_fluid,:);
                            MBC_Grt(ind_fluid,:)  = ME_fluid_solid(ind_fluid,4:3+nTE).*Kd_Grt(ind_fluid,:);
                            MBC_Spl(ind_fluid,:)  = ME_fluid_solid(ind_fluid,4:3+nTE).*Kd_Spl(ind_fluid,:);
                            MBC_Plg(ind_fluid,:)  = ME_fluid_solid(ind_fluid,4:3+nTE).*Kd_Plg(ind_fluid,:);
    
                        % Mass balance constraint - upper and lower bounds for melt's TE concentration
                            auxE_Oli   = repmat(ME_solid_TP_v(:,1)>0,1,nTE);
                            auxE_Cpx  = repmat(ME_solid_TP_v(:,2)>0,1,nTE);
                            auxE_Opx  = repmat(ME_solid_TP_v(:,3)>0,1,nTE);
                            auxE_Grt   = repmat(ME_solid_TP_v(:,4)>0,1,nTE);
                            auxE_Spl   = repmat(ME_solid_TP_v(:,5)>0,1,nTE);
                            auxE_Plg   = repmat(ME_solid_TP_v(:,6)>0,1,nTE);
                            auxE_Melt = repmat(ME_solid_TP_v(:,7)>0,1,nTE);
    
                        % Mass TP / vol T
                            ME_solid_RhoPhi_Oli  = repmat(ME_solid_TP_v(:,1).*ME_solid_Rho_TP(:,1),1,nTE);
                            ME_solid_RhoPhi_Cpx = repmat(ME_solid_TP_v(:,2).*ME_solid_Rho_TP(:,2),1,nTE);
                            ME_solid_RhoPhi_Opx = repmat(ME_solid_TP_v(:,3).*ME_solid_Rho_TP(:,3),1,nTE);
                            ME_solid_RhoPhi_Grt  = repmat(ME_solid_TP_v(:,4).*ME_solid_Rho_TP(:,4),1,nTE);
                            ME_solid_RhoPhi_Spl  = repmat(ME_solid_TP_v(:,5).*ME_solid_Rho_TP(:,5),1,nTE);
                            ME_solid_RhoPhi_Plg  = repmat(ME_solid_TP_v(:,6).*ME_solid_Rho_TP(:,6),1,nTE);
                            ME_solid_RhoPhi_Melt  = repmat(ME_solid_TP_v(:,7).*ME_solid_Rho_TP(:,7),1,nTE);
    
                        % Mass T / vol T
                            ME_Rho_solid = ME_solid_RhoPhi_Oli(:,1) + ME_solid_RhoPhi_Cpx(:,1) + ME_solid_RhoPhi_Opx(:,1) + ME_solid_RhoPhi_Grt(:,1) + ME_solid_RhoPhi_Spl(:,1) + ME_solid_RhoPhi_Plg(:,1) + ME_solid_RhoPhi_Melt(:,1);
    
                        % Mass TP / vol T
                            ME_fluid_RhoPhi_Oli  = repmat(ME_fluid_TP_v(:,1).*ME_fluid_Rho_TP(:,1),1,nTE);
                            ME_fluid_RhoPhi_Cpx = repmat(ME_fluid_TP_v(:,2).*ME_fluid_Rho_TP(:,2),1,nTE);
                            ME_fluid_RhoPhi_Opx = repmat(ME_fluid_TP_v(:,3).*ME_fluid_Rho_TP(:,3),1,nTE);
                            ME_fluid_RhoPhi_Grt  = repmat(ME_fluid_TP_v(:,4).*ME_fluid_Rho_TP(:,4),1,nTE);
                            ME_fluid_RhoPhi_Spl  = repmat(ME_fluid_TP_v(:,5).*ME_fluid_Rho_TP(:,5),1,nTE);
                            ME_fluid_RhoPhi_Plg  = repmat(ME_fluid_TP_v(:,6).*ME_fluid_Rho_TP(:,6),1,nTE);
                            ME_fluid_RhoPhi_Melt  = repmat(ME_fluid_TP_v(:,7).*ME_fluid_Rho_TP(:,7),1,nTE);
    
                        % Mass T / vol T
                            ME_Rho_fluid = ME_fluid_RhoPhi_Oli(:,1) + ME_fluid_RhoPhi_Cpx(:,1) + ME_fluid_RhoPhi_Opx(:,1) + ME_fluid_RhoPhi_Grt(:,1) + ME_fluid_RhoPhi_Spl(:,1) + ME_fluid_RhoPhi_Plg(:,1) + ME_fluid_RhoPhi_Melt(:,1);
    
                        % Mass TP / mass T
                            ME_solid_w_Oli   = repmat(ME_solid_TP_w(:,1),1,nTE); ME_solid_RhoPhi_Oli = ME_solid_w_Oli;
                            ME_solid_w_Cpx  = repmat(ME_solid_TP_w(:,2),1,nTE); ME_solid_RhoPhi_Cpx = ME_solid_w_Cpx;
                            ME_solid_w_Opx  = repmat(ME_solid_TP_w(:,3),1,nTE); ME_solid_RhoPhi_Opx = ME_solid_w_Opx;
                            ME_solid_w_Grt   = repmat(ME_solid_TP_w(:,4),1,nTE); ME_solid_RhoPhi_Grt = ME_solid_w_Grt;
                            ME_solid_w_Spl   = repmat(ME_solid_TP_w(:,5),1,nTE); ME_solid_RhoPhi_Spl = ME_solid_w_Spl;
                            ME_solid_w_Plg   = repmat(ME_solid_TP_w(:,6),1,nTE); ME_solid_RhoPhi_Plg = ME_solid_w_Plg;
                            ME_solid_w_Melt = repmat(ME_solid_TP_w(:,7),1,nTE); ME_solid_RhoPhi_Melt = ME_solid_w_Melt;
    
                        % Solid [mass solid / mass T]
                            ME_solid_w_solid = ME_solid_w_Oli + ME_solid_w_Cpx + ME_solid_w_Opx + ME_solid_w_Grt + ME_solid_w_Spl + ME_solid_w_Plg;  % normalize to 1 every TP contribution
    
                        % Factor [virtual total mass / total mass]
                            ME_solid_w_total = ME_solid_w_solid+ME_solid_w_Melt;
    
                            Melt_TE = ME_TE./(auxE_Oli.*Kd_Oli.*ME_solid_w_Oli./ME_solid_w_total    ...
                                                + auxE_Cpx.*Kd_Cpx.*ME_solid_w_Cpx./ME_solid_w_total ...
                                                + auxE_Opx.*Kd_Opx.*ME_solid_w_Opx./ME_solid_w_total ...
                                                + auxE_Grt.*Kd_Grt.*ME_solid_w_Grt./ME_solid_w_total    ...
                                                + auxE_Spl.*Kd_Spl.*ME_solid_w_Spl./ME_solid_w_total    ...
                                                + auxE_Plg.*Kd_Plg.*ME_solid_w_Plg./ME_solid_w_total    ...
                                                + auxE_Melt.*ME_solid_w_Melt./ME_solid_w_total);    % [melt] = [solid]/Kd_bulk    (normalized to solid)
    
                            Melt_TE = ME_TE./(auxE_Oli.*Kd_Oli.*ME_solid_w_Oli./ME_solid_w_solid    ...
                                                + auxE_Cpx.*Kd_Cpx.*ME_solid_w_Cpx./ME_solid_w_solid ...
                                                + auxE_Opx.*Kd_Opx.*ME_solid_w_Opx./ME_solid_w_solid ...
                                                + auxE_Grt.*Kd_Grt.*ME_solid_w_Grt./ME_solid_w_solid    ...
                                                + auxE_Spl.*Kd_Spl.*ME_solid_w_Spl./ME_solid_w_solid    ...
                                                + auxE_Plg.*Kd_Plg.*ME_solid_w_Plg./ME_solid_w_solid);    % [melt] = [solid]/Kd_bulk    (normalized to solid)
    
                            ME_melt_eq  = zeros(size(ME_solid,1),nTE);                 
                            ME_melt_eq(ME_solid_TP_v(:,7)>0,:) = Melt_TE(ME_solid_TP_v(:,7)>0,:);  % TE mass / fluid mass
    
                            % TE mass / TP mass
                            ME_Oli_eq   =  Melt_TE.*auxE_Oli.*Kd_Oli;           
                            ME_Cpx_eq  =  Melt_TE.*auxE_Cpx.*Kd_Cpx ;   
                            ME_Opx_eq  =  Melt_TE.*auxE_Opx.*Kd_Opx;       
                            ME_Grt_eq   =  Melt_TE.*auxE_Grt.*Kd_Grt;         
                            ME_Spl_eq   =  Melt_TE.*auxE_Spl.*Kd_Spl;        
                            ME_Plg_eq   =  Melt_TE.*auxE_Plg.*Kd_Plg;  
    
                            correct_TE_met = logical(repmat(ME_solid_TP_v(:,7)<tol_NTE,1,nTE));
                            correct_TE_met_aux = logical(0*correct_TE_met);
    
                        % TE compositions per mineral 
                                for TE_pos = 1:nTE
                                    TE = TE_list(TE_pos);
                                    TEprofile_Oli.(TE) = repmat(ME_Oli_eq(:,TE_pos)',nc,1)';
                                    TEprofile_Cpx.(TE) = repmat(ME_Cpx_eq(:,TE_pos)',nc,1)';
                                    TEprofile_Opx.(TE) = repmat(ME_Opx_eq(:,TE_pos)',nc,1)';
                                    TEprofile_Grt.(TE) = repmat(ME_Grt_eq(:,TE_pos)',nc,1)';
                                    TEprofile_Spl.(TE) = repmat(ME_Spl_eq(:,TE_pos)',nc,1)';
                                    TEprofile_Plg.(TE) = repmat(ME_Plg_eq(:,TE_pos)',nc,1)';
                                    if TE_pos == 1
                                        Mat_Oli_aux = TEprofile_Oli.(TE);
                                        Mat_Cpx_aux = TEprofile_Cpx.(TE);
                                        Mat_Opx_aux = TEprofile_Opx.(TE);
                                        Mat_Grt_aux = TEprofile_Grt.(TE);
                                        Mat_Spl_aux = TEprofile_Spl.(TE);
                                        Mat_Plg_aux = TEprofile_Plg.(TE);
                                    else
                                        Mat_Oli_aux = cat(2,Mat_Oli_aux,TEprofile_Oli.(TE));
                                        Mat_Cpx_aux = cat(2,Mat_Cpx_aux,TEprofile_Cpx.(TE));
                                        Mat_Opx_aux = cat(2,Mat_Opx_aux,TEprofile_Opx.(TE));
                                        Mat_Grt_aux = cat(2,Mat_Grt_aux,TEprofile_Grt.(TE));
                                        Mat_Spl_aux = cat(2,Mat_Spl_aux,TEprofile_Spl.(TE));
                                        Mat_Plg_aux = cat(2,Mat_Plg_aux,TEprofile_Plg.(TE));
                                    end
                                end
    
                        % Prepare for iterations
                            Error_Oli = [];
                            Error_Cpx = [];
                            Error_Opx = [];
                            Error_Grt = [];
                            Error_Spl = [];
                            Error_Plg = [];
    
                            MBC_Oli_iter = [];
                            MBC_Cpx_iter = [];
                            MBC_Opx_iter = [];
                            MBC_Grt_iter = [];
                            MBC_Spl_iter = [];
                            MBC_Plg_iter = [];
    
                            MGamma_e_solid_iter = zeros(size(ME_solid,1),nTE);
                            MGamma_e_solid_iter_aux = [];
                            MGamma_e_solid_max =  zeros(size(ME_solid,1),nTE);
                            MGamma_e_solid_min =  zeros(size(ME_solid,1),nTE);
                            MGamma_e_solid_max_aux =  zeros(size(ME_solid,1),nTE);
    
                            ME_fluid_aux_iter_aux = [];
                            ind_gain_loss = logical(zeros(size(ME_solid,1),nTE));
                            ind_loss_gain = logical(zeros(size(ME_solid,1),nTE));
    
                            TE_Solid_iter = [];
                            TE_Fluid_iter = [];
                            index_double = [];
                            ind_gain = logical(zeros(size(ME_solid,1),nTE));
                            ind_loss = logical(zeros(size(ME_solid,1),nTE));
                            ind_negative_total = logical(zeros(size(ME_solid,1),nTE));
    
                            damp = 0.5;
                            correct_TE_ind = ones(size(MBC_Oli))>0;
                            correct_TE_ind = correct_TE_ind & correct_TE_met==0;

                            tol_noChange = tol_TE*tol_NTE/timestep;

                            ME_Oli_noChange = cell2mat(ME_Oli);
                            ME_Cpx_noChange = cell2mat(ME_Cpx);
                            ME_Opx_noChange = cell2mat(ME_Opx);
                            ME_Grt_noChange = cell2mat(ME_Grt);
                            ME_Spl_noChange = cell2mat(ME_Spl);
                            ME_Plg_noChange = cell2mat(ME_Plg);

                            TE_Oli_noChange = TE_Oli;
                            TE_Cpx_noChange = TE_Cpx; 
                            TE_Opx_noChange = TE_Opx; 
                            TE_Grt_noChange = TE_Grt; 
                            TE_Spl_noChange = TE_Spl; 
                            TE_Plg_noChange = TE_Plg;

                            ME_solid_noChange = ME_solid(:,4:3+nTE);
                            aux_negative_ind = [];
    
                        % TE iterations
                            iter_TE = 0; exit_TE = 0;
                            while exit_TE == 0
                                iter_TE = iter_TE +1;
    
                                % Compute diffusion
                                    % Input:
                                        %   -> Boundary conditions t_n and t_{n+1}. 
                                        %               E.g. NOli_eq(nodes) and BC_Oli_0(particles)
                                        %   -> X,Y coordinates of particles (MR_solid) 
                                        %   -> Grain densities (MR_solid)
                                        %   -> Proportion (vol%) of TP inside solid (MR_solid_TP_w)
                                        %   -> Profiles. E.g. MR_Oli
                                        %   -> Diffusion and timestep
                                    % Output:
                                        %   -> Gradients of every oxide inside TP
                                        %   -> New profiles
                                        %   -> Boundary conditions for next time step.
    
                                        tic
                                            correct_TE_met(correct_TE_met_aux) = true;
                                            ind_correct_TE_met_aux = correct_TE_met_aux;
    
                                            [Gamma_Oli,Gamma_Cpx,Gamma_Opx,Gamma_Grt,Gamma_Spl,Gamma_Plg,ME_Oli,ME_Cpx,ME_Opx,ME_Grt,ME_Spl,ME_Plg,TE_Oli,TE_Cpx,TE_Opx,TE_Grt,TE_Spl,TE_Plg,MRadi,MRadi_0,ME_solid] = ...
                                                diffusion1D(ME_Oli_eq,ME_Cpx_eq,ME_Opx_eq,ME_Grt_eq,ME_Spl_eq,ME_Plg_eq,ME_melt_eq,MBC_Oli_0,MBC_Cpx_0,MBC_Opx_0,MBC_Grt_0,MBC_Spl_0,MBC_Plg_0,MBC_Oli,MBC_Cpx,MBC_Opx,MBC_Grt,MBC_Spl,MBC_Plg,ME_solid,ME_solid0,ME_solid_TP_v,ME_solid_TP_v0,ME_solid_TP_w,ME_Oli_0,ME_Cpx_0,ME_Opx_0,ME_Grt_0,ME_Spl_0,ME_Plg_0,ME_solid_Rho_TP,ME_solid_Rho_TP0,Mdiff_Oli,Mdiff_Cpx,Mdiff_Opx,Mdiff_Grt,Mdiff_Spl,Mdiff_Plg,timestep,nc,nTE,diff_nTP,diff_nTE,fix_radi,correct_TE_ind,correct_TE_met);
                                        toc
    
                                % Leave unchanged where nothing happens
                                    if iter_TE == 1
                                        ind_noChange = abs(Gamma_Oli+Gamma_Cpx+Gamma_Opx+Gamma_Grt+Gamma_Spl+Gamma_Plg) < tol_noChange; % ind_noChange(:,end-1:end) = [];
                                        ind_noChange = logical(0*ind_noChange);
                                    end
    
                                    mat_ind_noChange = reshape(ind_noChange',[],1);
    
                                    mat_ME_Oli = cell2mat(ME_Oli);
                                    mat_ME_Cpx = cell2mat(ME_Cpx);
                                    mat_ME_Opx = cell2mat(ME_Opx);
                                    mat_ME_Grt = cell2mat(ME_Grt);
                                    mat_ME_Spl = cell2mat(ME_Spl);
                                    mat_ME_Plg = cell2mat(ME_Plg);
    
                                    row_cell = nTE*ones(size(mat_ME_Oli,1)/nTE,1);
    
                                    mat_ME_Oli(mat_ind_noChange,:)    = ME_Oli_noChange(mat_ind_noChange,:);    mat_ME_Cpx(mat_ind_noChange,:)   = ME_Cpx_noChange(mat_ind_noChange,:);    mat_ME_Opx(mat_ind_noChange,:)    = ME_Opx_noChange(mat_ind_noChange,:);    mat_ME_Grt(mat_ind_noChange,:)    = ME_Grt_noChange(mat_ind_noChange,:);    mat_ME_Spl(mat_ind_noChange,:)    = ME_Spl_noChange(mat_ind_noChange,:);    mat_ME_Plg(mat_ind_noChange,:)    = ME_Plg_noChange(mat_ind_noChange,:);
                                    ME_Oli    = mat2cell(mat_ME_Oli,row_cell,nc);    ME_Cpx    = mat2cell(mat_ME_Cpx,row_cell,nc);    ME_Opx    = mat2cell(mat_ME_Opx,row_cell,nc);    ME_Grt    = mat2cell(mat_ME_Grt,row_cell,nc);    ME_Spl    = mat2cell(mat_ME_Spl,row_cell,nc);    ME_Plg    = mat2cell(mat_ME_Plg,row_cell,nc);
                                    TE_Oli(ind_noChange) = TE_Oli_noChange(ind_noChange);   TE_Cpx(ind_noChange) = TE_Cpx_noChange(ind_noChange);   TE_Opx(ind_noChange) = TE_Opx_noChange(ind_noChange);   TE_Grt(ind_noChange) = TE_Grt_noChange(ind_noChange);   TE_Spl(ind_noChange) = TE_Spl_noChange(ind_noChange);   TE_Plg(ind_noChange) = TE_Plg_noChange(ind_noChange);
                                    ME_solid_noChange_aux = ME_solid(:,4:3+nTE);
                                    ME_solid_noChange_aux(ind_noChange) = ME_solid_noChange(ind_noChange);
                                    ME_solid(:,4:3+nTE) = ME_solid_noChange_aux;
    
                                % Update where something happens
                                    if iter_TE==1
                                        mat_ME_Oli = cell2mat(ME_Oli);
                                        mat_ME_Cpx = cell2mat(ME_Cpx);
                                        mat_ME_Opx = cell2mat(ME_Opx);
                                        mat_ME_Grt = cell2mat(ME_Grt);
                                        mat_ME_Spl = cell2mat(ME_Spl);
                                        mat_ME_Plg = cell2mat(ME_Plg);
    
                                        Gamma_Oli_iter = Gamma_Oli; Gamma_Cpx_iter = Gamma_Cpx; Gamma_Opx_iter = Gamma_Opx; Gamma_Grt_iter = Gamma_Grt; Gamma_Spl_iter = Gamma_Spl; Gamma_Plg_iter = Gamma_Plg;
                                        ME_Oli_iter    = mat_ME_Oli;    ME_Cpx_iter    = mat_ME_Cpx;    ME_Opx_iter    = mat_ME_Opx;    ME_Grt_iter    = mat_ME_Grt;    ME_Spl_iter    = mat_ME_Spl;    ME_Plg_iter    = mat_ME_Plg;    
                                        TE_Oli_iter   = TE_Oli;   TE_Cpx_iter   = TE_Cpx;   TE_Opx_iter   = TE_Opx;   TE_Grt_iter   = TE_Grt;   TE_Spl_iter   = TE_Spl;   TE_Plg_iter   = TE_Plg;   
                                        ME_solid_iter = ME_solid(:,4:3+nTE);
                                    else
                                        mat_ME_Oli = cell2mat(ME_Oli);
                                        mat_ME_Cpx = cell2mat(ME_Cpx);
                                        mat_ME_Opx = cell2mat(ME_Opx);
                                        mat_ME_Grt = cell2mat(ME_Grt);
                                        mat_ME_Spl = cell2mat(ME_Spl);
                                        mat_ME_Plg = cell2mat(ME_Plg);
                                        ME_solid_aux = ME_solid(:,4:3+nTE);
    
                                        % Correct this iteration
                                            mat_correct_TE_ind = reshape(correct_TE_ind',[],1);
                                            Gamma_Oli_correct = Gamma_Oli_iter; Gamma_Cpx_correct = Gamma_Cpx_iter; Gamma_Opx_correct = Gamma_Opx_iter; Gamma_Grt_correct = Gamma_Grt_iter; Gamma_Spl_correct = Gamma_Spl_iter; Gamma_Plg_correct = Gamma_Plg_iter;
                                            ME_Oli_correct    = ME_Oli_iter;    ME_Cpx_correct   = ME_Cpx_iter;    ME_Opx_correct    = ME_Opx_iter;    ME_Grt_correct    = ME_Grt_iter;    ME_Spl_correct    = ME_Spl_iter;    ME_Plg_correct    = ME_Plg_iter;    
                                            TE_Oli_correct   = TE_Oli_iter;   TE_Cpx_correct  = TE_Cpx_iter;   TE_Opx_correct   = TE_Opx_iter;   TE_Grt_correct   = TE_Grt_iter;   TE_Spl_correct   = TE_Spl_iter;   TE_Plg_correct   = TE_Plg_iter;   
                                            ME_solid_correct = ME_solid_iter;
                                            Gamma_Oli_correct(correct_TE_ind) = Gamma_Oli(correct_TE_ind); Gamma_Cpx_correct(correct_TE_ind) = Gamma_Cpx(correct_TE_ind); Gamma_Opx_correct(correct_TE_ind) = Gamma_Opx(correct_TE_ind); Gamma_Grt_correct(correct_TE_ind) = Gamma_Grt(correct_TE_ind); Gamma_Spl_correct(correct_TE_ind) = Gamma_Spl(correct_TE_ind); Gamma_Plg_correct(correct_TE_ind) = Gamma_Plg(correct_TE_ind);
                                            ME_Oli_correct(mat_correct_TE_ind,:)    = mat_ME_Oli(mat_correct_TE_ind,:);    ME_Cpx_correct(mat_correct_TE_ind,:)   = mat_ME_Cpx(mat_correct_TE_ind,:);    ME_Opx_correct(mat_correct_TE_ind,:)    = mat_ME_Opx(mat_correct_TE_ind,:);    ME_Grt_correct(mat_correct_TE_ind,:)    = mat_ME_Grt(mat_correct_TE_ind,:);    ME_Spl_correct(mat_correct_TE_ind,:)    = mat_ME_Spl(mat_correct_TE_ind,:);    ME_Plg_correct(mat_correct_TE_ind,:)    = mat_ME_Plg(mat_correct_TE_ind,:);    
                                            TE_Oli_correct(correct_TE_ind)   = TE_Oli(correct_TE_ind);   TE_Cpx_correct(correct_TE_ind)   = TE_Cpx(correct_TE_ind);   TE_Opx_correct(correct_TE_ind)   = TE_Opx(correct_TE_ind);   TE_Grt_correct(correct_TE_ind)   = TE_Grt(correct_TE_ind);   TE_Spl_correct(correct_TE_ind)   = TE_Spl(correct_TE_ind);   TE_Plg_correct(correct_TE_ind)   = TE_Plg(correct_TE_ind);   
                                            ME_solid_correct(correct_TE_ind) = ME_solid_aux(correct_TE_ind);
    
                                        % Update for next iteration
                                            Gamma_Oli_iter = Gamma_Oli_correct; Gamma_Cpx_iter = Gamma_Cpx_correct; Gamma_Opx_iter = Gamma_Opx_correct; Gamma_Grt_iter = Gamma_Grt_correct; Gamma_Spl_iter = Gamma_Spl_correct; Gamma_Plg_iter = Gamma_Plg_correct;
                                            ME_Oli_iter    = ME_Oli_correct;    ME_Cpx_iter    = ME_Cpx_correct;    ME_Opx_iter    = ME_Opx_correct;    ME_Grt_iter    = ME_Grt_correct;    ME_Spl_iter    = ME_Spl_correct;    ME_Plg_iter    = ME_Plg_correct;    
                                            TE_Oli_iter   = TE_Oli_correct;   TE_Cpx_iter   = TE_Cpx_correct;   TE_Opx_iter   = TE_Opx_correct;   TE_Grt_iter   = TE_Grt_correct;   TE_Spl_iter   = TE_Spl_correct;   TE_Plg_iter   = TE_Plg_correct;   
                                            ME_solid_iter = ME_solid_correct;
    
                                        % Save data
                                            Gamma_Oli = Gamma_Oli_correct; Gamma_Cpx = Gamma_Cpx_correct; Gamma_Opx = Gamma_Opx_correct; Gamma_Grt = Gamma_Grt_correct; Gamma_Spl = Gamma_Spl_correct; Gamma_Plg = Gamma_Plg_correct;
                                            ME_Oli    = mat2cell(ME_Oli_correct,row_cell,nc);    ME_Cpx    = mat2cell(ME_Cpx_correct,row_cell,nc);    ME_Opx    = mat2cell(ME_Opx_correct,row_cell,nc);    ME_Grt    = mat2cell(ME_Grt_correct,row_cell,nc);    ME_Spl    = mat2cell(ME_Spl_correct,row_cell,nc);    ME_Plg    = mat2cell(ME_Plg_correct,row_cell,nc);    
                                            TE_Oli   = TE_Oli_correct;   TE_Cpx   = TE_Cpx_correct;   TE_Opx   = TE_Opx_correct;   TE_Grt   = TE_Grt_correct;   TE_Spl   = TE_Spl_correct;   TE_Plg   = TE_Plg_correct;   
                                            ME_solid(:,4:3+nTE) = ME_solid_correct;
                                    end
    
                                    MGamma_e_solid = Gamma_Oli+Gamma_Cpx+Gamma_Opx+Gamma_Grt+Gamma_Spl+Gamma_Plg;  MGamma_e_solid(ind_noChange)=0;
                                    MGamma_e_solid_save = MGamma_e_solid;
    
                                    % Solid gains TE
                                    ind_gain = MGamma_e_solid-MGamma_e_solid_iter>0;
                                    % Solid losses TE
                                    ind_loss = MGamma_e_solid-MGamma_e_solid_iter<0;
    
                                    if iter_TE>1
                                        ind_gain_loss(ind_gain_loss==0 & (ind_gain & ind_loss_iter)) = 1;  % before loss now switch to gain
                                        ind_loss_gain(ind_loss_gain==0 & (ind_loss & ind_gain_iter)) = 1;  % before gain now switch to loss
    
                                        MGamma_e_solid_min(ind_loss & ind_loss_iter & ind_gain_loss==0) = MGamma_e_solid(ind_loss & ind_loss_iter & ind_gain_loss==0);
                                        MGamma_e_solid_max(ind_loss & ind_loss_iter) = MGamma_e_solid_iter(ind_loss & ind_loss_iter);
    
                                        MGamma_e_solid_min(ind_gain & ind_loss_iter) = MGamma_e_solid_iter(ind_gain & ind_loss_iter); 
                                        MGamma_e_solid_max(ind_loss & ind_gain_iter) = MGamma_e_solid_iter(ind_loss & ind_gain_iter);
    
                                        MGamma_e_solid_min(ind_gain & ind_gain_iter) = MGamma_e_solid_iter(ind_gain & ind_gain_iter);
                                        MGamma_e_solid_max(ind_gain & ind_gain_iter & ind_loss_gain==0) = MGamma_e_solid(ind_gain & ind_gain_iter & ind_loss_gain==0);
    
                                        MGamma_e_solid_max(MGamma_e_solid_max>MGamma_e_solid_max_abs) = MGamma_e_solid_max_abs(MGamma_e_solid_max>MGamma_e_solid_max_abs);
    
                                        MGamma_e_solid(ME_solid_TP_v(:,7)>0,:) = MGamma_e_solid_min(ME_solid_TP_v(:,7)>0,:)+damp*(MGamma_e_solid_max(ME_solid_TP_v(:,7)>0,:)-MGamma_e_solid_min(ME_solid_TP_v(:,7)>0,:));
                                    end
    
                                    % Solid TE gain-loss
                                        ind_gain_iter = ind_gain;
                                        ind_loss_iter = ind_loss;
    
                                    MGamma_e_fluid = MGamma_e_solid;
                                    damping=1;
                                    % Update ME_fluid value
                                            a_fluid = 0*MEGamma_solid./(repmat(ME_fluid_solid(:,3).*MERho_fluid,1,nTE));          a_fluid(isnan(a_fluid))=0;
                                            a_fluid = 0*MEGamma_solid./(repmat(ME_Rho_solid.*ME_solid_TP_w(:,7),1,nTE));          a_fluid(isnan(a_fluid))=0;
                                            b_fluid = -damping*MGamma_e_fluid./(repmat(ME_fluid_solid(:,3).*MERho_fluid,1,nTE));  b_fluid(isnan(b_fluid))=0;
                                            b_fluid = -damping*MGamma_e_fluid./(repmat(ME_Rho_solid.*ME_solid_TP_w(:,7),1,nTE));  b_fluid(isnan(b_fluid))=0;
                                            c_fluid_aux = ME_fluid_solid0(:,4:3+nTE);
                                            c_fluid = (repmat(ME_RhoTP_solid0_adv,1,nTE)).*ME_fluid_solid_adv(:,4:3+nTE)./(repmat(ME_Rho_solid.*ME_solid_TP_w(:,7),1,nTE));  c_fluid(isnan(c_fluid))=0;  c_fluid(isinf(c_fluid))=0;
                                            c_fluid0 = (repmat(ME_Rho_solid0.*ME_solid_TP_w0(:,7),1,nTE)).*ME_fluid_solid0(:,4:3+nTE)./(repmat(ME_Rho_solid.*ME_solid_TP_w(:,7),1,nTE));  c_fluid0(isnan(c_fluid0))=0;  c_fluid0(isinf(c_fluid0))=0;
                                            dt_fluid = timestep;
    
                                            ME_fluid_aux = (b_fluid - exp(-a_fluid*dt_fluid).*(b_fluid - a_fluid.*c_fluid))./a_fluid;
                                            ME_fluid_aux(a_fluid==0) = c_fluid(a_fluid==0) + dt_fluid*b_fluid(a_fluid==0);
                                            ME_fluid_aux = c_fluid + dt_fluid*b_fluid;
                                            MGamma_e_solid_max_abs = c_fluid.*(repmat(ME_Rho_solid.*ME_solid_TP_w(:,7),1,nTE))/dt_fluid;
    
                                            ME_fluid_aux(isnan(ME_fluid_aux))=c_fluid(isnan(ME_fluid_aux)); 
                                            ME_fluid_aux(isinf(ME_fluid_aux))=c_fluid(isinf(ME_fluid_aux));
    
                                    if iter_TE==1
                                        correct_TE_met_aux = ME_fluid_aux<0 & c_fluid==0;
                                        correct_TE_met_aux = logical(repmat(any(correct_TE_met_aux,2),1,nTE));
                                    end
    
                                    ME_fluid_aux(ME_fluid_aux<0 & c_fluid==0)=c_fluid(ME_fluid_aux<0 & c_fluid==0);
    
                                    % check if there are negative values and correct them
                                        iter_negative = 0;
                                        iter_negative_max = 30;
                                        ind_negative = ME_fluid_aux<0 & correct_TE_ind;
                                        ind_negative_total = ind_negative | ind_negative_total;
                                        aux_negative_ind(iter_TE) = sum(sum(ind_negative));
    
                                        while sum(sum(ind_negative))
                                            iter_negative = iter_negative + 1;
    
                                            MGamma_e_solid_max_aux(ind_negative) = MGamma_e_solid(ind_negative); MGamma_e_solid_max_aux(ind_negative & (MGamma_e_solid_max_abs<MGamma_e_fluid)) = MGamma_e_solid_max_abs(ind_negative & (MGamma_e_solid_max_abs<MGamma_e_fluid));
                                            MGamma_e_solid(ind_negative) = MGamma_e_solid_min(ind_negative)+damp/2*(MGamma_e_solid_max_aux(ind_negative)-MGamma_e_solid_min(ind_negative));
                                            MGamma_e_fluid = MGamma_e_solid;
    
                                            b_fluid = -damping*MGamma_e_fluid./(repmat(ME_fluid_solid0(:,3).*MERho_fluid,1,nTE));  b_fluid(isnan(b_fluid))=0;  b_fluid(isinf(b_fluid))=0;
                                            b_fluid = -damping*MGamma_e_fluid./(repmat(ME_Rho_solid.*ME_solid_TP_w(:,7),1,nTE));  b_fluid(isnan(b_fluid))=0;  b_fluid(isinf(b_fluid))=0;
    
                                            ME_fluid_aux(a_fluid==0) = c_fluid(a_fluid==0) + dt_fluid*b_fluid(a_fluid==0);
    
                                            if iter_negative == iter_negative_max
                                                ME_fluid_aux(ind_negative) = ME_fluid_aux_iter(ind_negative);
                                            end
    
                                            ind_negative = ME_fluid_aux<0 & correct_TE_ind;
    
                                        end
    
                                        if iter_negative>0
                                            fprintf('Negative fluid composition reached at iteration %d. And exit in %d inner iterations\n',iter_TE,iter_negative);
                                        end
    
                                    if iter_TE == 1
                                        MGamma_e_solid_min(ind_loss) = MGamma_e_solid(ind_loss);
                                        MGamma_e_solid_max(ind_gain) = MGamma_e_solid(ind_gain);
                                    end
    
                                    ME_fluid_aux(isnan(ME_fluid_aux))=c_fluid(isnan(ME_fluid_aux));
                                    ME_fluid_aux(isinf(ME_fluid_aux))=c_fluid(isinf(ME_fluid_aux));
                                    ME_fluid_aux(correct_TE_met)=0;
    
                                    if iter_TE==1
                                        ME_fluid_aux_iter = ME_fluid_aux;
                                    else
                                        ME_fluid_correct = ME_fluid_aux_iter;
                                        ME_fluid_correct(correct_TE_ind) = ME_fluid_aux(correct_TE_ind); 
                                        ME_fluid_correct(correct_TE_met)=0;
                                        ME_fluid_aux_iter = ME_fluid_correct;
                                    end
    
                                    MGamma_e_solid_iter_aux = [MGamma_e_solid_iter_aux MGamma_e_solid]  ;  
                                    MGamma_e_solid_iter = MGamma_e_solid;
                                    ME_fluid_aux_iter_aux = [ME_fluid_aux_iter_aux ME_fluid_aux_iter];  
    
                                ME_fluid_solid(:,4:3+nTE) = ME_fluid_aux_iter; 
    
                                % Prepare BC for next iteration
                                    MBC_Oli(ME_solid(:,3)<1-tol_NTE,:)  = ME_fluid_solid((ME_solid(:,3)<1-tol_NTE),4:3+nTE).*Kd_Oli((ME_solid(:,3)<1-tol_NTE),:); 
                                    MBC_Cpx(ME_solid(:,3)<1-tol_NTE,:) = ME_fluid_solid((ME_solid(:,3)<1-tol_NTE),4:3+nTE).*Kd_Cpx((ME_solid(:,3)<1-tol_NTE),:);
                                    MBC_Opx(ME_solid(:,3)<1-tol_NTE,:) = ME_fluid_solid((ME_solid(:,3)<1-tol_NTE),4:3+nTE).*Kd_Opx((ME_solid(:,3)<1-tol_NTE),:);
                                    MBC_Grt(ME_solid(:,3)<1-tol_NTE,:)  = ME_fluid_solid((ME_solid(:,3)<1-tol_NTE),4:3+nTE).*Kd_Grt((ME_solid(:,3)<1-tol_NTE),:);
                                    MBC_Spl(ME_solid(:,3)<1-tol_NTE,:)  = ME_fluid_solid((ME_solid(:,3)<1-tol_NTE),4:3+nTE).*Kd_Spl((ME_solid(:,3)<1-tol_NTE),:);
                                    MBC_Plg(ME_solid(:,3)<1-tol_NTE,:)  = ME_fluid_solid((ME_solid(:,3)<1-tol_NTE),4:3+nTE).*Kd_Plg((ME_solid(:,3)<1-tol_NTE),:);
    
                                    ME_Oli_mat  = cell2mat(ME_Oli); 
                                    ME_Cpx_mat = cell2mat(ME_Cpx);
                                    ME_Opx_mat = cell2mat(ME_Opx);
                                    ME_Grt_mat  = cell2mat(ME_Grt);
                                    ME_Spl_mat  = cell2mat(ME_Spl); 
                                    ME_Plg_mat  = cell2mat(ME_Plg);
    
                                    % TE compositions per mineral 
                                        for TE_pos = 1:nTE
    
                                            TE = TE_list(TE_pos,:);
    
                                            ME_Oli_mat_TE.(TE)  = ME_Oli_mat(TE_pos:nTE:end,:); 
                                            ME_Cpx_mat_TE.(TE) = ME_Cpx_mat(TE_pos:nTE:end,:);
                                            ME_Opx_mat_TE.(TE) = ME_Opx_mat(TE_pos:nTE:end,:);
                                            ME_Grt_mat_TE.(TE)  = ME_Grt_mat(TE_pos:nTE:end,:);
                                            ME_Spl_mat_TE.(TE)  = ME_Spl_mat(TE_pos:nTE:end,:); 
                                            ME_Plg_mat_TE.(TE)  = ME_Plg_mat(TE_pos:nTE:end,:);  
                                        
                                            if TE_pos == 1
                                                MBC_Oli_aux = ME_Oli_mat_TE.(TE)(:,end);
                                                MBC_Cpx_aux = ME_Cpx_mat_TE.(TE)(:,end);
                                                MBC_Opx_aux = ME_Opx_mat_TE.(TE)(:,end);
                                                MBC_Grt_aux = ME_Grt_mat_TE.(TE)(:,end);
                                                MBC_Spl_aux = ME_Spl_mat_TE.(TE)(:,end);
                                                MBC_Plg_aux = ME_Plg_mat_TE.(TE)(:,end);
                                            else
                                                MBC_Oli_aux = cat(2,MBC_Oli_aux,ME_Oli_mat_TE.(TE)(:,end));
                                                MBC_Cpx_aux = cat(2,MBC_Cpx_aux,ME_Cpx_mat_TE.(TE)(:,end));
                                                MBC_Opx_aux = cat(2,MBC_Opx_aux,ME_Opx_mat_TE.(TE)(:,end));
                                                MBC_Grt_aux = cat(2,MBC_Grt_aux,ME_Grt_mat_TE.(TE)(:,end));
                                                MBC_Spl_aux = cat(2,MBC_Spl_aux,ME_Spl_mat_TE.(TE)(:,end));
                                                MBC_Plg_aux = cat(2,MBC_Plg_aux,ME_Plg_mat_TE.(TE)(:,end));
                                            end
    
                                        end
                                    
                                % Correct for solid without fluid without TE 
                                    MBC_Oli(ME_solid(:,3)>=1-tol_NTE,:)  = MBC_Oli_aux(ME_solid(:,3)>=1-tol_NTE,:);
                                    MBC_Cpx(ME_solid(:,3)>=1-tol_NTE,:) = MBC_Cpx_aux(ME_solid(:,3)>=1-tol_NTE,:);
                                    MBC_Opx(ME_solid(:,3)>=1-tol_NTE,:) = MBC_Opx_aux(ME_solid(:,3)>=1-tol_NTE,:);
                                    MBC_Grt(ME_solid(:,3)>=1-tol_NTE,:)  = MBC_Grt_aux(ME_solid(:,3)>=1-tol_NTE,:);
                                    MBC_Spl(ME_solid(:,3)>=1-tol_NTE,:)  = MBC_Spl_aux(ME_solid(:,3)>=1-tol_NTE,:);
                                    MBC_Plg(ME_solid(:,3)>=1-tol_NTE,:)  = MBC_Plg_aux(ME_solid(:,3)>=1-tol_NTE,:);
    
                                % Correct for solid with fluid without TE 
                                    MBC_Oli(ind_correct_TE_met_aux)  = MBC_Oli_aux(ind_correct_TE_met_aux);
                                    MBC_Cpx(ind_correct_TE_met_aux) = MBC_Cpx_aux(ind_correct_TE_met_aux);
                                    MBC_Opx(ind_correct_TE_met_aux) = MBC_Opx_aux(ind_correct_TE_met_aux);
                                    MBC_Grt(ind_correct_TE_met_aux)  = MBC_Grt_aux(ind_correct_TE_met_aux);
                                    MBC_Spl(ind_correct_TE_met_aux)  = MBC_Spl_aux(ind_correct_TE_met_aux);
                                    MBC_Plg(ind_correct_TE_met_aux)  = MBC_Plg_aux(ind_correct_TE_met_aux);
    
                                MRadi_plot   = MRadi;   MRadi_plot(MRadi_plot==0) = NaN;
                                MRadi_0_plot = MRadi_0; MRadi_0_plot(MRadi_0_plot==0) = NaN;
    
                                TE_Solid = TE_Oli+TE_Cpx+TE_Opx+TE_Grt+TE_Spl+TE_Plg; 
    
                                index_double = [index_double; (iter_TE-1)*size(TE_Solid,1)+ind_fluid];
    
                                mat_ME_Oli = cell2mat(ME_Oli); 
                                mat_ME_Cpx = cell2mat(ME_Cpx); 
                                mat_ME_Opx = cell2mat(ME_Opx); 
                                mat_ME_Grt = cell2mat(ME_Grt); 
                                mat_ME_Spl = cell2mat(ME_Spl); 
                                mat_ME_Plg = cell2mat(ME_Plg); 
    
                                Error_Oli_iter = (MBC_Oli  - reshape(mat_ME_Oli(:,end),nTE,[])')./MBC_Oli; Error_Oli_iter(MRadi(:,1)==0,:)=0;       Error_Oli_iter(correct_TE_met) = 0;    %Error_Oli_iter(abs(MBC_Oli  - reshape(mat_ME_Oli(:,end),nTE,[])')<tol_TE_abs) = 0;
                                Error_Cpx_iter = (MBC_Cpx  - reshape(mat_ME_Cpx(:,end),nTE,[])')./MBC_Cpx; Error_Cpx_iter(MRadi(:,2)==0,:)=0;  Error_Cpx_iter(correct_TE_met) = 0;   %Error_Cpx_iter(abs(MBC_Cpx  - reshape(mat_ME_Cpx(:,end),nTE,[])')<tol_TE_abs) = 0;
                                Error_Opx_iter = (MBC_Opx  - reshape(mat_ME_Opx(:,end),nTE,[])')./MBC_Opx; Error_Opx_iter(MRadi(:,3)==0,:)=0;  Error_Opx_iter(correct_TE_met) = 0;   %Error_Opx_iter(abs(MBC_Opx  - reshape(mat_ME_Opx(:,end),nTE,[])')<tol_TE_abs) = 0;
                                Error_Grt_iter = (MBC_Grt  - reshape(mat_ME_Grt(:,end),nTE,[])')./MBC_Grt; Error_Grt_iter(MRadi(:,4)==0,:)=0;       Error_Grt_iter(correct_TE_met) = 0;    %Error_Grt_iter(abs(MBC_Grt  - reshape(mat_ME_Grt(:,end),nTE,[])')<tol_TE_abs) = 0;
                                Error_Spl_iter = (MBC_Spl  - reshape(mat_ME_Spl(:,end),nTE,[])')./MBC_Spl; Error_Spl_iter(MRadi(:,5)==0,:)=0;       Error_Spl_iter(correct_TE_met) = 0;    %Error_Spl_iter(abs(MBC_Spl  - reshape(mat_ME_Spl(:,end),nTE,[])')<tol_TE_abs) = 0;
                                Error_Plg_iter = (MBC_Plg  - reshape(mat_ME_Plg(:,end),nTE,[])')./MBC_Plg; Error_Plg_iter(MRadi(:,6)==0,:)=0;       Error_Plg_iter(correct_TE_met) = 0;    %Error_Plg_iter(abs(MBC_Plg  - reshape(mat_ME_Plg(:,end),nTE,[])')<tol_TE_abs) = 0;
    
                                correct_TE_Oli_ind  = (abs(Error_Oli_iter)>tol_TE)  & (abs(MBC_Oli  - reshape(mat_ME_Oli(:,end),nTE,[])')>tol_TE_abs);
                                correct_TE_Cpx_ind = (abs(Error_Cpx_iter)>tol_TE) & (abs(MBC_Cpx - reshape(mat_ME_Cpx(:,end),nTE,[])')>tol_TE_abs);
                                correct_TE_Opx_ind = (abs(Error_Opx_iter)>tol_TE) & (abs(MBC_Opx - reshape(mat_ME_Opx(:,end),nTE,[])')>tol_TE_abs);
                                correct_TE_Grt_ind  = (abs(Error_Grt_iter)>tol_TE)  & (abs(MBC_Grt  - reshape(mat_ME_Grt(:,end),nTE,[])')>tol_TE_abs);
                                correct_TE_Spl_ind  = (abs(Error_Spl_iter)>tol_TE)  & (abs(MBC_Spl  - reshape(mat_ME_Spl(:,end),nTE,[])')>tol_TE_abs);
                                correct_TE_Plg_ind  = (abs(Error_Plg_iter)>tol_TE)  & (abs(MBC_Plg  - reshape(mat_ME_Plg(:,end),nTE,[])')>tol_TE_abs);
    
                                correct_TE_ind = correct_TE_Oli_ind | correct_TE_Cpx_ind | correct_TE_Opx_ind | correct_TE_Grt_ind | correct_TE_Spl_ind | correct_TE_Plg_ind; 
    
                                Kd_appa_Oli  = TE_Oli./ME_fluid_solid(:,4:end);  Kd_appa_Oli(isinf(Kd_appa_Oli))=NaN;     Kd_appa_Oli(ME_fluid_solid(:,4:end)<tol_TE)=NaN;
                                Kd_appa_Cpx = TE_Cpx./ME_fluid_solid(:,4:end); Kd_appa_Cpx(isinf(Kd_appa_Cpx))=NaN;   Kd_appa_Cpx(ME_fluid_solid(:,4:end)<tol_TE)=NaN;
                                Kd_appa_Opx = TE_Opx./ME_fluid_solid(:,4:end); Kd_appa_Opx(isinf(Kd_appa_Opx))=NaN;   Kd_appa_Opx(ME_fluid_solid(:,4:end)<tol_TE)=NaN;
                                Kd_appa_Grt  = TE_Grt./ME_fluid_solid(:,4:end);  Kd_appa_Grt(isinf(Kd_appa_Grt))=NaN;     Kd_appa_Grt(ME_fluid_solid(:,4:end)<tol_TE)=NaN;
                                Kd_appa_Spl  = TE_Spl./ME_fluid_solid(:,4:end);  Kd_appa_Spl(isinf(Kd_appa_Spl))=NaN;     Kd_appa_Spl(ME_fluid_solid(:,4:end)<tol_TE)=NaN;
                                Kd_appa_Plg  = TE_Plg./ME_fluid_solid(:,4:end);  Kd_appa_Plg(isinf(Kd_appa_Plg))=NaN;     Kd_appa_Plg(ME_fluid_solid(:,4:end)<tol_TE)=NaN;
    
                                Error_Oli  = [ Error_Oli  Error_Oli_iter];
                                Error_Cpx = [ Error_Cpx Error_Cpx_iter];
                                Error_Opx = [ Error_Opx Error_Opx_iter];
                                Error_Grt  = [ Error_Grt  Error_Grt_iter];
                                Error_Spl  = [ Error_Spl  Error_Spl_iter];
                                Error_Plg  = [ Error_Plg  Error_Plg_iter];
    
                                % Check convergence
                                    max_iter_TE = 40;
    
                                    if any(any(correct_TE_ind~=0))==0
                                            exit_TE = 1;
                                            fprintf('Time step %d (time %d [My]). Convergence reached at iteration %d.\n',iter,time/1e6/365/24/60/60,iter_TE);
                                    end
                                    if  iter_TE > max_iter_TE
                                        exit_TE = 1;
                                        fprintf('Time step %d (time %d [My]). Convergence NOT reached. Exit at iteration %d.\n',iter,time/1e6/365/24/60/60,iter_TE);
                                    end
                            end
    
                        for TE_pos = 1:nTE
                            [~,ia,~] = unique(ME_fluid_solid(:,2));
                            ME_fluid(:,3+TE_pos) = interp1(ME_fluid_solid(ia,2),ME_fluid_aux_iter(ia,TE_pos),ME_fluid(:,2),'nearest','extrap');
                        end 
    
                    % CLOSING
    
                        % Insert solid & fluid
                            [ME_solid,ME_fluid_solid,ME_solid_TP_w,ME_solid_TP_v,ME_solid_Rho_TP,ME_Oli,ME_Cpx,ME_Opx,ME_Grt,ME_Spl,ME_Plg,TE_Oli,TE_Cpx,TE_Opx,TE_Grt,TE_Spl,TE_Plg,MRadi] = ...
                                insertReactiveREE_1D_solid(ME_solid,ME_fluid_solid,ME_solid_TP_w,ME_solid_TP_v,ME_solid_Rho_TP,MRadi,ME_Oli,ME_Cpx,ME_Opx,ME_Grt,ME_Spl,ME_Plg,TE_Oli,TE_Cpx,TE_Opx,TE_Grt,TE_Spl,TE_Plg,input_TP_v,input_TP_w,gridxt,gridyt);
                            [ME_fluid,ME_fluid_TP_w,ME_fluid_TP_v,ME_fluid_Rho_TP] = ...
                                insertReactiveREE_1D_fluid(ME_fluid,ME_fluid_TP_w,ME_fluid_TP_v,ME_fluid_Rho_TP,gridxt,gridyt,TE_input_liquid,input_TP_v);
    
                        % mat_ME
                            mat_ME_Oli  = cell2mat(ME_Oli);  MBC_Oli_0  = reshape(mat_ME_Oli(:,end),nTE,[])';
                            mat_ME_Cpx = cell2mat(ME_Cpx); MBC_Cpx_0 = reshape(mat_ME_Cpx(:,end),nTE,[])';
                            mat_ME_Opx = cell2mat(ME_Opx); MBC_Opx_0 = reshape(mat_ME_Opx(:,end),nTE,[])';
                            mat_ME_Grt  = cell2mat(ME_Grt);  MBC_Grt_0  = reshape(mat_ME_Grt(:,end),nTE,[])';
                            mat_ME_Spl  = cell2mat(ME_Spl);  MBC_Spl_0  = reshape(mat_ME_Spl(:,end),nTE,[])';
                            mat_ME_Plg  = cell2mat(ME_Plg);  MBC_Plg_0  = reshape(mat_ME_Plg(:,end),nTE,[])';

            end

            %% COMPUTE BULK WITH NEW PARTICLES
                
                % mass TP / mass T
                    ME_solid_w_Oli   = repmat(ME_solid_TP_w(:,1),1,nTE);
                    ME_solid_w_Cpx  = repmat(ME_solid_TP_w(:,2),1,nTE);
                    ME_solid_w_Opx  = repmat(ME_solid_TP_w(:,3),1,nTE);
                    ME_solid_w_Grt   = repmat(ME_solid_TP_w(:,4),1,nTE);
                    ME_solid_w_Spl   = repmat(ME_solid_TP_w(:,5),1,nTE);
                    ME_solid_w_Plg   = repmat(ME_solid_TP_w(:,6),1,nTE);
                    ME_solid_w_Melt = repmat(ME_solid_TP_w(:,7),1,nTE);

                % solid [mass solid / mass T]
                    ME_solid_w_solid = ME_solid_w_Oli + ME_solid_w_Cpx + ME_solid_w_Opx + ME_solid_w_Grt + ME_solid_w_Spl + ME_solid_w_Plg;  % normalize every TP contribution
                
                % Factor [virtual total mass / total mass]
                    ME_solid_w_total = ME_solid_w_solid+ME_solid_w_Melt;

                % bulk TE [mass TE/ mass T]
                    ME_TE  =       (ME_solid_w_Oli.*TE_Oli ...
                                + ME_solid_w_Cpx.*TE_Cpx ...
                                + ME_solid_w_Opx.*TE_Opx ...
                                + ME_solid_w_Grt.*TE_Grt ...
                                + ME_solid_w_Spl.*TE_Spl ...
                                + ME_solid_w_Plg.*TE_Plg ...
                                + ME_solid_w_Melt.*ME_fluid_solid(:,4:3+nTE))./ME_solid_w_total;

                % mass TP / vol T
                    ME_solid_RhoPhi_Oli  = repmat(ME_solid_TP_v(:,1).*ME_solid_Rho_TP(:,1),1,nTE);
                    ME_solid_RhoPhi_Cpx = repmat(ME_solid_TP_v(:,2).*ME_solid_Rho_TP(:,2),1,nTE);
                    ME_solid_RhoPhi_Opx = repmat(ME_solid_TP_v(:,3).*ME_solid_Rho_TP(:,3),1,nTE);
                    ME_solid_RhoPhi_Grt  = repmat(ME_solid_TP_v(:,4).*ME_solid_Rho_TP(:,4),1,nTE);
                    ME_solid_RhoPhi_Spl  = repmat(ME_solid_TP_v(:,5).*ME_solid_Rho_TP(:,5),1,nTE);
                    ME_solid_RhoPhi_Plg  = repmat(ME_solid_TP_v(:,6).*ME_solid_Rho_TP(:,6),1,nTE);
                    ME_solid_RhoPhi_Melt  = repmat(ME_solid_TP_v(:,7).*ME_solid_Rho_TP(:,7),1,nTE);

                % mass T / vol T
                    ME_Rho_solid = ME_solid_RhoPhi_Oli(:,1) + ME_solid_RhoPhi_Cpx(:,1) + ME_solid_RhoPhi_Opx(:,1) + ME_solid_RhoPhi_Grt(:,1) + ME_solid_RhoPhi_Spl(:,1) + ME_solid_RhoPhi_Plg(:,1) + ME_solid_RhoPhi_Melt(:,1);

                % mass TP / vol T
                    ME_fluid_RhoPhi_Oli  = repmat(ME_fluid_TP_v(:,1).*ME_fluid_Rho_TP(:,1),1,nTE);
                    ME_fluid_RhoPhi_Cpx = repmat(ME_fluid_TP_v(:,2).*ME_fluid_Rho_TP(:,2),1,nTE);
                    ME_fluid_RhoPhi_Opx = repmat(ME_fluid_TP_v(:,3).*ME_fluid_Rho_TP(:,3),1,nTE);
                    ME_fluid_RhoPhi_Grt  = repmat(ME_fluid_TP_v(:,4).*ME_fluid_Rho_TP(:,4),1,nTE);
                    ME_fluid_RhoPhi_Spl  = repmat(ME_fluid_TP_v(:,5).*ME_fluid_Rho_TP(:,5),1,nTE);
                    ME_fluid_RhoPhi_Plg  = repmat(ME_fluid_TP_v(:,6).*ME_fluid_Rho_TP(:,6),1,nTE);
                    ME_fluid_RhoPhi_Melt  = repmat(ME_fluid_TP_v(:,7).*ME_fluid_Rho_TP(:,7),1,nTE);

                % mass T / vol T
                    ME_Rho_fluid = ME_fluid_RhoPhi_Oli(:,1) + ME_fluid_RhoPhi_Cpx(:,1) + ME_fluid_RhoPhi_Opx(:,1) + ME_fluid_RhoPhi_Grt(:,1) + ME_fluid_RhoPhi_Spl(:,1) + ME_fluid_RhoPhi_Plg(:,1) + ME_fluid_RhoPhi_Melt(:,1);

            %% SAVE & UPDATE TIME
            
                time = time + timestep;
                timeMy = time/1e6/365/24/60/60;
                timeMy_mat = [timeMy_mat timeMy];
                
                % Sorting
                    ME_solid_sorted = ME_solid(:,2:end);
                    ME_fluid_sorted = ME_fluid_solid(:,2:end);
                    TE_solid = ME_solid_sorted(:,3:2+nTE);
                    TE_fluid = ME_fluid_sorted(:,3:2+nTE);

                % Saving
                    iter_TE_save = single(iter_TE);
                    NT_save = single(NT);
                    NPhi_save = single(NPhi);
                    NRho_save = single(NRho);
                    NRho_tp_save = single(NRho_tp);
                    ME_solid_TP_v_save = single(ME_solid_TP_v);
                    ME_solid_TP_w_save = single(ME_solid_TP_w);
                    ME_solid_Rho_TP_save = single(ME_solid_Rho_TP);
                    ME_TE_save = single(ME_TE);
                    NV_mat_save = single(NV_mat);
                    XT_save = single(XT);
                    ME_solid_save = single(ME_solid);
                    ME_fluid_save = single(ME_fluid);
                    ME_fluid_solid_save = single(ME_fluid_solid);
                    MRadi_save = single(MRadi);
                    TE_Oli_save = single(TE_Oli);
                    TE_Cpx_save = single(TE_Cpx);
                    TE_Opx_save = single(TE_Opx);
                    TE_Grt_save = single(TE_Grt);
                    TE_Spl_save = single(TE_Spl);
                    TE_Plg_save = single(TE_Plg);
                    Gamma_Oli_save = single(Gamma_Oli);
                    Gamma_Cpx_save = single(Gamma_Cpx);
                    Gamma_Opx_save = single(Gamma_Opx);
                    Gamma_Grt_save = single(Gamma_Grt);
                    Gamma_Spl_save = single(Gamma_Spl);
                    Gamma_Plg_save = single(Gamma_Plg);
                    Error_Oli_save = single(Error_Oli);
                    Error_Cpx_save = single(Error_Cpx);
                    Error_Opx_save = single(Error_Opx);
                    Error_Grt_save = single(Error_Grt);
                    Error_Spl_save = single(Error_Spl);
                    Error_Plg_save = single(Error_Plg);
                    timestep_save = single(timestep);
                    time_save = single(time);
                    timeMy_save = single(timeMy);
                    timeMy_mat_save = single(timeMy_mat);
                    
                    iter_save = single(iter);
                    gridxt_save = single(gridxt);
                    gridyt_save = single(gridyt);

                    ME_Oli_save  = single(cell2mat(ME_Oli));
                    ME_Cpx_save = single(cell2mat(ME_Cpx));
                    ME_Opx_save = single(cell2mat(ME_Opx));
                    ME_Grt_save  = single(cell2mat(ME_Grt));
                    ME_Spl_save  = single(cell2mat(ME_Spl));
                    ME_Plg_save  = single(cell2mat(ME_Plg));

                    exp_factor = 0.3;
                    for index_mesh=1:nc; aux_mesh(index_mesh)=(1-(1/index_mesh)^exp_factor)/(1-(1/nc)^exp_factor); end
                    xmesh_Cpx = MRadi(1,2)*aux_mesh;

                    if ~exist('output')
                        mkdir output
                    end

                    if ismember(iter,save_list)
                        baseFileName = sprintf('savetime_%d', save_iter);
                        fullFileName = fullfile('output', baseFileName);
                        save(fullFileName,'save_list','xmesh_Cpx','ME_solid_TP_v_save','NV_mat_save','XT_save','ME_solid_save','ME_fluid_save','ME_fluid_solid_save','ME_Cpx_save','MRadi_save','timestep_save','time_save','timeMy_save','iter_save','gridxt_save','gridyt_save','timeMy_mat_save');
                        save_iter = save_iter+1;
                    end

                    ME_solid_sorted = ME_solid_save(:,2:end);
                    ME_fluid_sorted = ME_fluid_solid_save(:,2:end);

                    disp('Saving completed')
                    iter=iter+1;

        end