function [Gamma_Ol,Gamma_Cpx,Gamma_Opx,Gamma_Gt,Gamma_Sp,Gamma_Pl,ME_Ol,ME_Cpx,ME_Opx,ME_Gt,ME_Sp,ME_Pl,TE_Ol,TE_Cpx,TE_Opx,TE_Gt,TE_Sp,TE_Pl,MRadi,MRadi_0,ME_solid] = ...
    diffusion1D_6TP(ME_Ol_eq,ME_Cpx_eq,ME_Opx_eq,ME_Gt_eq,ME_Sp_eq,ME_Pl_eq,ME_melt_eq,BC_Ol_0,BC_Cpx_0,BC_Opx_0,BC_Gt_0,BC_Sp_0,BC_Pl_0,BC_Ol,BC_Cpx,BC_Opx,BC_Gt,BC_Sp,BC_Pl,ME_solid,ME_solid0,ME_solid_TP_v,ME_solid_TP_v0,ME_solid_TP_w,ME_Ol,ME_Cpx,ME_Opx,ME_Gt,ME_Sp,ME_Pl,ME_solid_Rho_TP,ME_solid_Rho_TP0,Mdiff_Ol,Mdiff_Cpx,Mdiff_Opx,Mdiff_Gt,Mdiff_Sp,Mdiff_Pl,timestep,nc,nREE,diff_nTP,diff_nREE,fix_radi,correct_TE_ind,correct_TE_met)

TP_list=lst.TP;
nTP=lst.nTP;
nTE=lst.nTE;

Mdiff.Oli = Mdiff_Ol;
Mdiff.Cpx = Mdiff_Cpx;
Mdiff.Opx = Mdiff_Opx;
Mdiff.Grt = Mdiff_Gt;
Mdiff.Spl = Mdiff_Sp;
Mdiff.Plg = Mdiff_Pl;

ME.Oli = ME_Ol;
ME.Cpx = ME_Cpx;
ME.Opx = ME_Opx;
ME.Grt = ME_Gt;
ME.Spl = ME_Sp;
ME.Plg = ME_Pl;

ME_eq.Oli = ME_Ol_eq;
ME_eq.Cpx = ME_Cpx_eq;
ME_eq.Opx = ME_Opx_eq;
ME_eq.Grt = ME_Gt_eq;
ME_eq.Spl = ME_Sp_eq;
ME_eq.Plg = ME_Pl_eq;

BC.Oli = BC_Ol;
BC.Cpx = BC_Cpx;
BC.Opx = BC_Opx;
BC.Grt = BC_Gt;
BC.Spl = BC_Sp;
BC.Plg = BC_Pl;

BC_0.Oli = BC_Ol_0;
BC_0.Cpx = BC_Cpx_0;
BC_0.Opx = BC_Opx_0;
BC_0.Grt = BC_Gt_0;
BC_0.Spl = BC_Sp_0;
BC_0.Plg = BC_Pl_0;

% Compute diffusion profiles inside thermodynamic phases (TP) for different
% oxides/trace elements (b)
% Compute the grain boundary gradient
%
% Oxide/trace elements are computed in %wt of solid (s). Thus:
%   ->           \sum_{b} C_s^{b,TP} \neq 1  (sum of oxides concentrations inside a TP is not = 1)
%   -> \sum_{TP} \sum_{b} C_s^{b,TP}    = 1  (sum of oxides concentrations over all TP is     = 1)
% Not sure about this. Need to check

% Compute boundary conditions for each oxide (b) inside each thermodynamic phaseare obtained
% BC_Ol  = zeros(size(ME_solid,1),5); % 5 is probably nTP
% 
% for index = 1:5
%     F_Ol  = scatteredInterpolant(XT(:,1),XT(:,2),NOl(:,index),'linear','nearest');
%     BC_Ol(:,index)  = F_Ol([ME_solid(:,1),ME_solid(:,2)]);
% end

% compute initial radious
% 
%         /             \^(1/3)
%        |     X_tp      |
% Radi = |---------------|
%        | (4/3)*pi*n_tp |
%         \             /
%
% X_tp as volume fraction of tp within solid
%

% Normalize, ME_solid_TP_w 
ME_solid_TP_v = ME_solid_TP_v./repmat(sum(ME_solid_TP_v,2),1,size(ME_solid_TP_v,2));

% Compute Radios
for i=1:nTP
    TP_str=TP_list(i);
    Rradi.(TP_str) = ((ME_solid_TP_v(:,i))./((4/3)*pi.*ME_solid(:,3+nTE+i))).^(1/3);
    if diff_nTP(i)==0; Rradi.(TP_str) = 0*Rradi.(TP_str); end
    if fix_radi(i)~=0; Rradi.(TP_str) = fix_radi(i)*ones(size(Rradi.(TP_str))); 
        ME_solid(:,3+nTE+i) = ME_solid_TP_v(:,i)./((4/3)*pi.*Rradi.(TP_str).^3);
    end
end

% Rradi.(TP_str) = 0.005*ones(size(Rradi.(TP_str)));

MRadi=zeros(size(ME_solid, 1), nTP);
for i = 1:nTP
    TP_str = TP_list(i);
    MRadi(:,i) = Rradi.(TP_str);
end

% Normalize, ME_solid_TP_w 
ME_solid_TP_v0 = ME_solid_TP_v0./repmat(sum(ME_solid_TP_v0,2),1,size(ME_solid_TP_v0,2));

% Compute Initial Radios
for i = 1:nTP
    TP_str = TP_list(i);
    Rradi_0.(TP_str)  = ((ME_solid_TP_v0(:,i))./((4/3)*pi.*ME_solid0(:,3+nTE+i))).^(1/3);
    if diff_nTP(i)==0;  Rradi_0.(TP_str)  = 0*Rradi_0.(TP_str); end
    if fix_radi(i)~=0; Rradi_0.(TP_str)   = fix_radi(i)*ones(size(Rradi_0.(TP_str))); 
        ME_solid0(:,3+nTE+i) = ME_solid_TP_v0(:,i)./((4/3)*pi.*Rradi_0.(TP_str).^3);
    end
end

MRadi_0=zeros(size(MRadi));
for i = 1:nTP
    TP_str = TP_list(i);
    MRadi_0(:,i) = Rradi_0.(TP_str);
    Gamma.(TP_str) = mat2cell(zeros(size(ME_solid,1),nTE),ones(size(ME_solid,1),1),nTE);
    TE.(TP_str) = Gamma.(TP_str);
end

% we approximate the initial profiles with the previous one but with the mesh changed
% loop over all the particles to compute 1D profiles [#TP x #oxides/trace elements]

% transient time stepping. Number of timesteps
nt = 2;
n1 = nt-1;
aux_BC = repmat(0:n1,nTE,1);
ind_TE = find(diff_nTE==1);
index = any(correct_TE_ind,2);
index_met = any(correct_TE_met,2); 
% index_met = index;
index_keep = any(correct_TE_ind==0,2); % WARNING unused
ind_Particles = 1:size(ME_solid,1);
vec_particles     = ind_Particles(index);
vec_particles_met = ind_Particles(index_met);
% index = any(correct_TE_ind,2)
% parfor nparticles = 1:size(MR_solid,1)
% parfor nparticles = 1:size(ME_solid,1);
% parfor nparticles_index = 1:length(vec_particles)

for i = 1:nTP
    TP_str = TP_list(i);
    BC_0_parfor.(TP_str) = BC_0.(TP_str)(index,:);
    BC_parfor.(TP_str)   = BC.(TP_str)(index,:);
    ME_parfor.(TP_str)   = ME.(TP_str)(index);
    Gamma_parfor.(TP_str) = Gamma.(TP_str)(index);
    Rradi_parfor.(TP_str) = Rradi.(TP_str)(index);
    Rradi_0_parfor.(TP_str) = Rradi_0.(TP_str)(index);
    Mdiff_parfor.(TP_str) = Mdiff.(TP_str)(index,:);
    TE_parfor.(TP_str) = TE.(TP_str)(index);
end

ME_solid_parfor  = ME_solid(index,:);
ME_solid0_parfor = ME_solid0(index,:);
ME_solid_Rho_TP_parfor  = ME_solid_Rho_TP(index,:);
ME_solid_Rho_TP0_parfor = ME_solid_Rho_TP0(index,:);
ME_solid_TP_v_parfor = ME_solid_TP_v(index,:); % WARNING unused
ME_solid_TP_v0_parfor = ME_solid_TP_v0(index,:); % WARNING unused

for i = 1:nTP
    TP_str = TP_list(i);
    % ME_parfor_met.(TP_str)   = ME.(TP_str)(index_met); WARNING check why
    ME_parfor_met.(TP_str)   = ME_eq.(TP_str)(index_met);
    Gamma_parfor_met.(TP_str) = Gamma.(TP_str)(index_met);
    Rradi_parfor_met.(TP_str) = Rradi.(TP_str)(index_met);
    TE_parfor_met.(TP_str) = TE.(TP_str)(index_met);
end

exp_factor = 0.3; % WARNING to be passed to settings
% exp_factor = -1;
for index_mesh=1:nc
    aux_mesh(index_mesh)=(1-(1/index_mesh)^exp_factor)/(1-(1/nc)^exp_factor); % WARNING to be initialized
end

tol_TE = 1e-15; % WARNING to be passed to settings

%parfor nparticles = 1:length(vec_particles) 
for nparticles = 1:length(vec_particles) % WARNING parfor loop to be benchmarked

    time_discrete = [0:timestep/nt:timestep];
    for i = 1:nTP
        TP_str = TP_list(i);
        if diff_nTP(i)==1
            BC_aux.(TP_str)  = repmat(BC_0_parfor.(TP_str)(nparticles,:)',1,nt)  + aux_BC.*repmat((BC_parfor.(TP_str)(nparticles,:)'  - BC_0_parfor.(TP_str)(nparticles,:)')/n1,1,nt); % WARNING unused
        end
    end
    
    % Solve
    for index_time = 1:(nt-1)
        % we need 3 steps for the solver
        tspan = [time_discrete(nt) 0.5*(time_discrete(nt+1)+time_discrete(nt)) time_discrete(nt+1)]; % WARNING unused
        % update initial condition for transient problem
        for i = 1:nTP
            TP_str = TP_list(i);
            MRini_0.(TP_str)  = ME_parfor.(TP_str){nparticles};
        end

        % update initial condition for PC problem
        for i = 1:nTP
            TP_str = TP_list(i);
            MRini_PC.(TP_str)  = ME_parfor.(TP_str){nparticles};
        end
        
        % Matrix Assembly
        % Espatial discretization
        for i = 1:nTP
            TP_str = TP_list(i);
            x1.(TP_str)    = Rradi_parfor.(TP_str)(nparticles); 
            xmesh.(TP_str)    = aux_mesh*x1.(TP_str);
            if x1.(TP_str)==0;    xmesh.(TP_str)=zeros(1,nc);    end
            x1_0.(TP_str)  = Rradi_0_parfor.(TP_str)(nparticles);
            xmesh_0.(TP_str)  = aux_mesh*x1_0.(TP_str);
            if x1_0.(TP_str)==0
                if x1.(TP_str)==0
                    xmesh_0.(TP_str)=zeros(1,nc);
                else
                    xmesh_0.(TP_str)  = xmesh.(TP_str);
                end
            end
        end
       
        xipg = [-1/sqrt(3) 1/sqrt(3)]';
        wpg = [1 1]';
        
        % Shape functions and its derivatives on the reference element
        N_mef   =  [(1-xipg)/2 (1+xipg)/2];
        Nxi_mef =  [-1/2 1/2; -1/2 1/2];

        % Matrices obtained by discretizing a convection-diffusion equation
        for i = 1:nTP
            TP_str = TP_list(i);
            if diff_nTP(i)==1
                [M.(TP_str),K.(TP_str),~]   = matrices_1D(xmesh.(TP_str),xipg,wpg,N_mef,Nxi_mef);
                M.(TP_str) = spdiags(sum(M.(TP_str),2), 0,length(M.(TP_str)),length(M.(TP_str)));
            end
        end

        for j = 1:length(ind_TE)

            % if Rradi_Ol(vec_particlesnparticles)~=0;  ME_Ol{vec_particlesnparticles}(ind_TE(index),:)  = diffusion_1D_PDE(Rradi_Ol(vec_particlesnparticles),MRini_Ol_0(ind_TE(index),:),Mdiff_Ol(vec_particlesnparticles,ind_TE(index)),BC_Ol_aux(ind_TE(index),index_time+1),tspan,nc); end           
            % if Rradi_Ol(vec_particlesnparticles)~=0;  ME_Ol{vec_particlesnparticles}(index,:)  = diffusion_1D_PDE(Rradi_Ol(vec_particlesnparticles),MRini_Ol_0(index,:),Mdiff_Ol(vec_particlesnparticles,index),BC_Ol_aux(index,index_time+1),tspan,nc); end
                        
            % W = 1/2;
            % w = 1; 
            % % Boundary conditions (lagrange multipliers method)
            % Accd = zeros(1,nc);
            % Accd(1,nc)=1; 
            % Accd = sparse(Accd);
            % bccd = [BC_Ol_aux(ind_TE(index),index_time+1)];
            % M_Ol  = spdiags(sum(M_Ol,2), 0,length(M_Ol),length(M_Ol));
            % Sol = Galerkin(W,w,0,Mdiff_Ol(vec_particles_parfor(nparticles),ind_TE(index)),0*xmesh_Ol',K_Ol,M_Ol,0,xmesh_Ol,timestep,1,MRini_Ol_0(ind_TE(index),:)',Accd,bccd*5);
            % keyboard
            % if vec_particles(nparticles) == 1188; keyboard; end

            % Phase change Part
            for i=1:nTP
                TP_str=TP_list(i);
                if x1.(TP_str) > x1_0.(TP_str)
                    MRini_PC.(TP_str)(j,:)  = interp1(xmesh_0.(TP_str),MRini_0.(TP_str)(j,:),xmesh.(TP_str),'linear');
                    MRini_PC.(TP_str)(j,xmesh.(TP_str)>=x1_0.(TP_str))    = BC_aux.(TP_str)(j,index_time+1);
                elseif  x1.(TP_str) < x1_0.(TP_str)
                    MRini_PC.(TP_str)(j,:)  = interp1(xmesh_0.(TP_str), MRini_0.(TP_str)(j,:),xmesh.(TP_str),'linear');
                end
            end

            % Diffusion Part
                for i=1:nTP
                    TP_str=TP_list(i);
                    % if Rradi_Ol_parfor(nparticles)~=0;
                    %     ME_Ol_parfor{nparticles}(j,:)  = solveWithLagrangeMultipliers(M_Ol + timestep*0.5*K_Ol*Mdiff_Ol_parfor(nparticles,j) , M_Ol*MRini_Ol_PC(j,:)' + 0.5*timestep*K_Ol*Mdiff_Ol_parfor(nparticles,j)*MRini_Ol_PC(j,:)', [nc BC_Ol_aux(j,index_time+1)]);
                    % end;
                    % if all(abs(ME_Ol_parfor{nparticles}(j,1:end-1)- MRini_Ol_0(j,1:end-1))./ ME_Ol_parfor{nparticles}(j,1:end-1)<tol_TE);
                    %     ME_Ol_parfor{nparticles}(j,:)  = MRini_Ol_0(j,:);
                    % end;
                    if Rradi_parfor.(TP_str)(nparticles)~=0
                        ME_parfor.(TP_str){nparticles}(j,:) = solveWithLagrangeMultipliers(M.(TP_str) + timestep*K.(TP_str)*Mdiff_parfor.(TP_str)(nparticles,j) , M.(TP_str)*MRini_PC.(TP_str)(j,:)', [nc BC_aux.(TP_str)(j,index_time+1)]);
                    end
                    if all(abs(ME_parfor.(TP_str){nparticles}(j,1:end-1)- MRini_0.(TP_str)(j,1:end-1))./ ME_parfor.(TP_str){nparticles}(j,1:end-1)<tol_TE)
                        ME_parfor.(TP_str){nparticles}(j,:) = MRini_0.(TP_str)(j,:);
                    end
                end

            % Mass balance
                for i=1:nTP
                    TP_str=TP_list(i);

                    integral0 = trapz(xmesh_0.(TP_str),4*pi*xmesh_0.(TP_str).^2.*MRini_0.(TP_str)(j,:));
                    integral1 = trapz(xmesh.(TP_str),4*pi*xmesh.(TP_str).^2.*ME_parfor.(TP_str){nparticles}(j,:));

                    % Gamma_Ol_parfor{nparticles}(j)   = (ME_solid_parfor(nparticles,3+nTE+1).*ME_solid_Rho_TP_parfor(nparticles,1).*integral1-ME_solid0_parfor(nparticles,3+nTE+1).*ME_solid_Rho_TP0_parfor(nparticles,1).*integral0)/timestep;
                    Gamma_parfor.(TP_str){nparticles}(j)   =  (integral1-integral0).*ME_solid0_parfor(nparticles,3+nTE+1).*ME_solid_Rho_TP0_parfor(nparticles,1)/timestep;
                    Gamma_parfor.(TP_str){nparticles}(j)   =  (integral1.*ME_solid_parfor(nparticles,3+nTE+1).*ME_solid_Rho_TP_parfor(nparticles,1)-integral0.*ME_solid0_parfor(nparticles,3+nTE+1).*ME_solid_Rho_TP0_parfor(nparticles,1))/timestep;
                    % Gamma_Ol_parfor{nparticles}(j)   =  ((integral1/trapz(xmesh_Ol,4*pi*xmesh_Ol.^2))*ME_solid_TP_v_parfor(nparticles,1)*ME_solid_Rho_TP_parfor(nparticles,1) - (integral0/trapz(xmesh_Ol_0,4*pi*xmesh_Ol_0.^2))*ME_solid_TP_v0_parfor(nparticles,1)*ME_solid_Rho_TP0_parfor(nparticles,1))/timestep;

                    % TE_Ol_parfor{nparticles}(index)  = ME_solid_parfor(nparticles,12).*integral1;   if xmesh_Ol(end)==0;  TE_Ol_parfor{nparticles}(index)=0;  end
                    % TE_Ol_parfor{nparticles}(ind_TE(index))  = integral1/(4/3*pi*xmesh_Ol(end).^3)/(1-ME_solid_TP_v_parfor(nparticles,7));   if xmesh_Ol(end)==0;  TE_Ol_parfor{nparticles}(ind_TE(index))=0;  end
                    
                    TE_parfor.(TP_str){nparticles}(j) = integral1/trapz(xmesh.(TP_str),4*pi*xmesh.(TP_str).^2);
                    if xmesh.(TP_str)(end)==0;  TE_parfor.(TP_str){nparticles}(j)=0;  end
                    % if Mdiff_Ol_parfor(nparticles,ind_TE(index))==0;   if integral1~=0; ME_Ol_parfor{nparticles}(ind_TE(index),:)  = integral0/integral1*ME_Ol_parfor{nparticles}(ind_TE(index),:);  end; end
                    
                end

            %             % Without n_tp
            %             integral0 = trapz(xmesh_Ol_0,4*pi*xmesh_Ol_0.^2.*MRini_Ol_0(index,:));      integral1 = trapz(xmesh_Ol,4*pi*xmesh_Ol.^2.*ME_Ol_parfor{nparticles}(index,:));      Gamma_Ol_parfor{nparticles}(index)   = (ME_solid_Rho_TP_parfor(nparticles,1).*integral1/(4/3*pi*xmesh_Ol(end).^3)-ME_solid_Rho_TP0_parfor(nparticles,1).*integral0/(4/3*pi*xmesh_Ol_0(end).^3))/timestep;
            %             TE_Ol_parfor{nparticles}(index)  = integral1/(4/3*pi*xmesh_Ol(end).^3);   if xmesh_Ol(end)==0;  TE_Ol_parfor{nparticles}(index)=0;  end
            %             if Mdiff_Ol_parfor(nparticles,index)==0;   ME_Ol_parfor{nparticles}(index,:)  = integral0/integral1*ME_Ol_parfor{nparticles}(index,:);  end

        end
        
        %         if ME_solid_parfor(nparticles,3)>1-1E-4
        %             Gamma_Ol_parfor{nparticles}=zeros(1,8);
        %         end

    end
    
    for i=1:nTP
        TP_str=TP_list(i);
        Grad_oxides_parfor.(TP_str)(nparticles,:) = (ME_parfor.(TP_str){nparticles}(:,end-1)-ME_parfor.(TP_str){nparticles}(:,end-1))/(Rradi_parfor.(TP_str)(nparticles)/(nt-1)); %WARNING unused
    end
    
end

%parfor nparticles = 1:length(vec_particles_met) 
for nparticles = 1:length(vec_particles_met) % WARNING parfor loop to be benchmarked
    
    % Spatial discretization
        for i=1:nTP
            TP_str=TP_list(i);
            x1.(TP_str) = Rradi_parfor_met.(TP_str)(nparticles);
            xmesh.(TP_str) = aux_mesh*x1.(TP_str);
            if x1.(TP_str)==0; xmesh.(TP_str)=zeros(1,nc); end
        end
    
    for j = 1:nTE
        
        % Mass balance
        for i=1:nTP
            TP_str=TP_list(i);
            integral1 = trapz(xmesh.(TP_str),4*pi*xmesh.(TP_str).^2.*ME_parfor_met.(TP_str){nparticles}(j,:));
            Gamma_parfor_met.(TP_str){nparticles}(j)   = 0;
            TE_parfor_met.(TP_str){nparticles}(j) = integral1/trapz(xmesh.(TP_str),4*pi*xmesh.(TP_str).^2);
            if xmesh.(TP_str)(end)==0
                TE_parfor_met.(TP_str){nparticles}(j)=0;  ME_parfor_met.(TP_str){nparticles}(j,:) = 0*ME_parfor_met.(TP_str){nparticles}(j,:);
            end
        end

    end
    
end
for i=1:nTP
    TP_str=TP_list(i);
    ME.(TP_str)(index) = ME_parfor.(TP_str);
    ME.(TP_str)(index_met)    = ME_parfor_met.(TP_str);
    Rradi.(TP_str)(index)     = Rradi_parfor.(TP_str);
    Rradi.(TP_str)(index_met) = Rradi_parfor_met.(TP_str);
    Mdiff.(TP_str)(index,:)   = Mdiff_parfor.(TP_str);
    TE.(TP_str)(index)       = TE_parfor.(TP_str);
    TE.(TP_str)(index_met)   = TE_parfor_met.(TP_str);
    Gamma.(TP_str)(index)     = Gamma_parfor.(TP_str);
    Gamma.(TP_str)(index_met) = Gamma_parfor_met.(TP_str);
end

for i=1:nTP
    TP_str=TP_list(i);
    ME_solid_w.(TP_str)   = repmat(ME_solid_TP_w(:,i),1,nTE);
end
ME_solid_w.melt = repmat(ME_solid_TP_w(:,nTP+1),1,nTE);

zero_aux = zeros(size(ME_solid,1),nTE);

for i=1:nTP
    TP_str=TP_list(i);
    % WARNING pass to settings
    Gamma.(TP_str) = cell2mat(Gamma.(TP_str));
    Gamma.(TP_str)(ME_solid(:,3)>1-2E-8,:)  = zero_aux(ME_solid(:,3)>1-2E-8,:);
    Gamma.(TP_str)(isnan(Gamma.(TP_str)))=0;
end

% TE mass / TP mass
for i=1:nTP
    TP_str=TP_list(i);
    TE.(TP_str) = cell2mat(TE.(TP_str));
end

% Factor TP mass / solid mass
ME_solid_w_solid=zeros(size(ME_solid,1),nTE);
for i=1:nTP
    TP_str=TP_list(i);
    ME_solid_w_solid = ME_solid_w_solid + ME_solid_w.(TP_str);
end

% TE mass / solid mass
ME_solid(:, 4:3+nTE) = zeros(size(ME_solid,1), nTE);
for i=1:nTP
    TP_str=TP_list(i);
    ME_solid(:,4:3+nTE) = ME_solid(:,4:3+nTE) + TE.(TP_str).*ME_solid_w.(TP_str)./ME_solid_w_solid; % TE mass / solid mass
end

ME_Ol  = ME.Oli;
ME_Cpx = ME.Cpx;
ME_Opx = ME.Opx;
ME_Gt  = ME.Grt;
ME_Sp  = ME.Spl;
ME_Pl  = ME.Plg;

TE_Ol  = TE.Oli;
TE_Cpx = TE.Cpx;
TE_Opx = TE.Opx;
TE_Gt  = TE.Grt;
TE_Sp  = TE.Spl;
TE_Pl  = TE.Plg;

Gamma_Ol  = Gamma.Oli;
Gamma_Cpx = Gamma.Cpx;
Gamma_Opx = Gamma.Opx;
Gamma_Gt  = Gamma.Grt;
Gamma_Sp  = Gamma.Spl;
Gamma_Pl  = Gamma.Plg;

end