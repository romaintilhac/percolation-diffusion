function [Gamma_TP, ME_TP, TE_TP, MRadi, MRadi_0, ME_solid] = ...
    diffusion1D_flexible(ME_TP_eq, ME_melt_eq, BC_TP_0, BC_TP, ...
    ME_solid, ME_solid0, ME_solid_TP_v, ME_solid_TP_v0, ...
    ME_solid_TP_w, ME_TP, ME_solid_Rho_TP, ME_solid_Rho_TP0, ...
    Mdiff_TP, correct_TE_ind, correct_TE_met, diff_nTP, diff_nTE, ...
    fix_radi, nTE, nc, timestep)
% Flexible, robust, phase-general 1D diffusion with defensive initialization.
% ME_TP_eq and ME_TP_0 must be cell arrays of numeric matrices (nParticles x nTE [x nMesh])
% If ME_TP_0 is not provided, equilibrium values will be used.

nTP = numel(ME_TP_eq); % Number of solid thermodynamic phases (melt is excluded)
nParticles = size(ME_solid,1);

% Normalize, ME_solid_TP_v
ME_solid_TP_v = ME_solid_TP_v ./ sum(ME_solid_TP_v,2);

% Compute Radii (current)
Rradi_TP = zeros(size(ME_solid_TP_v));
for i=1:nTP
    Rradi_TP(:,i) = (ME_solid_TP_v(:,i) ./ ((4/3)*pi.*ME_solid(:,3+nTE+i))).^(1/3);
    if diff_nTP(i)==0
        Rradi_TP(:,i) = 0*Rradi_TP(:,i);
    end
    if fix_radi(i)~=0
        Rradi_TP(:,i) = fix_radi(i)*ones(size(Rradi_TP(:,i)));
        ME_solid(:,3+nTE+i) = ME_solid_TP_v(:,i)./((4/3)*pi.*Rradi_TP(:,i).^3);
    end
end
MRadi = Rradi_TP;

% Normalize, ME_solid_TP_v0
ME_solid_TP_v0 = ME_solid_TP_v0 ./ sum(ME_solid_TP_v0,2);

% Compute Radii (initial)
Rradi_TP_0 = zeros(size(ME_solid_TP_v0));
for i=1:nTP
    Rradi_TP_0(:,i) = (ME_solid_TP_v0(:,i)./((4/3)*pi.*ME_solid0(:,3+nTE+i))).^(1/3);
    if diff_nTP(i)==0
        Rradi_TP_0(:,i) = 0*Rradi_TP_0(:,i);
    end
    if fix_radi(i)~=0
        Rradi_TP_0(:,i) = fix_radi(i)*ones(size(Rradi_TP_0(:,i)));
        ME_solid0(:,3+nTE+i) = ME_solid_TP_v0(:,i)./((4/3)*pi.*Rradi_TP_0(:,i).^3);
    end
end
MRadi_0 = Rradi_TP_0;


for i=1:nTP
    Gamma_TP{i} = mat2cell(zeros(nParticles,nTE),ones(nParticles,1),nTE);
    TE_TP{i} = Gamma_TP{i};
end

nt = 2;
n1 = nt-1;
aux_BC = repmat(0:n1,nTE,1);
ind_TE = find(diff_nTE==1);
index = any(correct_TE_ind,2);
index_met = any(correct_TE_met,2);
ind_Particles = 1:nParticles;
vec_particles     = ind_Particles(index);
vec_particles_met = ind_Particles(index_met);

% % Precut BC, ME, Gamma, Rradi, Mdiff, TE, Rho, TP_v for loop
% BC_TP_0_parfor = cell(1,nTP); BC_TP_parfor = cell(1,nTP);
% ME_TP_parfor = cell(1,nTP); Gamma_TP_parfor = cell(1,nTP); Rradi_TP_parfor = cell(1,nTP);
% Rradi_TP_0_parfor = cell(1,nTP); Mdiff_TP_parfor = cell(1,nTP); TE_TP_parfor = cell(1,nTP);

for i=1:nTP
    BC_TP_0_parfor{i} = BC_TP_0{i}(index,:);
    BC_TP_parfor{i} = BC_TP{i}(index,:);
    ME_TP_parfor{i} = ME_TP{i}(index);
    Gamma_TP_parfor{i} = Gamma_TP{i}(index);
    Rradi_TP_parfor{i} = Rradi_TP(:,i);
    Rradi_TP_0_parfor{i} = Rradi_TP_0(:,i);
    TE_TP_parfor{i} = TE_TP{i}(index);
    Mdiff_TP_parfor{i} = squeeze(Mdiff_TP(:,:,i)); %from 3d matrix Mdiff_TP
end

ME_solid_parfor  = ME_solid(index,:);
ME_solid0_parfor = ME_solid0(index,:);
ME_solid_Rho_TP_parfor  = ME_solid_Rho_TP(index,:);
ME_solid_Rho_TP0_parfor = ME_solid_Rho_TP0(index,:);
ME_solid_TP_v_parfor = ME_solid_TP_v(index,:);
ME_solid_TP_v0_parfor = ME_solid_TP_v0(index,:);

exp_factor = 0.3;
aux_mesh = zeros(1,nc);
for index_mesh=1:nc
    aux_mesh(index_mesh) = (1-(1/index_mesh)^exp_factor)/(1-(1/nc)^exp_factor);
end

tol_TE = 1e-15;

%% Parfor loop over particules

parfor np = 1:length(vec_particles)

    nparticles = vec_particles(np);
    time_discrete = [0:timestep/nt:timestep];
    BC_TP_aux = cell(1,nTP);
    for i=1:nTP
        if diff_nTP(i)==1
            BC_TP_aux{i} = repmat(BC_TP_0_parfor{i}(np,:)',1,nt) + ...
                aux_BC.*repmat((BC_TP_parfor{i}(np,:)' - BC_TP_0_parfor{i}(np,:)')/n1,1,nt);
        end
    end

    % Solve
    for index_time = 1:(nt-1)
        tspan = [time_discrete(nt) 0.5*(time_discrete(nt+1)+time_discrete(nt)) time_discrete(nt+1)];
        MRini_TP_0 = cell(1,nTP);
        MRini_TP_PC = cell(1,nTP);
        for i=1:nTP
            MRini_TP_0{i} = ME_TP_parfor{i}{np};
            MRini_TP_PC{i} = MRini_TP_0{i};
        end

        % Matrix assembly and spatial discretization
        xmesh_TP = cell(1,nTP); xmesh_TP_0 = cell(1,nTP);
        for i=1:nTP
            x1 = Rradi_TP_parfor{i}(np);
            xmesh_TP{i} = aux_mesh*x1;
            if x1==0
                xmesh_TP{i}=zeros(1,nc);
            end
            x1_0 = Rradi_TP_0_parfor{i}(np);
            xmesh_TP_0{i} = aux_mesh*x1_0;
            if x1_0==0
                if x1==0
                    xmesh_TP_0{i}=zeros(1,nc);
                else
                    xmesh_TP_0{i} = xmesh_TP{i};
                end
            end
        end
        xipg = [-1/sqrt(3) 1/sqrt(3)]';
        wpg = [1 1]';
        N_mef   =  [(1-xipg)/2 (1+xipg)/2];
        Nxi_mef =  [-1/2 1/2; -1/2 1/2];

        M_TP = cell(1,nTP); K_TP = cell(1,nTP);
        for i=1:nTP
            if diff_nTP(i)==1
                [Mtmp,Ktmp,~] = matrices_1D(xmesh_TP{i},xipg,wpg,N_mef,Nxi_mef);
                M_TP{i} = spdiags(sum(Mtmp,2), 0, length(Mtmp), length(Mtmp));
                K_TP{i} = Ktmp;
            end
        end

        for index_TE = 1:length(ind_TE)
            % Phase change interpolation
            for i=1:nTP
                x1 = Rradi_TP_parfor{i}(np);
                x1_0 = Rradi_TP_0_parfor{i}(np);
                if x1 > x1_0 && ~isempty(MRini_TP_0{i}) && size(MRini_TP_0{i},1)>=max(ind_TE(index_TE))
                    MRini_TP_PC{i}(ind_TE(index_TE),:) = interp1(xmesh_TP_0{i},MRini_TP_0{i}(ind_TE(index_TE),:),xmesh_TP{i},'linear');
                    MRini_TP_PC{i}(ind_TE(index_TE),xmesh_TP{i}>=x1_0) = MRini_TP_0{i}(ind_TE(index_TE),end);
                end
            end
            % Diffusion part
            for i=1:nTP
                if diff_nTP(i)==1 && Rradi_TP_parfor{i}(np)~=0 && ...
                   ~isempty(MRini_TP_PC{i}) && size(MRini_TP_PC{i},1)>=max(ind_TE(index_TE))
                    % Prepare boundary conditions for this TP
                    C = [1 BC_TP_aux{i}(ind_TE(index_TE),index_time+1); nc BC_TP_aux{i}(ind_TE(index_TE),index_time+1)];
                    % Fix: pass per-TE diffusion coefficients
                    D_TE = Mdiff_TP_parfor{i}(np, ind_TE(index_TE));
                    ME_TP_parfor{i}{np}(ind_TE(index_TE),:) = solveWithLagrangeMultipliers(...
                        M_TP{i} + timestep*K_TP{i}*D_TE, ...
                        M_TP{i} * MRini_TP_PC{i}(ind_TE(index_TE),:)', ...
                        C, 0); % WARNING was originally transposed by copilot
                else
                    if isempty(ME_TP_parfor{i}{np})
                        ME_TP_parfor{i}{np} = zeros(length(ind_TE),nc);
                    elseif size(ME_TP_parfor{i}{np},1)<max(ind_TE)
                        ME_TP_parfor{i}{np}(max(ind_TE),nc) = 0;
                    end
                end
            end
            % Mass balance and TE calculation
            for i=1:nTP
                xmesh = xmesh_TP{i};
                xmesh_0 = xmesh_TP_0{i};
                MRini = MRini_TP_0{i};
                ME_now = ME_TP_parfor{i}{np};
                if ~isempty(MRini) && size(MRini,1)>=max(ind_TE(index_TE)) && ~isempty(ME_now) && size(ME_now,1)>=max(ind_TE(index_TE))
                    integral0 = trapz(xmesh_0,4*pi*xmesh_0.^2.*MRini(ind_TE(index_TE),:));
                    integral1 = trapz(xmesh,4*pi*xmesh.^2.*ME_now(ind_TE(index_TE),:));
                    % Gamma
                    Gamma_TP_parfor{i}{np}(ind_TE(index_TE)) = ...
                        (integral1.*ME_solid_parfor(np,3+nTE+i).*ME_solid_Rho_TP_parfor(np,i) - ...
                        integral0.*ME_solid0_parfor(np,3+nTE+i).*ME_solid_Rho_TP0_parfor(np,i))/timestep;
                    % TE
                    TE_TP_parfor{i}{np}(ind_TE(index_TE)) = integral1/trapz(xmesh,4*pi*xmesh.^2);
                    if xmesh(end)==0
                        TE_TP_parfor{i}{np}(ind_TE(index_TE))=0;
                    end
                else
                    if isempty(Gamma_TP_parfor{i}{np})
                        Gamma_TP_parfor{i}{np} = zeros(1,length(ind_TE));
                    end
                    if isempty(TE_TP_parfor{i}{np})
                        TE_TP_parfor{i}{np} = zeros(1,length(ind_TE));
                    end
                end
            end
        end
    end
end

% Output assignment, assembling back full arrays
for i=1:nTP
    % Prepopulate with zeros
    ME_TP{i} = cell(nParticles,1);
    for p = 1:nParticles
        ME_TP{i}{p} = zeros(nTE, nc);
    end
    TE_TP{i} = cell(nParticles,1);
    for p = 1:nParticles
        TE_TP{i}{p} = zeros(1, nTE); % or whatever your TE array shape is
    end
    Gamma_TP{i} = cell(nParticles,1);
    for p = 1:nParticles
        Gamma_TP{i}{p} = zeros(1, nTE);
    end

    % Now assign only to selected indices
    ME_TP{i}(index) = ME_TP_parfor{i};
    TE_TP{i}(index) = TE_TP_parfor{i};
    Gamma_TP{i}(index) = Gamma_TP_parfor{i};
end

ME_solid_w_TP = cell(1,nTP);
for i=1:nTP
    ME_solid_w_TP{i} = repmat(ME_solid_TP_w(:,i),1,nTE);
end

zero_aux = zeros(size(ME_solid,1),nTE);

for i=1:nTP
    Gamma_TP{i} = cell2mat(Gamma_TP{i});
    Gamma_TP{i}(ME_solid(:,3)>1-2E-8,:) = zero_aux(ME_solid(:,3)>1-2E-8,:); % WARNING move to settings
    Gamma_TP{i}(isnan(Gamma_TP{i}))=0;
    TE_TP{i} = cell2mat(TE_TP{i});
end

ME_solid_w_solid = zeros(size(ME_solid,1),nTE);
for i=1:nTP
    ME_solid_w_solid = ME_solid_w_solid + ME_solid_w_TP{i};
end

ME_solid(:,4:3+nTE) = 0;
for i=1:nTP
    ME_solid(:,4:3+nTE) = ME_solid(:,4:3+nTE) + TE_TP{i}.*ME_solid_w_TP{i}./ME_solid_w_solid;
end

end