function [Gamma_Ol,Gamma_Cpx,Gamma_Opx,Gamma_Gt,Gamma_Sp,Gamma_Pl,ME_Ol,ME_Cpx,ME_Opx,ME_Gt,ME_Sp,ME_Pl,REE_Ol,REE_Cpx,REE_Opx,REE_Gt,REE_Sp,REE_Pl,MRadi,MRadi_0,ME_solid] = ...
    diffusion1D(ME_Ol_eq,ME_Cpx_eq,ME_Opx_eq,ME_Gt_eq,ME_Sp_eq,ME_Pl_eq,ME_melt_eq,BC_Ol_0,BC_Cpx_0,BC_Opx_0,BC_Gt_0,BC_Sp_0,BC_Pl_0,BC_Ol,BC_Cpx,BC_Opx,BC_Gt,BC_Sp,BC_Pl,ME_solid,ME_solid0,ME_solid_TP_v,ME_solid_TP_v0,ME_solid_TP_w,ME_Ol,ME_Cpx,ME_Opx,ME_Gt,ME_Sp,ME_Pl,ME_solid_Rho_TP,ME_solid_Rho_TP0,Mdiff_Ol,Mdiff_Cpx,Mdiff_Opx,Mdiff_Gt,Mdiff_Sp,Mdiff_Pl,timestep,nc,nREE,diff_nTP,diff_nREE,fix_radi,correct_REE_ind,correct_REE_met)
% Compute diffusion profiles inside thermodynamic phases (TP) for different
% oxides/trace elements (b)
% Compute the grain boundary gradient
%
% Oxide/trace elements are computed in %wt of solid (s). Thus:
%   ->           \sum_{b} C_s^{b,TP} \neq 1  (sum of oxides concentrations inside a TP is not = 1)
%   -> \sum_{TP} \sum_{b} C_s^{b,TP}    = 1  (sum of oxides concentrations over all TP is     = 1)
% Not sure about this. Need to check



% Compute boundary conditions for each oxide (b) inside each thermodynamic phase
% BC_Ol  = zeros(size(ME_solid,1),5);
% BC_Cpx = zeros(size(ME_solid,1),5);
% BC_Opx = zeros(size(ME_solid,1),5);
% BC_Gt  = zeros(size(ME_solid,1),5);
% BC_Sp  = zeros(size(ME_solid,1),5);
% 
% for index = 1:5
%     F_Ol  = scatteredInterpolant(XT(:,1),XT(:,2),NOl(:,index),'linear','nearest');
%     F_Cpx = scatteredInterpolant(XT(:,1),XT(:,2),NCpx(:,index),'linear','nearest');
%     F_Opx = scatteredInterpolant(XT(:,1),XT(:,2),NOpx(:,index),'linear','nearest');
%     F_Gt  = scatteredInterpolant(XT(:,1),XT(:,2),NGt(:,index),'linear','nearest');
%     F_Sp  = scatteredInterpolant(XT(:,1),XT(:,2),NSp(:,index),'linear','nearest');
%     BC_Ol(:,index)  = F_Ol([ME_solid(:,1),ME_solid(:,2)]);
%     BC_Cpx(:,index) = F_Cpx([ME_solid(:,1),ME_solid(:,2)]);
%     BC_Opx(:,index) = F_Opx([ME_solid(:,1),ME_solid(:,2)]);
%     BC_Gt(:,index)  = F_Gt([ME_solid(:,1),ME_solid(:,2)]);
%     BC_Sp(:,index)  = F_Sp([ME_solid(:,1),ME_solid(:,2)]);
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
ME_solid_TP_v        = ME_solid_TP_v./repmat(sum(ME_solid_TP_v,2),1,size(ME_solid_TP_v,2));

% Compute Radios
Rradi_Ol  = ((ME_solid_TP_v(:,1))./((4/3)*pi.*ME_solid(:,3+nREE+1))).^(1/3); if diff_nTP(1)==0;  Rradi_Ol  = 0*Rradi_Ol; end;   if fix_radi(1)~=0; Rradi_Ol   = fix_radi(1)*ones(size(Rradi_Ol));  ME_solid(:,3+nREE+1) = ME_solid_TP_v(:,1)./((4/3)*pi.*Rradi_Ol.^3); end;
Rradi_Cpx = ((ME_solid_TP_v(:,2))./((4/3)*pi.*ME_solid(:,3+nREE+2))).^(1/3); if diff_nTP(2)==0;  Rradi_Cpx = 0*Rradi_Cpx; end;  if fix_radi(2)~=0; Rradi_Cpx  = fix_radi(2)*ones(size(Rradi_Cpx)); ME_solid(:,3+nREE+2) = ME_solid_TP_v(:,2)./((4/3)*pi.*Rradi_Cpx.^3); end;
Rradi_Opx = ((ME_solid_TP_v(:,3))./((4/3)*pi.*ME_solid(:,3+nREE+3))).^(1/3); if diff_nTP(3)==0;  Rradi_Opx = 0*Rradi_Opx; end;  if fix_radi(3)~=0; Rradi_Opx  = fix_radi(3)*ones(size(Rradi_Opx)); ME_solid(:,3+nREE+3) = ME_solid_TP_v(:,3)./((4/3)*pi.*Rradi_Opx.^3); end;
Rradi_Gt  = ((ME_solid_TP_v(:,4))./((4/3)*pi.*ME_solid(:,3+nREE+4))).^(1/3); if diff_nTP(4)==0;  Rradi_Gt  = 0*Rradi_Gt; end;   if fix_radi(4)~=0; Rradi_Gt   = fix_radi(4)*ones(size(Rradi_Gt));  ME_solid(:,3+nREE+4) = ME_solid_TP_v(:,4)./((4/3)*pi.*Rradi_Gt.^3); end;
Rradi_Sp  = ((ME_solid_TP_v(:,5))./((4/3)*pi.*ME_solid(:,3+nREE+5))).^(1/3); if diff_nTP(5)==0;  Rradi_Sp  = 0*Rradi_Sp; end;   if fix_radi(5)~=0; Rradi_Sp   = fix_radi(5)*ones(size(Rradi_Sp));  ME_solid(:,3+nREE+5) = ME_solid_TP_v(:,5)./((4/3)*pi.*Rradi_Sp.^3); end;
Rradi_Pl  = ((ME_solid_TP_v(:,6))./((4/3)*pi.*ME_solid(:,3+nREE+6))).^(1/3); if diff_nTP(6)==0;  Rradi_Pl  = 0*Rradi_Pl; end;   if fix_radi(6)~=0; Rradi_Pl   = fix_radi(6)*ones(size(Rradi_Pl));  ME_solid(:,3+nREE+6) = ME_solid_TP_v(:,6)./((4/3)*pi.*Rradi_Pl.^3); end;
% Rradi_Ol = 0.005*ones(size(Rradi_Ol));
% Rradi_Cpx = 0.005*ones(size(Rradi_Cpx));
% Rradi_Opx = 0.005*ones(size(Rradi_Opx));
% Rradi_Gt = 0.005*ones(size(Rradi_Gt));
% Rradi_Sp = 0.005*ones(size(Rradi_Sp));

MRadi = [Rradi_Ol Rradi_Cpx Rradi_Opx Rradi_Gt Rradi_Sp Rradi_Pl];

% Normalize, ME_solid_TP_w 
ME_solid_TP_v0        = ME_solid_TP_v0./repmat(sum(ME_solid_TP_v0,2),1,size(ME_solid_TP_v0,2));
% Compute Initial Radios
Rradi_Ol_0  = ((ME_solid_TP_v0(:,1))./((4/3)*pi.*ME_solid0(:,3+nREE+1))).^(1/3); if diff_nTP(1)==0;  Rradi_Ol_0  = 0*Rradi_Ol_0; end;   if fix_radi(1)~=0; Rradi_Ol_0   = fix_radi(1)*ones(size(Rradi_Ol_0));  ME_solid0(:,3+nREE+1) = ME_solid_TP_v0(:,1)./((4/3)*pi.*Rradi_Ol_0.^3); end;
Rradi_Cpx_0 = ((ME_solid_TP_v0(:,2))./((4/3)*pi.*ME_solid0(:,3+nREE+2))).^(1/3); if diff_nTP(2)==0;  Rradi_Cpx_0 = 0*Rradi_Cpx_0; end;  if fix_radi(2)~=0; Rradi_Cpx_0  = fix_radi(2)*ones(size(Rradi_Cpx_0)); ME_solid0(:,3+nREE+2) = ME_solid_TP_v0(:,2)./((4/3)*pi.*Rradi_Cpx_0.^3); end;
Rradi_Opx_0 = ((ME_solid_TP_v0(:,3))./((4/3)*pi.*ME_solid0(:,3+nREE+3))).^(1/3); if diff_nTP(3)==0;  Rradi_Opx_0 = 0*Rradi_Opx_0; end;  if fix_radi(3)~=0; Rradi_Opx_0  = fix_radi(3)*ones(size(Rradi_Opx_0)); ME_solid0(:,3+nREE+3) = ME_solid_TP_v0(:,3)./((4/3)*pi.*Rradi_Opx_0.^3); end;
Rradi_Gt_0  = ((ME_solid_TP_v0(:,4))./((4/3)*pi.*ME_solid0(:,3+nREE+4))).^(1/3); if diff_nTP(4)==0;  Rradi_Gt_0  = 0*Rradi_Gt_0; end;   if fix_radi(4)~=0; Rradi_Gt_0   = fix_radi(4)*ones(size(Rradi_Gt_0));  ME_solid0(:,3+nREE+4) = ME_solid_TP_v0(:,4)./((4/3)*pi.*Rradi_Gt_0.^3); end;
Rradi_Sp_0  = ((ME_solid_TP_v0(:,5))./((4/3)*pi.*ME_solid0(:,3+nREE+5))).^(1/3); if diff_nTP(5)==0;  Rradi_Sp_0  = 0*Rradi_Sp_0; end;   if fix_radi(5)~=0; Rradi_Sp_0   = fix_radi(5)*ones(size(Rradi_Sp_0));  ME_solid0(:,3+nREE+5) = ME_solid_TP_v0(:,5)./((4/3)*pi.*Rradi_Sp_0.^3); end;
Rradi_Pl_0  = ((ME_solid_TP_v0(:,6))./((4/3)*pi.*ME_solid0(:,3+nREE+6))).^(1/3); if diff_nTP(6)==0;  Rradi_Pl_0  = 0*Rradi_Pl_0; end;   if fix_radi(6)~=0; Rradi_Pl_0   = fix_radi(6)*ones(size(Rradi_Pl_0));  ME_solid0(:,3+nREE+6) = ME_solid_TP_v0(:,6)./((4/3)*pi.*Rradi_Pl_0.^3); end;

MRadi_0 = [Rradi_Ol_0 Rradi_Cpx_0 Rradi_Opx_0 Rradi_Gt_0 Rradi_Sp_0 Rradi_Pl_0];

Gamma_Ol  = mat2cell(zeros(size(ME_solid,1),nREE),ones(size(ME_solid,1),1),nREE); REE_Ol  = Gamma_Ol;
Gamma_Cpx = mat2cell(zeros(size(ME_solid,1),nREE),ones(size(ME_solid,1),1),nREE); REE_Cpx = Gamma_Cpx;
Gamma_Opx = mat2cell(zeros(size(ME_solid,1),nREE),ones(size(ME_solid,1),1),nREE); REE_Opx = Gamma_Opx;
Gamma_Gt  = mat2cell(zeros(size(ME_solid,1),nREE),ones(size(ME_solid,1),1),nREE); REE_Gt  = Gamma_Gt;
Gamma_Sp  = mat2cell(zeros(size(ME_solid,1),nREE),ones(size(ME_solid,1),1),nREE); REE_Sp  = Gamma_Sp;
Gamma_Pl  = mat2cell(zeros(size(ME_solid,1),nREE),ones(size(ME_solid,1),1),nREE); REE_Pl  = Gamma_Pl;

% we approximate the initial profiles with the previous one but with the mesh changed
% loop over all the particles to compute 1D profiles [#TP x #oxides/trace elements]

% transient time stepping. Number of timesteps
nt = 2;
n1 = nt-1;
aux_BC = repmat(0:n1,nREE,1);
ind_REE = find(diff_nREE==1);
index = any(correct_REE_ind,2);
index_met = any(correct_REE_met,2); 
% index_met = index;
index_keep = any(correct_REE_ind==0,2);
ind_Particles = 1:size(ME_solid,1);
vec_particles     = ind_Particles(index);
vec_particles_met = ind_Particles(index_met);
% index = any(correct_REE_ind,2)
% parfor nparticles = 1:size(MR_solid,1)
% parfor nparticles = 1:size(ME_solid,1);
% parfor nparticles_index = 1:length(vec_particles)

BC_Ol_0_parfor = BC_Ol_0(index,:);  BC_Cpx_0_parfor = BC_Cpx_0(index,:);  BC_Opx_0_parfor = BC_Opx_0(index,:);  BC_Gt_0_parfor = BC_Gt_0(index,:);  BC_Sp_0_parfor = BC_Sp_0(index,:);  BC_Pl_0_parfor = BC_Pl_0(index,:);
BC_Ol_parfor   = BC_Ol(index,:);    BC_Cpx_parfor   = BC_Cpx(index,:);    BC_Opx_parfor   = BC_Opx(index,:);    BC_Gt_parfor   = BC_Gt(index,:);    BC_Sp_parfor   = BC_Sp(index,:);    BC_Pl_parfor   = BC_Pl(index,:);    
ME_Ol_parfor   = ME_Ol(index);      ME_Cpx_parfor   = ME_Cpx(index);      ME_Opx_parfor   = ME_Opx(index);      ME_Gt_parfor   = ME_Gt(index);      ME_Sp_parfor   = ME_Sp(index);      ME_Pl_parfor   = ME_Pl(index);        
Gamma_Ol_parfor = Gamma_Ol(index);  Gamma_Cpx_parfor= Gamma_Cpx(index);   Gamma_Opx_parfor= Gamma_Opx(index);   Gamma_Gt_parfor= Gamma_Gt(index);   Gamma_Sp_parfor= Gamma_Sp(index);   Gamma_Pl_parfor= Gamma_Pl(index);        
Rradi_Ol_parfor = Rradi_Ol(index);  Rradi_Cpx_parfor = Rradi_Cpx(index);  Rradi_Opx_parfor = Rradi_Opx(index);  Rradi_Gt_parfor = Rradi_Gt(index);  Rradi_Sp_parfor = Rradi_Sp(index);  Rradi_Pl_parfor = Rradi_Pl(index);  
Rradi_Ol_0_parfor = Rradi_Ol_0(index);  Rradi_Cpx_0_parfor = Rradi_Cpx_0(index);  Rradi_Opx_0_parfor = Rradi_Opx_0(index);  Rradi_Gt_0_parfor = Rradi_Gt_0(index);  Rradi_Sp_0_parfor = Rradi_Sp_0(index);  Rradi_Pl_0_parfor = Rradi_Pl_0(index);  
Mdiff_Ol_parfor = Mdiff_Ol(index,:);    Mdiff_Cpx_parfor = Mdiff_Cpx(index,:);    Mdiff_Opx_parfor = Mdiff_Opx(index,:);    Mdiff_Gt_parfor = Mdiff_Gt(index,:);    Mdiff_Sp_parfor = Mdiff_Sp(index,:);    Mdiff_Pl_parfor = Mdiff_Pl(index,:);    
REE_Ol_parfor = REE_Ol(index);          REE_Cpx_parfor = REE_Cpx(index);          REE_Opx_parfor = REE_Opx(index);          REE_Gt_parfor = REE_Gt(index);          REE_Sp_parfor = REE_Sp(index);          REE_Pl_parfor = REE_Pl(index);          
ME_solid_parfor  = ME_solid(index,:);
ME_solid0_parfor = ME_solid0(index,:);
ME_solid_Rho_TP_parfor  = ME_solid_Rho_TP(index,:);
ME_solid_Rho_TP0_parfor = ME_solid_Rho_TP0(index,:);
ME_solid_TP_v_parfor = ME_solid_TP_v(index,:);
ME_solid_TP_v0_parfor = ME_solid_TP_v0(index,:);


% ME_Ol_parfor_met   = ME_Ol(index_met);      ME_Cpx_parfor_met   = ME_Cpx(index_met);      ME_Opx_parfor_met   = ME_Opx(index_met);      ME_Gt_parfor_met   = ME_Gt(index_met);      ME_Sp_parfor_met   = ME_Sp(index_met);      ME_Pl_parfor_met   = ME_Pl(index_met);        
ME_Ol_parfor_met   = ME_Ol_eq(index_met);   ME_Cpx_parfor_met   = ME_Cpx_eq(index_met);   ME_Opx_parfor_met   = ME_Opx_eq(index_met);   ME_Gt_parfor_met   = ME_Gt_eq(index_met);   ME_Sp_parfor_met   = ME_Sp_eq(index_met);   ME_Pl_parfor_met   = ME_Pl_eq(index_met);        
Gamma_Ol_parfor_met = Gamma_Ol(index_met);  Gamma_Cpx_parfor_met= Gamma_Cpx(index_met);   Gamma_Opx_parfor_met= Gamma_Opx(index_met);   Gamma_Gt_parfor_met= Gamma_Gt(index_met);   Gamma_Sp_parfor_met= Gamma_Sp(index_met);   Gamma_Pl_parfor_met= Gamma_Pl(index_met);        
Rradi_Ol_parfor_met = Rradi_Ol(index_met);  Rradi_Cpx_parfor_met = Rradi_Cpx(index_met);  Rradi_Opx_parfor_met = Rradi_Opx(index_met);  Rradi_Gt_parfor_met = Rradi_Gt(index_met);  Rradi_Sp_parfor_met = Rradi_Sp(index_met);  Rradi_Pl_parfor_met = Rradi_Pl(index_met);  
REE_Ol_parfor_met = REE_Ol(index_met);          REE_Cpx_parfor_met = REE_Cpx(index_met);          REE_Opx_parfor_met = REE_Opx(index_met);          REE_Gt_parfor_met = REE_Gt(index_met);          REE_Sp_parfor_met = REE_Sp(index_met);          REE_Pl_parfor_met = REE_Pl(index_met);          

exp_factor = 0.3;
% exp_factor = -1;
for index_mesh=1:nc; aux_mesh(index_mesh)=(1-(1/index_mesh)^exp_factor)/(1-(1/nc)^exp_factor); end

tol_REE = 1e-15;

parfor nparticles = 1:length(vec_particles)
% for nparticles = 1:length(vec_particles)
    time_discrete = [0:timestep/nt:timestep];
    if diff_nTP(1)==1; BC_Ol_aux  = repmat(BC_Ol_0_parfor(nparticles,:)',1,nt)  + aux_BC.*repmat((BC_Ol_parfor(nparticles,:)'  - BC_Ol_0_parfor(nparticles,:)')/n1,1,nt); end
    if diff_nTP(2)==1; BC_Cpx_aux = repmat(BC_Cpx_0_parfor(nparticles,:)',1,nt) + aux_BC.*repmat((BC_Cpx_parfor(nparticles,:)' - BC_Cpx_0_parfor(nparticles,:)')/n1,1,nt); end
    if diff_nTP(3)==1; BC_Opx_aux = repmat(BC_Opx_0_parfor(nparticles,:)',1,nt) + aux_BC.*repmat((BC_Opx_parfor(nparticles,:)' - BC_Opx_0_parfor(nparticles,:)')/n1,1,nt); end
    if diff_nTP(4)==1; BC_Gt_aux  = repmat(BC_Gt_0_parfor(nparticles,:)',1,nt)  + aux_BC.*repmat((BC_Gt_parfor(nparticles,:)'  - BC_Gt_0_parfor(nparticles,:)')/n1,1,nt); end
    if diff_nTP(5)==1; BC_Sp_aux  = repmat(BC_Sp_0_parfor(nparticles,:)',1,nt)  + aux_BC.*repmat((BC_Sp_parfor(nparticles,:)'  - BC_Sp_0_parfor(nparticles,:)')/n1,1,nt); end
    if diff_nTP(6)==1; BC_Pl_aux  = repmat(BC_Pl_0_parfor(nparticles,:)',1,nt)  + aux_BC.*repmat((BC_Pl_parfor(nparticles,:)'  - BC_Pl_0_parfor(nparticles,:)')/n1,1,nt); end

    
    % Solve
    for index_time = 1:(nt-1)
        % we need 3 steps for the solver
        tspan = [time_discrete(nt) 0.5*(time_discrete(nt+1)+time_discrete(nt)) time_discrete(nt+1)];
        % update initial condition for transient problem
        MRini_Ol_0  = ME_Ol_parfor{nparticles};
        MRini_Cpx_0 = ME_Cpx_parfor{nparticles};
        MRini_Opx_0 = ME_Opx_parfor{nparticles};
        MRini_Gt_0  = ME_Gt_parfor{nparticles};
        MRini_Sp_0  = ME_Sp_parfor{nparticles};
        MRini_Pl_0  = ME_Pl_parfor{nparticles};
        
        % update initial condition for PC problem
        MRini_Ol_PC  = ME_Ol_parfor{nparticles};
        MRini_Cpx_PC = ME_Cpx_parfor{nparticles};
        MRini_Opx_PC = ME_Opx_parfor{nparticles};
        MRini_Gt_PC  = ME_Gt_parfor{nparticles};
        MRini_Sp_PC  = ME_Sp_parfor{nparticles};
        MRini_Pl_PC  = ME_Pl_parfor{nparticles};
        
        % Matrix Assembly
        % Espatial discretization
        x1_Ol    = Rradi_Ol_parfor(nparticles);       xmesh_Ol    = aux_mesh*x1_Ol;      if x1_Ol==0;    xmesh_Ol=zeros(1,nc);    end
        x1_Ol_0  = Rradi_Ol_0_parfor(nparticles);     xmesh_Ol_0  = aux_mesh*x1_Ol_0;    if x1_Ol_0==0;  if x1_Ol==0; xmesh_Ol_0=zeros(1,nc);   else; xmesh_Ol_0  = xmesh_Ol;  end;  end
        x1_Cpx    = Rradi_Cpx_parfor(nparticles);     xmesh_Cpx    = aux_mesh*x1_Cpx;    if x1_Cpx==0;   xmesh_Cpx=zeros(1,nc);   end
        x1_Cpx_0  = Rradi_Cpx_0_parfor(nparticles);   xmesh_Cpx_0  = aux_mesh*x1_Cpx_0;  if x1_Cpx_0==0; if x1_Cpx==0; xmesh_Cpx_0=zeros(1,nc); else; xmesh_Cpx_0 = xmesh_Cpx; end;  end
        x1_Opx    = Rradi_Opx_parfor(nparticles);     xmesh_Opx    = aux_mesh*x1_Opx;    if x1_Opx==0;   xmesh_Opx=zeros(1,nc);   end
        x1_Opx_0  = Rradi_Opx_0_parfor(nparticles);   xmesh_Opx_0  = aux_mesh*x1_Opx_0;  if x1_Opx_0==0; if x1_Opx==0; xmesh_Opx_0=zeros(1,nc); else; xmesh_Opx_0 = xmesh_Opx; end;  end
        x1_Gt    = Rradi_Gt_parfor(nparticles);       xmesh_Gt    = aux_mesh*x1_Gt;      if x1_Gt==0;    xmesh_Gt=zeros(1,nc);    end
        x1_Gt_0  = Rradi_Gt_0_parfor(nparticles);     xmesh_Gt_0  = aux_mesh*x1_Gt_0;    if x1_Gt_0==0;  if x1_Gt==0; xmesh_Gt_0=zeros(1,nc);   else; xmesh_Gt_0  = xmesh_Gt;  end;  end
        x1_Sp    = Rradi_Sp_parfor(nparticles);       xmesh_Sp    = aux_mesh*x1_Sp;      if x1_Sp==0;    xmesh_Sp=zeros(1,nc);    end
        x1_Sp_0  = Rradi_Sp_0_parfor(nparticles);     xmesh_Sp_0  = aux_mesh*x1_Sp_0;    if x1_Sp_0==0;  if x1_Sp==0; xmesh_Sp_0=zeros(1,nc);   else; xmesh_Sp_0  = xmesh_Sp;  end;  end
        x1_Pl    = Rradi_Pl_parfor(nparticles);       xmesh_Pl    = aux_mesh*x1_Pl;      if x1_Pl==0;    xmesh_Pl=zeros(1,nc);    end
        x1_Pl_0  = Rradi_Pl_0_parfor(nparticles);     xmesh_Pl_0  = aux_mesh*x1_Pl_0;    if x1_Pl_0==0;  if x1_Pl==0; xmesh_Pl_0=zeros(1,nc);   else; xmesh_Pl_0  = xmesh_Pl;  end;  end
       
        xipg = [-1/sqrt(3) 1/sqrt(3)]';
        wpg = [1 1]';
        
        % Shape functions and its derivatives on the reference element
        N_mef   =  [(1-xipg)/2 (1+xipg)/2];
        Nxi_mef =  [-1/2 1/2; -1/2 1/2];

        % Matrices obtained by discretizng a convection-diffusion equation
        if diff_nTP(1)==1; [M_Ol,K_Ol,~]   = matrices_1D(xmesh_Ol,xipg,wpg,N_mef,Nxi_mef);   M_Ol   = spdiags(sum(M_Ol,2), 0,length(M_Ol),length(M_Ol));    end
        if diff_nTP(2)==1; [M_Cpx,K_Cpx,~] = matrices_1D(xmesh_Cpx,xipg,wpg,N_mef,Nxi_mef);  M_Cpx  = spdiags(sum(M_Cpx,2), 0,length(M_Cpx),length(M_Cpx)); end
        if diff_nTP(3)==1; [M_Opx,K_Opx,~] = matrices_1D(xmesh_Opx,xipg,wpg,N_mef,Nxi_mef);  M_Opx  = spdiags(sum(M_Opx,2), 0,length(M_Opx),length(M_Opx)); end
        if diff_nTP(4)==1; [M_Gt,K_Gt,~]   = matrices_1D(xmesh_Gt,xipg,wpg,N_mef,Nxi_mef);   M_Gt   = spdiags(sum(M_Gt,2), 0,length(M_Gt),length(M_Gt));    end
        if diff_nTP(5)==1; [M_Sp,K_Sp,~]   = matrices_1D(xmesh_Sp,xipg,wpg,N_mef,Nxi_mef);   M_Sp   = spdiags(sum(M_Sp,2), 0,length(M_Sp),length(M_Sp));    end
        if diff_nTP(6)==1; [M_Pl,K_Pl,~]   = matrices_1D(xmesh_Pl,xipg,wpg,N_mef,Nxi_mef);   M_Pl   = spdiags(sum(M_Pl,2), 0,length(M_Pl),length(M_Pl));    end
        
        for index_REE = 1:length(ind_REE)
%                         if Rradi_Ol(vec_particlesnparticles)~=0;  ME_Ol{vec_particlesnparticles}(ind_REE(index),:)  = diffusion_1D_PDE(Rradi_Ol(vec_particlesnparticles),MRini_Ol_0(ind_REE(index),:),Mdiff_Ol(vec_particlesnparticles,ind_REE(index)),BC_Ol_aux(ind_REE(index),index_time+1),tspan,nc); end
%                         if Rradi_Cpx(vec_particlesnparticles)~=0; ME_Cpx{vec_particlesnparticles}(ind_REE(index),:) = diffusion_1D_PDE(Rradi_Cpx(vec_particlesnparticles),MRini_Cpx_0(ind_REE(index),:),Mdiff_Cpx(vec_particlesnparticles,ind_REE(index)),BC_Cpx_aux(ind_REE(index),index_time+1),tspan,nc); end
%                         if Rradi_Opx(vec_particlesnparticles)~=0; ME_Opx{vec_particlesnparticles}(ind_REE(index),:) = diffusion_1D_PDE(Rradi_Opx(vec_particlesnparticles),MRini_Opx_0(ind_REE(index),:),Mdiff_Opx(vec_particlesnparticles,ind_REE(index)),BC_Opx_aux(ind_REE(index),index_time+1),tspan,nc); end
%                         if Rradi_Gt(vec_particlesnparticles)~=0;  ME_Gt{vec_particlesnparticles}(ind_REE(index),:)  = diffusion_1D_PDE(Rradi_Gt(vec_particlesnparticles),MRini_Gt_0(ind_REE(index),:),Mdiff_Gt(vec_particlesnparticles,ind_REE(index)),BC_Gt_aux(ind_REE(index),index_time+1),tspan,nc); end
%                         if Rradi_Sp(vec_particlesnparticles)~=0;  ME_Sp{vec_particlesnparticles}(ind_REE(index),:)  = diffusion_1D_PDE(Rradi_Sp(vec_particlesnparticles),MRini_Sp_0(ind_REE(index),:),Mdiff_Sp(vec_particlesnparticles,ind_REE(index)),BC_Sp_aux(ind_REE(index),index_time+1),tspan,nc); end
%                         if Rradi_Pl(vec_particlesnparticles)~=0;  ME_Pl{vec_particlesnparticles}(ind_REE(index),:)  = diffusion_1D_PDE(Rradi_Pl(vec_particlesnparticles),MRini_Pl_0(ind_REE(index),:),Mdiff_Pl(vec_particlesnparticles,ind_REE(index)),BC_Pl_aux(ind_REE(index),index_time+1),tspan,nc); end
            
            %             if Rradi_Ol(vec_particlesnparticles)~=0;  ME_Ol{vec_particlesnparticles}(index,:)  = diffusion_1D_PDE(Rradi_Ol(vec_particlesnparticles),MRini_Ol_0(index,:),Mdiff_Ol(vec_particlesnparticles,index),BC_Ol_aux(index,index_time+1),tspan,nc); end
            %             if Rradi_Cpx(vec_particlesnparticles)~=0; ME_Cpx{vec_particlesnparticles}(index,:) = diffusion_1D_PDE(Rradi_Cpx(vec_particlesnparticles),MRini_Cpx_0(index,:),Mdiff_Cpx(vec_particlesnparticles,index),BC_Cpx_aux(index,index_time+1),tspan,nc); end
            %             if Rradi_Opx(vec_particlesnparticles)~=0; ME_Opx{vec_particlesnparticles}(index,:) = diffusion_1D_PDE(Rradi_Opx(vec_particlesnparticles),MRini_Opx_0(index,:),Mdiff_Opx(vec_particlesnparticles,index),BC_Opx_aux(index,index_time+1),tspan,nc); end
            %             if Rradi_Gt(vec_particlesnparticles)~=0;  ME_Gt{vec_particlesnparticles}(index,:)  = diffusion_1D_PDE(Rradi_Gt(vec_particlesnparticles),MRini_Gt_0(index,:),Mdiff_Gt(vec_particlesnparticles,index),BC_Gt_aux(index,index_time+1),tspan,nc); end
            %             if Rradi_Sp(vec_particlesnparticles)~=0;  ME_Sp{vec_particlesnparticles}(index,:)  = diffusion_1D_PDE(Rradi_Sp(vec_particlesnparticles),MRini_Sp_0(index,:),Mdiff_Sp(vec_particlesnparticles,index),BC_Sp_aux(index,index_time+1),tspan,nc); end
            %             if Rradi_Pl(vec_particlesnparticles)~=0;  ME_Pl{vec_particlesnparticles}(index,:)  = diffusion_1D_PDE(Rradi_Pl(vec_particlesnparticles),MRini_Pl_0(index,:),Mdiff_Pl(vec_particlesnparticles,index),BC_Pl_aux(index,index_time+1),tspan,nc); end
            
%     W = 1/2;
%     w = 1; 
%     % Boundary conditions (lagrange multipliers method)
% Accd = zeros(1,nc);
% Accd(1,nc)=1; 
% Accd = sparse(Accd);
% bccd = [BC_Ol_aux(ind_REE(index),index_time+1)];
% M_Ol  = spdiags(sum(M_Ol,2), 0,length(M_Ol),length(M_Ol));
%             Sol = Galerkin(W,w,0,Mdiff_Ol(vec_particles_parfor(nparticles),ind_REE(index)),0*xmesh_Ol',K_Ol,M_Ol,0,xmesh_Ol,timestep,1,MRini_Ol_0(ind_REE(index),:)',Accd,bccd*5);
% keyboard
% if vec_particles(nparticles) == 1188; keyboard; end
            % Phase change Part
            if x1_Ol > x1_Ol_0;   MRini_Ol_PC(ind_REE(index_REE),:)  = interp1(xmesh_Ol_0,MRini_Ol_0(ind_REE(index_REE),:),xmesh_Ol,'linear');    MRini_Ol_PC(ind_REE(index_REE),xmesh_Ol>=x1_Ol_0)    = BC_Ol_aux(ind_REE(index_REE),index_time+1);   elseif  x1_Ol < x1_Ol_0;   MRini_Ol_PC(ind_REE(index_REE),:)  = interp1(xmesh_Ol_0,MRini_Ol_0(ind_REE(index_REE),:),xmesh_Ol,'linear'); end
            if x1_Cpx > x1_Cpx_0; MRini_Cpx_PC(ind_REE(index_REE),:) = interp1(xmesh_Cpx_0,MRini_Cpx_0(ind_REE(index_REE),:),xmesh_Cpx,'linear'); MRini_Cpx_PC(ind_REE(index_REE),xmesh_Cpx>=x1_Cpx_0) = BC_Cpx_aux(ind_REE(index_REE),index_time+1);  elseif  x1_Cpx < x1_Cpx_0; MRini_Cpx_PC(ind_REE(index_REE),:) = interp1(xmesh_Cpx_0,MRini_Cpx_0(ind_REE(index_REE),:),xmesh_Cpx,'linear'); end
            if x1_Opx > x1_Opx_0; MRini_Opx_PC(ind_REE(index_REE),:) = interp1(xmesh_Opx_0,MRini_Opx_0(ind_REE(index_REE),:),xmesh_Opx,'linear'); MRini_Opx_PC(ind_REE(index_REE),xmesh_Opx>=x1_Opx_0) = BC_Opx_aux(ind_REE(index_REE),index_time+1);  elseif  x1_Opx < x1_Opx_0; MRini_Opx_PC(ind_REE(index_REE),:) = interp1(xmesh_Opx_0,MRini_Opx_0(ind_REE(index_REE),:),xmesh_Opx,'linear'); end
            if x1_Gt > x1_Gt_0;   MRini_Gt_PC(ind_REE(index_REE),:)  = interp1(xmesh_Gt_0,MRini_Gt_0(ind_REE(index_REE),:),xmesh_Gt,'linear');    MRini_Gt_PC(ind_REE(index_REE),xmesh_Gt>=x1_Gt_0)    = BC_Gt_aux(ind_REE(index_REE),index_time+1);   elseif  x1_Gt < x1_Gt_0;   MRini_Gt_PC(ind_REE(index_REE),:)  = interp1(xmesh_Gt_0,MRini_Gt_0(ind_REE(index_REE),:),xmesh_Gt,'linear'); end
            if x1_Sp > x1_Sp_0;   MRini_Sp_PC(ind_REE(index_REE),:)  = interp1(xmesh_Sp_0,MRini_Sp_0(ind_REE(index_REE),:),xmesh_Sp,'linear');    MRini_Sp_PC(ind_REE(index_REE),xmesh_Sp>=x1_Sp_0)    = BC_Sp_aux(ind_REE(index_REE),index_time+1);   elseif  x1_Sp < x1_Sp_0;   MRini_Sp_PC(ind_REE(index_REE),:)  = interp1(xmesh_Sp_0,MRini_Sp_0(ind_REE(index_REE),:),xmesh_Sp,'linear'); end
            if x1_Pl > x1_Pl_0;   MRini_Pl_PC(ind_REE(index_REE),:)  = interp1(xmesh_Pl_0,MRini_Pl_0(ind_REE(index_REE),:),xmesh_Pl,'linear');    MRini_Pl_PC(ind_REE(index_REE),xmesh_Pl>=x1_Pl_0)    = BC_Pl_aux(ind_REE(index_REE),index_time+1);   elseif  x1_Pl < x1_Pl_0;   MRini_Pl_PC(ind_REE(index_REE),:)  = interp1(xmesh_Pl_0,MRini_Pl_0(ind_REE(index_REE),:),xmesh_Pl,'linear'); end
            
            % Diffusion Part
%             if Rradi_Ol_parfor(nparticles)~=0;  ME_Ol_parfor{nparticles}(ind_REE(index_REE),:)  = solveWithLagrangeMultipliers(M_Ol + timestep*0.5*K_Ol*Mdiff_Ol_parfor(nparticles,ind_REE(index_REE)) , M_Ol*MRini_Ol_PC(ind_REE(index_REE),:)' + 0.5*timestep*K_Ol*Mdiff_Ol_parfor(nparticles,ind_REE(index_REE))*MRini_Ol_PC(ind_REE(index_REE),:)', [nc BC_Ol_aux(ind_REE(index_REE),index_time+1)]); end;           if all(abs(ME_Ol_parfor{nparticles}(ind_REE(index_REE),1:end-1)- MRini_Ol_0(ind_REE(index_REE),1:end-1))./ ME_Ol_parfor{nparticles}(ind_REE(index_REE),1:end-1)<tol_REE);     ME_Ol_parfor{nparticles}(ind_REE(index_REE),:)  = MRini_Ol_0(ind_REE(index_REE),:);   end;
            if Rradi_Ol_parfor(nparticles)~=0;  ME_Ol_parfor{nparticles}(ind_REE(index_REE),:)  = solveWithLagrangeMultipliers(M_Ol + timestep*K_Ol*Mdiff_Ol_parfor(nparticles,ind_REE(index_REE)) , M_Ol*MRini_Ol_PC(ind_REE(index_REE),:)', [nc BC_Ol_aux(ind_REE(index_REE),index_time+1)]); end;           if all(abs(ME_Ol_parfor{nparticles}(ind_REE(index_REE),1:end-1)- MRini_Ol_0(ind_REE(index_REE),1:end-1))./ ME_Ol_parfor{nparticles}(ind_REE(index_REE),1:end-1)<tol_REE);     ME_Ol_parfor{nparticles}(ind_REE(index_REE),:)  = MRini_Ol_0(ind_REE(index_REE),:);   end;
%             if Rradi_Cpx_parfor(nparticles)~=0; ME_Cpx_parfor{nparticles}(ind_REE(index_REE),:) = solveWithLagrangeMultipliers(M_Cpx + timestep*0.5*K_Cpx*Mdiff_Cpx_parfor(nparticles,ind_REE(index_REE)) , M_Cpx*MRini_Cpx_PC(ind_REE(index_REE),:)' + 0.5*timestep*K_Cpx*Mdiff_Cpx_parfor(nparticles,ind_REE(index_REE))*MRini_Cpx_PC(ind_REE(index_REE),:)', [nc BC_Cpx_aux(ind_REE(index_REE),index_time+1)]); end;  if all(abs(ME_Cpx_parfor{nparticles}(ind_REE(index_REE),1:end-1)- MRini_Cpx_0(ind_REE(index_REE),1:end-1))./ ME_Cpx_parfor{nparticles}(ind_REE(index_REE),1:end-1)<tol_REE);  ME_Cpx_parfor{nparticles}(ind_REE(index_REE),:) = MRini_Cpx_0(ind_REE(index_REE),:);  end;
            if Rradi_Cpx_parfor(nparticles)~=0; ME_Cpx_parfor{nparticles}(ind_REE(index_REE),:) = solveWithLagrangeMultipliers(M_Cpx + timestep*K_Cpx*Mdiff_Cpx_parfor(nparticles,ind_REE(index_REE)) , M_Cpx*MRini_Cpx_PC(ind_REE(index_REE),:)', [nc BC_Cpx_aux(ind_REE(index_REE),index_time+1)]); end;  if all(abs(ME_Cpx_parfor{nparticles}(ind_REE(index_REE),1:end-1)- MRini_Cpx_0(ind_REE(index_REE),1:end-1))./ ME_Cpx_parfor{nparticles}(ind_REE(index_REE),1:end-1)<tol_REE);  ME_Cpx_parfor{nparticles}(ind_REE(index_REE),:) = MRini_Cpx_0(ind_REE(index_REE),:);  end;
%             if Rradi_Opx_parfor(nparticles)~=0; ME_Opx_parfor{nparticles}(ind_REE(index_REE),:) = solveWithLagrangeMultipliers(M_Opx + timestep*0.5*K_Opx*Mdiff_Opx_parfor(nparticles,ind_REE(index_REE)) , M_Opx*MRini_Opx_PC(ind_REE(index_REE),:)' + 0.5*timestep*K_Opx*Mdiff_Opx_parfor(nparticles,ind_REE(index_REE))*MRini_Opx_PC(ind_REE(index_REE),:)', [nc BC_Opx_aux(ind_REE(index_REE),index_time+1)]); end;  if all(abs(ME_Opx_parfor{nparticles}(ind_REE(index_REE),1:end-1)- MRini_Opx_0(ind_REE(index_REE),1:end-1))./ ME_Opx_parfor{nparticles}(ind_REE(index_REE),1:end-1)<tol_REE);  ME_Opx_parfor{nparticles}(ind_REE(index_REE),:) = MRini_Opx_0(ind_REE(index_REE),:);  end;
            if Rradi_Opx_parfor(nparticles)~=0; ME_Opx_parfor{nparticles}(ind_REE(index_REE),:) = solveWithLagrangeMultipliers(M_Opx + timestep*K_Opx*Mdiff_Opx_parfor(nparticles,ind_REE(index_REE)) , M_Opx*MRini_Opx_PC(ind_REE(index_REE),:)', [nc BC_Opx_aux(ind_REE(index_REE),index_time+1)]); end;  if all(abs(ME_Opx_parfor{nparticles}(ind_REE(index_REE),1:end-1)- MRini_Opx_0(ind_REE(index_REE),1:end-1))./ ME_Opx_parfor{nparticles}(ind_REE(index_REE),1:end-1)<tol_REE);  ME_Opx_parfor{nparticles}(ind_REE(index_REE),:) = MRini_Opx_0(ind_REE(index_REE),:);  end;
%             if Rradi_Gt_parfor(nparticles)~=0;  ME_Gt_parfor{nparticles}(ind_REE(index_REE),:)  = solveWithLagrangeMultipliers(M_Gt + timestep*0.5*K_Gt*Mdiff_Gt_parfor(nparticles,ind_REE(index_REE)) , M_Gt*MRini_Gt_PC(ind_REE(index_REE),:)' + 0.5*timestep*K_Gt*Mdiff_Gt_parfor(nparticles,ind_REE(index_REE))*MRini_Gt_PC(ind_REE(index_REE),:)', [nc BC_Gt_aux(ind_REE(index_REE),index_time+1)]); end;           if all(abs(ME_Gt_parfor{nparticles}(ind_REE(index_REE),1:end-1)- MRini_Gt_0(ind_REE(index_REE),1:end-1))./ ME_Gt_parfor{nparticles}(ind_REE(index_REE),1:end-1)<tol_REE);     ME_Gt_parfor{nparticles}(ind_REE(index_REE),:)  = MRini_Gt_0(ind_REE(index_REE),:);   end;
            if Rradi_Gt_parfor(nparticles)~=0;  ME_Gt_parfor{nparticles}(ind_REE(index_REE),:)  = solveWithLagrangeMultipliers(M_Gt + timestep*K_Gt*Mdiff_Gt_parfor(nparticles,ind_REE(index_REE)) , M_Gt*MRini_Gt_PC(ind_REE(index_REE),:)', [nc BC_Gt_aux(ind_REE(index_REE),index_time+1)]); end;           if all(abs(ME_Gt_parfor{nparticles}(ind_REE(index_REE),1:end-1)- MRini_Gt_0(ind_REE(index_REE),1:end-1))./ ME_Gt_parfor{nparticles}(ind_REE(index_REE),1:end-1)<tol_REE);     ME_Gt_parfor{nparticles}(ind_REE(index_REE),:)  = MRini_Gt_0(ind_REE(index_REE),:);   end;
%             if Rradi_Sp_parfor(nparticles)~=0;  ME_Sp_parfor{nparticles}(ind_REE(index_REE),:)  = solveWithLagrangeMultipliers(M_Sp + timestep*0.5*K_Sp*Mdiff_Sp_parfor(nparticles,ind_REE(index_REE)) , M_Sp*MRini_Sp_PC(ind_REE(index_REE),:)' + 0.5*timestep*K_Sp*Mdiff_Sp_parfor(nparticles,ind_REE(index_REE))*MRini_Sp_PC(ind_REE(index_REE),:)', [nc BC_Sp_aux(ind_REE(index_REE),index_time+1)]); end;           if all(abs(ME_Sp_parfor{nparticles}(ind_REE(index_REE),1:end-1)- MRini_Sp_0(ind_REE(index_REE),1:end-1))./ ME_Sp_parfor{nparticles}(ind_REE(index_REE),1:end-1)<tol_REE);     ME_Sp_parfor{nparticles}(ind_REE(index_REE),:)  = MRini_Sp_0(ind_REE(index_REE),:);   end;
            if Rradi_Sp_parfor(nparticles)~=0;  ME_Sp_parfor{nparticles}(ind_REE(index_REE),:)  = solveWithLagrangeMultipliers(M_Sp + timestep*K_Sp*Mdiff_Sp_parfor(nparticles,ind_REE(index_REE)) , M_Sp*MRini_Sp_PC(ind_REE(index_REE),:)', [nc BC_Sp_aux(ind_REE(index_REE),index_time+1)]); end;           if all(abs(ME_Sp_parfor{nparticles}(ind_REE(index_REE),1:end-1)- MRini_Sp_0(ind_REE(index_REE),1:end-1))./ ME_Sp_parfor{nparticles}(ind_REE(index_REE),1:end-1)<tol_REE);     ME_Sp_parfor{nparticles}(ind_REE(index_REE),:)  = MRini_Sp_0(ind_REE(index_REE),:);   end;
%             if Rradi_Pl_parfor(nparticles)~=0;  ME_Pl_parfor{nparticles}(ind_REE(index_REE),:)  = solveWithLagrangeMultipliers(M_Pl + timestep*0.5*K_Pl*Mdiff_Pl_parfor(nparticles,ind_REE(index_REE)) , M_Pl*MRini_Pl_PC(ind_REE(index_REE),:)' + 0.5*timestep*K_Pl*Mdiff_Pl_parfor(nparticles,ind_REE(index_REE))*MRini_Pl_PC(ind_REE(index_REE),:)', [nc BC_Pl_aux(ind_REE(index_REE),index_time+1)]); end;           if all(abs(ME_Pl_parfor{nparticles}(ind_REE(index_REE),1:end-1)- MRini_Pl_0(ind_REE(index_REE),1:end-1))./ ME_Pl_parfor{nparticles}(ind_REE(index_REE),1:end-1)<tol_REE);     ME_Pl_parfor{nparticles}(ind_REE(index_REE),:)  = MRini_Pl_0(ind_REE(index_REE),:);   end;
            if Rradi_Pl_parfor(nparticles)~=0;  ME_Pl_parfor{nparticles}(ind_REE(index_REE),:)  = solveWithLagrangeMultipliers(M_Pl + timestep*K_Pl*Mdiff_Pl_parfor(nparticles,ind_REE(index_REE)) , M_Pl*MRini_Pl_PC(ind_REE(index_REE),:)', [nc BC_Pl_aux(ind_REE(index_REE),index_time+1)]); end;           if all(abs(ME_Pl_parfor{nparticles}(ind_REE(index_REE),1:end-1)- MRini_Pl_0(ind_REE(index_REE),1:end-1))./ ME_Pl_parfor{nparticles}(ind_REE(index_REE),1:end-1)<tol_REE);     ME_Pl_parfor{nparticles}(ind_REE(index_REE),:)  = MRini_Pl_0(ind_REE(index_REE),:);   end;

            % Mass balance
            integral0 = trapz(xmesh_Ol_0,4*pi*xmesh_Ol_0.^2.*MRini_Ol_0(ind_REE(index_REE),:));      integral1 = trapz(xmesh_Ol,4*pi*xmesh_Ol.^2.*ME_Ol_parfor{nparticles}(ind_REE(index_REE),:));      % Gamma_Ol_parfor{nparticles}(ind_REE(index_REE))   = (ME_solid_parfor(nparticles,3+nREE+1).*ME_solid_Rho_TP_parfor(nparticles,1).*integral1-ME_solid0_parfor(nparticles,3+nREE+1).*ME_solid_Rho_TP0_parfor(nparticles,1).*integral0)/timestep;
            Gamma_Ol_parfor{nparticles}(ind_REE(index_REE))   =  (integral1-integral0).*ME_solid0_parfor(nparticles,3+nREE+1).*ME_solid_Rho_TP0_parfor(nparticles,1)/timestep;
            Gamma_Ol_parfor{nparticles}(ind_REE(index_REE))   =  (integral1.*ME_solid_parfor(nparticles,3+nREE+1).*ME_solid_Rho_TP_parfor(nparticles,1)-integral0.*ME_solid0_parfor(nparticles,3+nREE+1).*ME_solid_Rho_TP0_parfor(nparticles,1))/timestep;
%             Gamma_Ol_parfor{nparticles}(ind_REE(index_REE))   =  ((integral1/trapz(xmesh_Ol,4*pi*xmesh_Ol.^2))*ME_solid_TP_v_parfor(nparticles,1)*ME_solid_Rho_TP_parfor(nparticles,1) ...
%                                                                -  (integral0/trapz(xmesh_Ol_0,4*pi*xmesh_Ol_0.^2))*ME_solid_TP_v0_parfor(nparticles,1)*ME_solid_Rho_TP0_parfor(nparticles,1))/timestep;
            %             REE_Ol_parfor{nparticles}(index)  = ME_solid_parfor(nparticles,12).*integral1;   if xmesh_Ol(end)==0;  REE_Ol_parfor{nparticles}(index)=0;  end
%             REE_Ol_parfor{nparticles}(ind_REE(index))  = integral1/(4/3*pi*xmesh_Ol(end).^3)/(1-ME_solid_TP_v_parfor(nparticles,7));   if xmesh_Ol(end)==0;  REE_Ol_parfor{nparticles}(ind_REE(index))=0;  end
            REE_Ol_parfor{nparticles}(ind_REE(index_REE)) = integral1/trapz(xmesh_Ol,4*pi*xmesh_Ol.^2);
            if xmesh_Ol(end)==0;  REE_Ol_parfor{nparticles}(ind_REE(index_REE))=0;  end
%             if Mdiff_Ol_parfor(nparticles,ind_REE(index))==0;   if integral1~=0; ME_Ol_parfor{nparticles}(ind_REE(index),:)  = integral0/integral1*ME_Ol_parfor{nparticles}(ind_REE(index),:);  end; end
            
            integral0 = trapz(xmesh_Cpx_0,4*pi*xmesh_Cpx_0.^2.*MRini_Cpx_0(ind_REE(index_REE),:));   integral1 = trapz(xmesh_Cpx,4*pi*xmesh_Cpx.^2.*ME_Cpx_parfor{nparticles}(ind_REE(index_REE),:));   % Gamma_Cpx_parfor{nparticles}(ind_REE(index_REE))  = (ME_solid_parfor(nparticles,3+nREE+2).*ME_solid_Rho_TP_parfor(nparticles,2).*integral1-ME_solid0_parfor(nparticles,3+nREE+2).*ME_solid_Rho_TP0_parfor(nparticles,2).*integral0)/timestep;
            Gamma_Cpx_parfor{nparticles}(ind_REE(index_REE))  =  (integral1-integral0).*ME_solid0_parfor(nparticles,3+nREE+2).*ME_solid_Rho_TP0_parfor(nparticles,2)/timestep;
            Gamma_Cpx_parfor{nparticles}(ind_REE(index_REE))  =  (integral1.*ME_solid_parfor(nparticles,3+nREE+2).*ME_solid_Rho_TP_parfor(nparticles,2)-integral0.*ME_solid0_parfor(nparticles,3+nREE+2).*ME_solid_Rho_TP0_parfor(nparticles,2))/timestep;
            %             REE_Cpx_parfor{nparticles}(index) = ME_solid_parfor(nparticles,13).*integral1;  if xmesh_Cpx(end)==0; REE_Cpx_parfor{nparticles}(index)=0; end
%             REE_Cpx_parfor{nparticles}(ind_REE(index)) = integral1/(4/3*pi*xmesh_Cpx(end).^3)/(1-ME_solid_TP_v_parfor(nparticles,7));  if xmesh_Cpx(end)==0; REE_Cpx_parfor{nparticles}(ind_REE(index))=0; end
            REE_Cpx_parfor{nparticles}(ind_REE(index_REE)) = integral1/ trapz(xmesh_Cpx,4*pi*xmesh_Cpx.^2);
            if xmesh_Cpx(end)==0; REE_Cpx_parfor{nparticles}(ind_REE(index_REE))=0; end
%             if Mdiff_Cpx_parfor(nparticles,ind_REE(index))==0;  if integral1~=0; ME_Cpx_parfor{nparticles}(ind_REE(index),:) = integral0/integral1*ME_Cpx_parfor{nparticles}(ind_REE(index),:);  end; end
            
            integral0 = trapz(xmesh_Opx_0,4*pi*xmesh_Opx_0.^2.*MRini_Opx_0(ind_REE(index_REE),:));   integral1 = trapz(xmesh_Opx,4*pi*xmesh_Opx.^2.*ME_Opx_parfor{nparticles}(ind_REE(index_REE),:));   % Gamma_Opx_parfor{nparticles}(ind_REE(index_REE))  = (ME_solid_parfor(nparticles,3+nREE+3).*ME_solid_Rho_TP_parfor(nparticles,3).*integral1-ME_solid0_parfor(nparticles,3+nREE+3).*ME_solid_Rho_TP0_parfor(nparticles,3).*integral0)/timestep;
            Gamma_Opx_parfor{nparticles}(ind_REE(index_REE))  =  (integral1-integral0).*ME_solid0_parfor(nparticles,3+nREE+3).*ME_solid_Rho_TP0_parfor(nparticles,3)/timestep;
            Gamma_Opx_parfor{nparticles}(ind_REE(index_REE))  =  (integral1.*ME_solid_parfor(nparticles,3+nREE+3).*ME_solid_Rho_TP_parfor(nparticles,3)-integral0.*ME_solid0_parfor(nparticles,3+nREE+3).*ME_solid_Rho_TP0_parfor(nparticles,3))/timestep;
            %             REE_Opx_parfor{nparticles}(index) = ME_solid_parfor(nparticles,14).*integral1;  if xmesh_Opx(end)==0; REE_Opx_parfor{nparticles}(index)=0; end
%             REE_Opx_parfor{nparticles}(ind_REE(index)) = integral1/(4/3*pi*xmesh_Opx(end).^3)/(1-ME_solid_TP_v_parfor(nparticles,7));  if xmesh_Opx(end)==0; REE_Opx_parfor{nparticles}(ind_REE(index))=0; end
            REE_Opx_parfor{nparticles}(ind_REE(index_REE)) = integral1/trapz(xmesh_Opx,4*pi*xmesh_Opx.^2);
            if xmesh_Opx(end)==0; REE_Opx_parfor{nparticles}(ind_REE(index_REE))=0; end
%             if Mdiff_Opx_parfor(nparticles,ind_REE(index))==0;  if integral1~=0; ME_Opx_parfor{nparticles}(ind_REE(index),:) = integral0/integral1*ME_Opx_parfor{nparticles}(ind_REE(index),:);  end; end
            
            integral0 = trapz(xmesh_Gt_0,4*pi*xmesh_Gt_0.^2.*MRini_Gt_0(ind_REE(index_REE),:));      integral1 = trapz(xmesh_Gt,4*pi*xmesh_Gt.^2.*ME_Gt_parfor{nparticles}(ind_REE(index_REE),:));      % Gamma_Gt_parfor{nparticles}(ind_REE(index_REE))   = (ME_solid_parfor(nparticles,3+nREE+4).*ME_solid_Rho_TP_parfor(nparticles,4).*integral1-ME_solid0_parfor(nparticles,3+nREE+4).*ME_solid_Rho_TP0_parfor(nparticles,4).*integral0)/timestep;
            Gamma_Gt_parfor{nparticles}(ind_REE(index_REE))   = (integral1-integral0).*ME_solid0_parfor(nparticles,3+nREE+4).*ME_solid_Rho_TP0_parfor(nparticles,4)/timestep;
            Gamma_Gt_parfor{nparticles}(ind_REE(index_REE))   = (integral1.*ME_solid_parfor(nparticles,3+nREE+4).*ME_solid_Rho_TP_parfor(nparticles,4)-integral0.*ME_solid0_parfor(nparticles,3+nREE+4).*ME_solid_Rho_TP0_parfor(nparticles,4))/timestep;
            %             REE_Gt_parfor{nparticles}(index)  = ME_solid_parfor(nparticles,15).*integral1;   if xmesh_Gt(end)==0;  REE_Gt_parfor{nparticles}(index)=0;  end
%             REE_Gt_parfor{nparticles}(ind_REE(index))  = integral1/(4/3*pi*xmesh_Gt(end).^3)/(1-ME_solid_TP_v_parfor(nparticles,7));   if xmesh_Gt(end)==0;  REE_Gt_parfor{nparticles}(ind_REE(index))=0;  end
            REE_Gt_parfor{nparticles}(ind_REE(index_REE))  = integral1/trapz(xmesh_Gt,4*pi*xmesh_Gt.^2);
            if xmesh_Gt(end)==0;  REE_Gt_parfor{nparticles}(ind_REE(index_REE))=0;  end
%             if Mdiff_Gt_parfor(nparticles,ind_REE(index))==0;   if integral1~=0; ME_Gt_parfor{nparticles}(ind_REE(index),:)  = integral0/integral1*ME_Gt_parfor{nparticles}(ind_REE(index),:);  end; end
            
            integral0 = trapz(xmesh_Sp_0,4*pi*xmesh_Sp_0.^2.*MRini_Sp_0(ind_REE(index_REE),:));      integral1 = trapz(xmesh_Sp,4*pi*xmesh_Sp.^2.*ME_Sp_parfor{nparticles}(ind_REE(index_REE),:));      % Gamma_Sp_parfor{nparticles}(ind_REE(index_REE))   = (ME_solid_parfor(nparticles,3+nREE+5).*ME_solid_Rho_TP_parfor(nparticles,5).*integral1-ME_solid0_parfor(nparticles,3+nREE+5).*ME_solid_Rho_TP0_parfor(nparticles,5).*integral0)/timestep;
            Gamma_Sp_parfor{nparticles}(ind_REE(index_REE))   =  (integral1-integral0).*ME_solid0_parfor(nparticles,3+nREE+5).*ME_solid_Rho_TP0_parfor(nparticles,5)/timestep;
            Gamma_Sp_parfor{nparticles}(ind_REE(index_REE))   =  (integral1.*ME_solid_parfor(nparticles,3+nREE+5).*ME_solid_Rho_TP_parfor(nparticles,5)-integral0.*ME_solid0_parfor(nparticles,3+nREE+5).*ME_solid_Rho_TP0_parfor(nparticles,5))/timestep;
            %             REE_Sp_parfor{nparticles}(index)  = ME_solid_parfor(nparticles,16).*integral1;   if xmesh_Sp(end)==0;  REE_Sp_parfor{nparticles}(index)=0;  end
%             REE_Sp_parfor{nparticles}(ind_REE(index))  = integral1/(4/3*pi*xmesh_Sp(end).^3)/(1-ME_solid_TP_v_parfor(nparticles,7));   if xmesh_Sp(end)==0;  REE_Sp_parfor{nparticles}(ind_REE(index))=0;  end
            REE_Sp_parfor{nparticles}(ind_REE(index_REE))  = integral1/trapz(xmesh_Sp,4*pi*xmesh_Sp.^2);
            if xmesh_Sp(end)==0;  REE_Sp_parfor{nparticles}(ind_REE(index_REE))=0;  end
%             if Mdiff_Sp_parfor(nparticles,ind_REE(index))==0;   if integral1~=0; ME_Sp_parfor{nparticles}(ind_REE(index),:)  = integral0/integral1*ME_Sp_parfor{nparticles}(ind_REE(index),:);  end; end
            
            integral0 = trapz(xmesh_Pl_0,4*pi*xmesh_Pl_0.^2.*MRini_Pl_0(ind_REE(index_REE),:));      integral1 = trapz(xmesh_Pl,4*pi*xmesh_Pl.^2.*ME_Pl_parfor{nparticles}(ind_REE(index_REE),:));      % Gamma_Pl_parfor{nparticles}(ind_REE(index_REE))   = (ME_solid_parfor(nparticles,3+nREE+6).*ME_solid_Rho_TP_parfor(nparticles,6).*integral1-ME_solid0_parfor(nparticles,3+nREE+6).*ME_solid_Rho_TP0_parfor(nparticles,6).*integral0)/timestep;
            Gamma_Pl_parfor{nparticles}(ind_REE(index_REE))   =  (integral1-integral0).*ME_solid0_parfor(nparticles,3+nREE+6).*ME_solid_Rho_TP0_parfor(nparticles,6)/timestep;
            Gamma_Pl_parfor{nparticles}(ind_REE(index_REE))   =  (integral1.*ME_solid_parfor(nparticles,3+nREE+6).*ME_solid_Rho_TP_parfor(nparticles,6)-integral0.*ME_solid0_parfor(nparticles,3+nREE+6).*ME_solid_Rho_TP0_parfor(nparticles,6))/timestep;
            %             REE_Pl_parfor{nparticles}(index)  = ME_solid_parfor(nparticles,17).*integral1;   if xmesh_Pl(end)==0;  REE_Pl_parfor{nparticles}(index)=0;  end
%             REE_Pl_parfor{nparticles}(ind_REE(index))  = integral1/(4/3*pi*xmesh_Pl(end).^3)/(1-ME_solid_TP_v_parfor(nparticles,7));   if xmesh_Pl(end)==0;  REE_Pl_parfor{nparticles}(ind_REE(index))=0;  end
            REE_Pl_parfor{nparticles}(ind_REE(index_REE))  = integral1/trapz(xmesh_Pl,4*pi*xmesh_Pl.^2);
            if xmesh_Pl(end)==0;  REE_Pl_parfor{nparticles}(ind_REE(index_REE))=0;  end
%             if Mdiff_Pl_parfor(nparticles,ind_REE(index))==0;   if integral1~=0; ME_Pl_parfor{nparticles}(ind_REE(index),:)  = integral0/integral1*ME_Pl_parfor{nparticles}(ind_REE(index),:);  end; end
            
            %             % Without n_tp
            %             integral0 = trapz(xmesh_Ol_0,4*pi*xmesh_Ol_0.^2.*MRini_Ol_0(index,:));      integral1 = trapz(xmesh_Ol,4*pi*xmesh_Ol.^2.*ME_Ol_parfor{nparticles}(index,:));      Gamma_Ol_parfor{nparticles}(index)   = (ME_solid_Rho_TP_parfor(nparticles,1).*integral1/(4/3*pi*xmesh_Ol(end).^3)-ME_solid_Rho_TP0_parfor(nparticles,1).*integral0/(4/3*pi*xmesh_Ol_0(end).^3))/timestep;
            %             REE_Ol_parfor{nparticles}(index)  = integral1/(4/3*pi*xmesh_Ol(end).^3);   if xmesh_Ol(end)==0;  REE_Ol_parfor{nparticles}(index)=0;  end
            %             if Mdiff_Ol_parfor(nparticles,index)==0;   ME_Ol_parfor{nparticles}(index,:)  = integral0/integral1*ME_Ol_parfor{nparticles}(index,:);  end
            %
            %             integral0 = trapz(xmesh_Cpx_0,4*pi*xmesh_Cpx_0.^2.*MRini_Cpx_0(index,:));   integral1 = trapz(xmesh_Cpx,4*pi*xmesh_Cpx.^2.*ME_Cpx_parfor{nparticles}(index,:));   Gamma_Cpx_parfor{nparticles}(index)  = (ME_solid_Rho_TP_parfor(nparticles,2).*integral1/(4/3*pi*xmesh_Cpx(end).^3)-ME_solid_Rho_TP0_parfor(nparticles,2).*integral0/(4/3*pi*xmesh_Cpx_0(end).^3))/timestep;
            %             REE_Cpx_parfor{nparticles}(index) = integral1/(4/3*pi*xmesh_Cpx(end).^3);  if xmesh_Cpx(end)==0; REE_Cpx_parfor{nparticles}(index)=0; end
            %             if Mdiff_Cpx_parfor(nparticles,index)==0;  ME_Cpx_parfor{nparticles}(index,:) = integral0/integral1*ME_Cpx_parfor{nparticles}(index,:); end
            %
            %             integral0 = trapz(xmesh_Opx_0,4*pi*xmesh_Opx_0.^2.*MRini_Opx_0(index,:));   integral1 = trapz(xmesh_Opx,4*pi*xmesh_Opx.^2.*ME_Opx_parfor{nparticles}(index,:));   Gamma_Opx_parfor{nparticles}(index)  = (ME_solid_Rho_TP_parfor(nparticles,3).*integral1/(4/3*pi*xmesh_Opx(end).^3)-ME_solid_Rho_TP0_parfor(nparticles,3).*integral0/(4/3*pi*xmesh_Opx_0(end).^3))/timestep;
            %             REE_Opx_parfor{nparticles}(index) = integral1/(4/3*pi*xmesh_Opx(end).^3);  if xmesh_Opx(end)==0; REE_Opx_parfor{nparticles}(index)=0; end
            %             if Mdiff_Opx_parfor(nparticles,index)==0;  ME_Opx_parfor{nparticles}(index,:) = integral0/integral1*ME_Opx_parfor{nparticles}(index,:); end
            %
            %             integral0 = trapz(xmesh_Gt_0,4*pi*xmesh_Gt_0.^2.*MRini_Gt_0(index,:));      integral1 = trapz(xmesh_Gt,4*pi*xmesh_Gt.^2.*ME_Gt_parfor{nparticles}(index,:));      Gamma_Gt_parfor{nparticles}(index)   = (ME_solid_Rho_TP_parfor(nparticles,4).*integral1/(4/3*pi*xmesh_Gt(end).^3)-ME_solid_Rho_TP0_parfor(nparticles,4).*integral0/(4/3*pi*xmesh_Gt_0(end).^3))/timestep;
            %             REE_Gt_parfor{nparticles}(index)  = integral1/(4/3*pi*xmesh_Gt(end).^3);   if xmesh_Gt(end)==0;  REE_Gt_parfor{nparticles}(index)=0;  end
            %             if Mdiff_Gt_parfor(nparticles,index)==0;   ME_Gt_parfor{nparticles}(index,:)  = integral0/integral1*ME_Gt_parfor{nparticles}(index,:);  end
            %
            %             integral0 = trapz(xmesh_Sp_0,4*pi*xmesh_Sp_0.^2.*MRini_Sp_0(index,:));      integral1 = trapz(xmesh_Sp,4*pi*xmesh_Sp.^2.*ME_Sp_parfor{nparticles}(index,:));      Gamma_Sp_parfor{nparticles}(index)   = (ME_solid_Rho_TP_parfor(nparticles,5).*integral1/(4/3*pi*xmesh_Sp(end).^3)-ME_solid_Rho_TP0_parfor(nparticles,5).*integral0/(4/3*pi*xmesh_Sp_0(end).^3))/timestep;
            %             REE_Sp_parfor{nparticles}(index)  = integral1/(4/3*pi*xmesh_Sp(end).^3);   if xmesh_Sp(end)==0;  REE_Sp_parfor{nparticles}(index)=0;  end
            %             if Mdiff_Sp_parfor(nparticles,index)==0;   ME_Sp_parfor{nparticles}(index,:)  = integral0/integral1*ME_Sp_parfor{nparticles}(index,:);  end
            %
            %             integral0 = trapz(xmesh_Pl_0,4*pi*xmesh_Pl_0.^2.*MRini_Pl_0(index,:));      integral1 = trapz(xmesh_Pl,4*pi*xmesh_Pl.^2.*ME_Pl_parfor{nparticles}(index,:));      Gamma_Pl_parfor{nparticles}(index)   = (ME_solid_Rho_TP_parfor(nparticles,6).*integral1/(4/3*pi*xmesh_Pl(end).^3)-ME_solid_Rho_TP0_parfor(nparticles,6).*integral0/(4/3*pi*xmesh_Pl_0(end).^3))/timestep;
            %             REE_Pl_parfor{nparticles}(index)  = integral1/(4/3*pi*xmesh_Pl(end).^3);   if xmesh_Pl(end)==0;  REE_Pl_parfor{nparticles}(index)=0;  end
            %             if Mdiff_Pl_parfor(nparticles,index)==0;   ME_Pl_parfor{nparticles}(index,:)  = integral0/integral1*ME_Pl_parfor{nparticles}(index,:);  end
        end
        
        %         if ME_solid_parfor(nparticles,3)>1-1E-4
        %             Gamma_Ol_parfor{nparticles}=zeros(1,8);
        %             Gamma_Cpx_parfor{nparticles}=zeros(1,8);
        %             Gamma_Opx_parfor{nparticles}=zeros(1,8);
        %             Gamma_Gt_parfor{nparticles}=zeros(1,8);
        %             Gamma_Sp_parfor{nparticles}=zeros(1,8);
        %             Gamma_Pl_parfor{nparticles}=zeros(1,8);
        %         end
    end
    
    
    Grad_oxides_Ol_parfor(nparticles,:) = (ME_Ol_parfor{nparticles}(:,end-1)-ME_Ol_parfor{nparticles}(:,end-1))/(Rradi_Ol_parfor(nparticles)/(nt-1));
    Grad_oxides_Cpx_parfor(nparticles,:) = (ME_Cpx_parfor{nparticles}(:,end-1)-ME_Cpx_parfor{nparticles}(:,end-1))/(Rradi_Cpx_parfor(nparticles)/(nt-1));
    Grad_oxides_Opx_parfor(nparticles,:) = (ME_Opx_parfor{nparticles}(:,end-1)-ME_Opx_parfor{nparticles}(:,end-1))/(Rradi_Opx_parfor(nparticles)/(nt-1));
    Grad_oxides_Gt_parfor(nparticles,:) = (ME_Gt_parfor{nparticles}(:,end-1)-ME_Gt_parfor{nparticles}(:,end-1))/(Rradi_Gt_parfor(nparticles)/(nt-1));
    Grad_oxides_Sp_parfor(nparticles,:) = (ME_Sp_parfor{nparticles}(:,end-1)-ME_Sp_parfor{nparticles}(:,end-1))/(Rradi_Sp_parfor(nparticles)/(nt-1));
    Grad_oxides_Pl_parfor(nparticles,:) = (ME_Pl_parfor{nparticles}(:,end-1)-ME_Pl_parfor{nparticles}(:,end-1))/(Rradi_Pl_parfor(nparticles)/(nt-1));
    
end

parfor nparticles = 1:length(vec_particles_met)
    
    % Espatial discretization
        x1_Ol    = Rradi_Ol_parfor_met(nparticles);       xmesh_Ol    = aux_mesh*x1_Ol;      if x1_Ol==0;    xmesh_Ol=zeros(1,nc);    end
        x1_Cpx    = Rradi_Cpx_parfor_met(nparticles);     xmesh_Cpx    = aux_mesh*x1_Cpx;    if x1_Cpx==0;   xmesh_Cpx=zeros(1,nc);   end
        x1_Opx    = Rradi_Opx_parfor_met(nparticles);     xmesh_Opx    = aux_mesh*x1_Opx;    if x1_Opx==0;   xmesh_Opx=zeros(1,nc);   end
        x1_Gt    = Rradi_Gt_parfor_met(nparticles);       xmesh_Gt    = aux_mesh*x1_Gt;      if x1_Gt==0;    xmesh_Gt=zeros(1,nc);    end
        x1_Sp    = Rradi_Sp_parfor_met(nparticles);       xmesh_Sp    = aux_mesh*x1_Sp;      if x1_Sp==0;    xmesh_Sp=zeros(1,nc);    end
        x1_Pl    = Rradi_Pl_parfor_met(nparticles);       xmesh_Pl    = aux_mesh*x1_Pl;      if x1_Pl==0;    xmesh_Pl=zeros(1,nc);    end       
    
    for index_REE = 1:length(ind_REE)
        
        % Mass balance
        integral1 = trapz(xmesh_Ol,4*pi*xmesh_Ol.^2.*ME_Ol_parfor_met{nparticles}(ind_REE(index_REE),:));      Gamma_Ol_parfor_met{nparticles}(ind_REE(index_REE))   = 0;
        REE_Ol_parfor_met{nparticles}(ind_REE(index_REE)) = integral1/trapz(xmesh_Ol,4*pi*xmesh_Ol.^2);
        if xmesh_Ol(end)==0;  REE_Ol_parfor_met{nparticles}(ind_REE(index_REE))=0;  ME_Ol_parfor_met{nparticles}(ind_REE(index_REE),:) = 0*ME_Ol_parfor_met{nparticles}(ind_REE(index_REE),:); end
        
        integral1 = trapz(xmesh_Cpx,4*pi*xmesh_Cpx.^2.*ME_Cpx_parfor_met{nparticles}(ind_REE(index_REE),:));   Gamma_Cpx_parfor_met{nparticles}(ind_REE(index_REE))  = 0;
        REE_Cpx_parfor_met{nparticles}(ind_REE(index_REE)) = integral1/ trapz(xmesh_Cpx,4*pi*xmesh_Cpx.^2);
        if xmesh_Cpx(end)==0; REE_Cpx_parfor_met{nparticles}(ind_REE(index_REE))=0;  ME_Cpx_parfor_met{nparticles}(ind_REE(index_REE),:) = 0*ME_Cpx_parfor_met{nparticles}(ind_REE(index_REE),:); end
        
        integral1 = trapz(xmesh_Opx,4*pi*xmesh_Opx.^2.*ME_Opx_parfor_met{nparticles}(ind_REE(index_REE),:));   Gamma_Opx_parfor_met{nparticles}(ind_REE(index_REE))  = 0;
        REE_Opx_parfor_met{nparticles}(ind_REE(index_REE)) = integral1/trapz(xmesh_Opx,4*pi*xmesh_Opx.^2);
        if xmesh_Opx(end)==0; REE_Opx_parfor_met{nparticles}(ind_REE(index_REE))=0;  ME_Opx_parfor_met{nparticles}(ind_REE(index_REE),:) = 0*ME_Opx_parfor_met{nparticles}(ind_REE(index_REE),:); end
        
        integral1 = trapz(xmesh_Gt,4*pi*xmesh_Gt.^2.*ME_Gt_parfor_met{nparticles}(ind_REE(index_REE),:));      Gamma_Gt_parfor_met{nparticles}(ind_REE(index_REE))   = 0;
        REE_Gt_parfor_met{nparticles}(ind_REE(index_REE))  = integral1/trapz(xmesh_Gt,4*pi*xmesh_Gt.^2);
        if xmesh_Gt(end)==0;  REE_Gt_parfor_met{nparticles}(ind_REE(index_REE))=0;  ME_Gt_parfor_met{nparticles}(ind_REE(index_REE),:) = 0*ME_Gt_parfor_met{nparticles}(ind_REE(index_REE),:);  end
        
        integral1 = trapz(xmesh_Sp,4*pi*xmesh_Sp.^2.*ME_Sp_parfor_met{nparticles}(ind_REE(index_REE),:));      Gamma_Sp_parfor_met{nparticles}(ind_REE(index_REE))   = 0;
        REE_Sp_parfor_met{nparticles}(ind_REE(index_REE))  = integral1/trapz(xmesh_Sp,4*pi*xmesh_Sp.^2);
        if xmesh_Sp(end)==0;  REE_Sp_parfor_met{nparticles}(ind_REE(index_REE))=0;  ME_Sp_parfor_met{nparticles}(ind_REE(index_REE),:) = 0*ME_Sp_parfor_met{nparticles}(ind_REE(index_REE),:);  end
        
        integral1 = trapz(xmesh_Pl,4*pi*xmesh_Pl.^2.*ME_Pl_parfor_met{nparticles}(ind_REE(index_REE),:));      Gamma_Pl_parfor_met{nparticles}(ind_REE(index_REE))   = 0;
        REE_Pl_parfor_met{nparticles}(ind_REE(index_REE))  = integral1/trapz(xmesh_Pl,4*pi*xmesh_Pl.^2);
        if xmesh_Pl(end)==0;  REE_Pl_parfor_met{nparticles}(ind_REE(index_REE))=0;  ME_Pl_parfor_met{nparticles}(ind_REE(index_REE),:) = 0*ME_Pl_parfor_met{nparticles}(ind_REE(index_REE),:);  end
    end
    
    
end

ME_Ol(index)        = ME_Ol_parfor;         ME_Cpx(index)        = ME_Cpx_parfor;         ME_Opx(index)        = ME_Opx_parfor;         ME_Gt(index)        = ME_Gt_parfor;         ME_Sp(index)        = ME_Sp_parfor;         ME_Pl(index)        = ME_Pl_parfor;        
ME_Ol(index_met)    = ME_Ol_parfor_met;     ME_Cpx(index_met)    = ME_Cpx_parfor_met;     ME_Opx(index_met)    = ME_Opx_parfor_met;     ME_Gt(index_met)    = ME_Gt_parfor_met;     ME_Sp(index_met)    = ME_Sp_parfor_met;     ME_Pl(index_met)    = ME_Pl_parfor_met;        
Rradi_Ol(index)     = Rradi_Ol_parfor;      Rradi_Cpx(index)     = Rradi_Cpx_parfor;      Rradi_Opx(index)     = Rradi_Opx_parfor;      Rradi_Gt(index)     = Rradi_Gt_parfor;      Rradi_Sp(index)     = Rradi_Sp_parfor;      Rradi_Pl(index)     = Rradi_Pl_parfor;  
Rradi_Ol(index_met) = Rradi_Ol_parfor_met;  Rradi_Cpx(index_met) = Rradi_Cpx_parfor_met;  Rradi_Opx(index_met) = Rradi_Opx_parfor_met;  Rradi_Gt(index_met) = Rradi_Gt_parfor_met;  Rradi_Sp(index_met) = Rradi_Sp_parfor_met;  Rradi_Pl(index_met) = Rradi_Pl_parfor_met;  
Mdiff_Ol(index,:)   = Mdiff_Ol_parfor;      Mdiff_Cpx(index,:)   = Mdiff_Cpx_parfor;      Mdiff_Opx(index,:)   = Mdiff_Opx_parfor;      Mdiff_Gt(index,:)   = Mdiff_Gt_parfor;      Mdiff_Sp(index,:)   = Mdiff_Sp_parfor;      Mdiff_Pl(index,:)   = Mdiff_Pl_parfor;    
REE_Ol(index)       = REE_Ol_parfor;        REE_Cpx(index)       = REE_Cpx_parfor;        REE_Opx(index)       = REE_Opx_parfor;        REE_Gt(index)       = REE_Gt_parfor;        REE_Sp(index)       = REE_Sp_parfor;        REE_Pl(index)       = REE_Pl_parfor;          
REE_Ol(index_met)   = REE_Ol_parfor_met;    REE_Cpx(index_met)   = REE_Cpx_parfor_met;    REE_Opx(index_met)   = REE_Opx_parfor_met;    REE_Gt(index_met)   = REE_Gt_parfor_met;    REE_Sp(index_met)   = REE_Sp_parfor_met;    REE_Pl(index_met)   = REE_Pl_parfor_met;          
Gamma_Ol(index)     = Gamma_Ol_parfor;      Gamma_Cpx(index)     = Gamma_Cpx_parfor;      Gamma_Opx(index)     = Gamma_Opx_parfor;      Gamma_Gt(index)     = Gamma_Gt_parfor;      Gamma_Sp(index)     = Gamma_Sp_parfor;      Gamma_Pl(index)     = Gamma_Pl_parfor;        
Gamma_Ol(index_met) = Gamma_Ol_parfor_met;  Gamma_Cpx(index_met) = Gamma_Cpx_parfor_met;  Gamma_Opx(index_met) = Gamma_Opx_parfor_met;  Gamma_Gt(index_met) = Gamma_Gt_parfor_met;  Gamma_Sp(index_met) = Gamma_Sp_parfor_met;  Gamma_Pl(index_met) = Gamma_Pl_parfor_met;        

ME_solid_w_Ol   = repmat(ME_solid_TP_w(:,1),1,nREE);
ME_solid_w_Cpx  = repmat(ME_solid_TP_w(:,2),1,nREE);
ME_solid_w_Opx  = repmat(ME_solid_TP_w(:,3),1,nREE);
ME_solid_w_Gt   = repmat(ME_solid_TP_w(:,4),1,nREE);
ME_solid_w_Sp   = repmat(ME_solid_TP_w(:,5),1,nREE);
ME_solid_w_Pl   = repmat(ME_solid_TP_w(:,6),1,nREE);
ME_solid_w_Melt = repmat(ME_solid_TP_w(:,7),1,nREE);

zero_aux = zeros(size(ME_solid,1),nREE);

Gamma_Ol  = cell2mat(Gamma_Ol);  Gamma_Ol(ME_solid(:,3)>1-2E-8,:)  = zero_aux(ME_solid(:,3)>1-2E-8,:); Gamma_Ol(isnan(Gamma_Ol))=0;
Gamma_Cpx = cell2mat(Gamma_Cpx); Gamma_Cpx(ME_solid(:,3)>1-2E-8,:) = zero_aux(ME_solid(:,3)>1-2E-8,:); Gamma_Cpx(isnan(Gamma_Cpx))=0;
Gamma_Opx = cell2mat(Gamma_Opx); Gamma_Opx(ME_solid(:,3)>1-2E-8,:) = zero_aux(ME_solid(:,3)>1-2E-8,:); Gamma_Opx(isnan(Gamma_Opx))=0;
Gamma_Gt  = cell2mat(Gamma_Gt);  Gamma_Gt(ME_solid(:,3)>1-2E-8,:)  = zero_aux(ME_solid(:,3)>1-2E-8,:); Gamma_Gt(isnan(Gamma_Gt))=0;
Gamma_Sp  = cell2mat(Gamma_Sp);  Gamma_Sp(ME_solid(:,3)>1-2E-8,:)  = zero_aux(ME_solid(:,3)>1-2E-8,:); Gamma_Sp(isnan(Gamma_Sp))=0;
Gamma_Pl  = cell2mat(Gamma_Pl);  Gamma_Pl(ME_solid(:,3)>1-2E-8,:)  = zero_aux(ME_solid(:,3)>1-2E-8,:); Gamma_Pl(isnan(Gamma_Pl))=0;

% REE mass / TP mass
REE_Ol  = cell2mat(REE_Ol);
REE_Cpx = cell2mat(REE_Cpx);
REE_Opx = cell2mat(REE_Opx);
REE_Gt  = cell2mat(REE_Gt);
REE_Sp  = cell2mat(REE_Sp);
REE_Pl  = cell2mat(REE_Pl); 

% Factor TP mass / solid mass
ME_solid_w_solid = ME_solid_w_Ol + ME_solid_w_Cpx + ME_solid_w_Opx + ME_solid_w_Gt + ME_solid_w_Sp + ME_solid_w_Pl;

% REE mass / solid mass
ME_solid(:,4:3+nREE) = REE_Ol.*ME_solid_w_Ol./ME_solid_w_solid ...
                     + REE_Cpx.*ME_solid_w_Cpx./ME_solid_w_solid ...
                     + REE_Opx.*ME_solid_w_Opx./ME_solid_w_solid ...
                     + REE_Gt.*ME_solid_w_Gt./ME_solid_w_solid ...
                     + REE_Sp.*ME_solid_w_Sp./ME_solid_w_solid ...
                     + REE_Pl.*ME_solid_w_Pl./ME_solid_w_solid;     % REE mass / solid mass 




