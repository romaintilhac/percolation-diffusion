function [ME,ME_fluid_solid,ME_solid_TP_w,ME_solid_TP_v,ME_solid_Rho_TP,ME_Ol,ME_Cpx,ME_Opx,ME_Gt,ME_Sp,ME_Pl,REE_Ol,REE_Cpx,REE_Opx,REE_Gt,REE_Sp,REE_Pl,MRadi,input_aux] = ...
    insertReactiveREE_1D_solid(ME,ME_fluid_solid,ME_solid_TP_w,ME_solid_TP_v,ME_solid_Rho_TP,MRadi,ME_Ol,ME_Cpx,ME_Opx,ME_Gt,ME_Sp,ME_Pl,REE_Ol,REE_Cpx,REE_Opx,REE_Gt,REE_Sp,REE_Pl,input_TP_v,input_TP_w,gridx,gridy)

input_aux = 0;
% Remove particle
aux_correct = ME(:,1)>gridx(end) | ME(:,1)<gridx(1);
ME(aux_correct,:) = [];
ME_fluid_solid(aux_correct,:) = [];
ME_solid_TP_w(aux_correct,:) = [];
ME_solid_TP_v(aux_correct,:) = [];
ME_solid_Rho_TP(aux_correct,:) = [];
MRadi(aux_correct,:) = [];
ME_Ol(aux_correct) = [];
ME_Cpx(aux_correct) = [];
ME_Opx(aux_correct) = [];
ME_Gt(aux_correct) = [];
ME_Sp(aux_correct) = [];
ME_Pl(aux_correct) = [];
REE_Ol(aux_correct,:) = [];
REE_Cpx(aux_correct,:) = [];
REE_Opx(aux_correct,:) = [];
REE_Gt(aux_correct,:) = [];
REE_Sp(aux_correct,:) = [];
REE_Pl(aux_correct,:) = [];

aux_correct = ME(:,2)>gridy(end) | ME(:,2)<gridy(1);
ME(aux_correct,:) = [];
ME_fluid_solid(aux_correct,:) = [];
ME_solid_TP_w(aux_correct,:) = [];
ME_solid_TP_v(aux_correct,:) = [];
ME_solid_Rho_TP(aux_correct,:) = [];
MRadi(aux_correct,:) = [];
ME_Ol(aux_correct) = [];
ME_Cpx(aux_correct) = [];
ME_Opx(aux_correct) = [];
ME_Gt(aux_correct) = [];
ME_Sp(aux_correct) = [];
ME_Pl(aux_correct) = [];
REE_Ol(aux_correct,:) = [];
REE_Cpx(aux_correct,:) = [];
REE_Opx(aux_correct,:) = [];
REE_Gt(aux_correct,:) = [];
REE_Sp(aux_correct,:) = [];
REE_Pl(aux_correct,:) = [];

ystp = gridy(2)-gridy(1);

% Updating temperature for markers
MS_aux = ME;
X_MS = MS_aux(:,1);
Y_MS = MS_aux(:,2);
nxt = length(gridx)-1;
nyt = length(gridy)-1;
xsize = gridx(end)-gridx(1);
ysize = gridy(end)-gridy(1);
ystp = gridy(2)-gridy(1);
xstp = gridx(2)-gridx(1);

% % xn = zeros(size(X_MS));
% % xn(X_MS<gridx(nxt+1))=double(int16(X_MS(X_MS<gridx(nxt+1))./xstp-0.5))+1;
% yn=double(int16(Y_MS./ystp-0.5))+1;
% yn(Y_MS<gridy(nyt+1))=double(int16(Y_MS(Y_MS<gridy(nyt+1))./ystp-0.5))+1;
yn = max(cumsum(Y_MS./gridy>1,2),[],2);

% if (xn<1)
%     xn  =   1;
% end
% if (xn>(nxt))
%     xn  =   (nxt);
% end
if (yn<1)
    yn  =   1;
end
if (yn>(nyt))
    yn  =   (nyt);
end


% remove_index = [];
% MS_add = [];
N_particles = accumarray(yn,ones(size(yn)));
correct_elem = (N_particles<8);

if sum(correct_elem)>0; input_aux = 1; end;

mxnum   =   1;         % total number of markers in horizontal direction
mynum_elem = 1;
mynum   =   (length(gridy)-1)*mynum_elem;         % total number of markers in vertical direction
mystep  =   ysize/mynum;        % step between markers in vertical direction
mystep  =   gridy(2:end)-gridy(1:end-1);        % step between markers in vertical direction
stringx = 0;
MX = repmat(stringx,mynum,1);
MX_rand = (rand(size(MX))-0.5)*0;
% stringy = mystep/5:mystep:mynum*mystep-mystep/2;
stringy = gridy(1:end-1)+(gridy(2:end)-gridy(1:end-1))/5;
MY = repmat(stringy',1,mxnum);
MY_rand = (rand(size(MY))-0.5).*mystep'/5;
MS_aux = [MX(:)+0*MX_rand(:) MY(:)+0*MY_rand(:)];


% find closest particle to copy in the desired value of stringy
insert_coord = stringy(correct_elem);
insert_coord_mat = repmat(insert_coord,length(Y_MS),1);
Y_MS_mat = repmat(Y_MS,1,length(insert_coord));
diff = insert_coord_mat-Y_MS_mat;
[~, ind_add] = min(abs(diff));

ME = [ME;...
      ME(ind_add,1) insert_coord' ME(ind_add,3:end)];
ME_fluid_solid = [ME_fluid_solid;...
                      ME_fluid_solid(ind_add,1) insert_coord' ME_fluid_solid(ind_add,3:end)];
ME_solid_TP_w = [ME_solid_TP_w;
                     ME_solid_TP_w(ind_add,:)];
ME_solid_TP_v = [ME_solid_TP_v;...
                     ME_solid_TP_v(ind_add,:)];
ME_solid_Rho_TP = [ME_solid_Rho_TP;...
                       ME_solid_Rho_TP(ind_add,:)];
MRadi = [MRadi;...
         MRadi(ind_add,:)];
ME_Ol = [ME_Ol;...
         ME_Ol(ind_add,:)];
ME_Cpx = [ME_Cpx;...
          ME_Cpx(ind_add,:)];
ME_Opx = [ME_Opx;...
          ME_Opx(ind_add,:)];
ME_Gt = [ME_Gt;...
         ME_Gt(ind_add,:)];
ME_Sp = [ME_Sp;...
         ME_Sp(ind_add,:)];
ME_Pl = [ME_Pl;...
         ME_Pl(ind_add,:)];
REE_Ol = [REE_Ol;...
          REE_Ol(ind_add,:)];
REE_Cpx = [REE_Cpx;...
           REE_Cpx(ind_add,:)];
REE_Opx = [REE_Opx;...
           REE_Opx(ind_add,:)];
REE_Gt = [REE_Gt;...
          REE_Gt(ind_add,:)];
REE_Sp = [REE_Sp;...
          REE_Sp(ind_add,:)];
REE_Pl = [REE_Pl;...
          REE_Pl(ind_add,:)];

% Remove the repeated values
ME_preDouble    =  ME(:,2);
[~,ia,~] = unique(ME_preDouble);
ME = ME(ia,:);
ME_fluid_solid  = ME_fluid_solid(ia,:);
ME_solid_TP_w   = ME_solid_TP_w(ia,:);
ME_solid_TP_v   = ME_solid_TP_v(ia,:);
ME_solid_Rho_TP = ME_solid_Rho_TP(ia,:);
MRadi   = MRadi(ia,:);
ME_Ol   = ME_Ol(ia,:);
ME_Cpx  = ME_Cpx(ia,:);
ME_Opx  = ME_Opx(ia,:);
ME_Gt   = ME_Gt(ia,:);
ME_Sp   = ME_Sp(ia,:);
ME_Pl   = ME_Pl(ia,:);
REE_Ol  = REE_Ol(ia,:);
REE_Cpx = REE_Cpx(ia,:);
REE_Opx = REE_Opx(ia,:);
REE_Gt  = REE_Gt(ia,:);
REE_Sp  = REE_Sp(ia,:);
REE_Pl  = REE_Pl(ia,:); 
      
% get the indexes in order
ME_preSorted    =  ME(:,2);
[~, ind_orden] = sort(ME_preSorted,1); 
ME = ME(ind_orden,:); 
ME_fluid_solid  = ME_fluid_solid(ind_orden,:);
ME_solid_TP_w   = ME_solid_TP_w(ind_orden,:);
ME_solid_TP_v   = ME_solid_TP_v(ind_orden,:);
ME_solid_Rho_TP = ME_solid_Rho_TP(ind_orden,:);
MRadi   = MRadi(ind_orden,:);
ME_Ol   = ME_Ol(ind_orden,:);
ME_Cpx  = ME_Cpx(ind_orden,:);
ME_Opx  = ME_Opx(ind_orden,:);
ME_Gt   = ME_Gt(ind_orden,:);
ME_Sp   = ME_Sp(ind_orden,:);
ME_Pl   = ME_Pl(ind_orden,:);
REE_Ol  = REE_Ol(ind_orden,:);
REE_Cpx = REE_Cpx(ind_orden,:);
REE_Opx = REE_Opx(ind_orden,:);
REE_Gt  = REE_Gt(ind_orden,:);
REE_Sp  = REE_Sp(ind_orden,:);
REE_Pl  = REE_Pl(ind_orden,:); 



