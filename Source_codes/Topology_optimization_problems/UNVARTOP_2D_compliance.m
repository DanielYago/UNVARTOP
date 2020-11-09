%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% A 200 line topology optimization code using the Unsmooth Variational Topology method %
%   for minimum compliance problems                                                    %
%   Paper's Title                                                                      %
%		Topology Optimization using the UNsmooth VARiational Topology OPtimization 	   %
%			(UNVARTOP) method: an educational implementation in Matlab                 %
%   Authors                                                                            %
%       Daniel Yago, Juan Carlos Cante, Oriol Lloberas-Valls, Javier Oliver            %
% Versions                                                                             %
%   26/3/2020 - Submission version (minimum compliance)                                %
% GitHub repository																	   %
%	This code and other modifications of it can be downloaded from the Website:		   %
%		https://github.com/DanielYago/UNVARTOP	. Please send your comments to		   %
%		the author: daniel.yago@upc.edu								   				   %
% Copyright																			   %
%	Copyright (C) 2020  Daniel Yago													   %
%																					   %
%	This program is free software: you can redistribute it and/or modify it under the  %
%	terms of the GNU General Public License as published by the Free Software 		   %
%	Foundation, either version 3 of the License, or (at your option) any later version.%
%																					   %
%	This program is distributed in the hope that it will be useful, but WITHOUT ANY    %
%	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    %
%	PARTICULAR PURPOSE.  See the GNU General Public License for more details.          %
%																					   %
%	You should have received a copy of the GNU General Public License along with this  %
%	program.  If not, see <https://www.gnu.org/licenses/>.							   %
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [iter,J] = UNVARTOP_2D_compliance (nelx,nely,nsteps,Vol0,Vol,k,tau)
n_dim = 2; n_unkn = 2; n_nodes = 4; n_gauss = 4; n = (nelx+1)*(nely+1); h_e = 1; alpha0 = 1e-3;
iter_max_step = 20; iter_min_step = 4; iter_max = 500;
opt = struct('Plot_top_iso',1,'Plot_vol_step',1,'EdgeColor','none','solver_Lap','direct');
%% Vector for assembling matrices
[X,Y] = meshgrid(0:nelx,nely:-1:0); coord = [X(:),Y(:)]; clear X Y
nodenrs = reshape(1:n,1+nely,1+nelx);
nodeVec = reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1); clear nodenrs;
connect = nodeVec+[0 nely+[1 0] -1]; clear nodeVec;
%% Loads and boundary setting for Cantilever beam
F = sparse(n_unkn*n,1);
U = zeros(n_unkn*n,1);
F(n_unkn*find(coord(:,2)==0 & coord(:,1)==nelx),1) = -0.01*nelx;
fixed_dofs = reshape(n_unkn*find(coord(:,1)==0)+(-n_unkn+1:0),1,[]);
active_node = []; passive_node = [];
free_dofs = setdiff(1:(n_unkn*n),fixed_dofs);
U(fixed_dofs,:) = 0;
%% Parameter definition
m = 5; E0 = 1; alpha = 1e-6; beta = nthroot(alpha,m); nu = 0.3;
%% Prepare animation
psi_vec = zeros(size(coord,1),nsteps+1);
chi_vec = zeros(size(connect,1),nsteps+1);
U_vec = zeros(n_unkn*size(coord,1),1,nsteps+1);
%% Finite element analysis preparation
[posgp4,W4] = gauss_points(n_gauss);
[posgp1,W1] = gauss_points(1);
[DE] = D_matrix_stress(E0,nu);
KE = zeros(n_nodes*n_unkn,n_nodes*n_unkn);
KE_i = zeros(n_nodes*n_unkn,n_nodes*n_unkn,n_gauss);
for i=1:n_gauss
	[BE,Det_Jacobian] = B_matrix(posgp4(:,i),n_unkn,n_nodes);
	KE_i(:,:,i) = BE'*DE*BE;
	KE = KE + KE_i(:,:,i)*Det_Jacobian*W4(i);
end
[BE_cut,Det_Jacobian_cut] = B_matrix(posgp1,n_unkn,n_nodes);
K_cut = BE_cut'*DE*BE_cut;
KE_cut = K_cut*Det_Jacobian_cut*W1(1);
edofMat = kron(connect,n_unkn*ones(1,n_unkn)) + repmat(1-n_unkn:0,1,n_nodes);
iK = reshape(kron(edofMat,ones(n_nodes*n_unkn,1))',(n_nodes*n_unkn)^2*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,n_nodes*n_unkn))',(n_nodes*n_unkn)^2*nelx*nely,1);
%% Laplacian filter preparation
KE_Lap = 1/6* [ 4 -1 -2 -1;-1 4 -1 -2;-2 -1 4 -1;-1 -2 -1 4];
ME_Lap = 1/36*[ 4  2  1  2; 2 4  2  1; 1  2 4  2; 2  1  2 4];
KE_Lap = ME_Lap + (tau*h_e).^2*KE_Lap;
i_KF = reshape(kron(connect,ones(n_nodes,1))',n_nodes^2*nelx*nely,1);
j_KF = reshape(kron(connect,ones(1,n_nodes))',n_nodes^2*nelx*nely,1);
s_KF = reshape(KE_Lap(:)*ones(1,nelx*nely),n_nodes^2*nelx*nely,1); clear KE_Lap ME_Lap;
K_Lap = sparse(i_KF,j_KF,s_KF);
if strcmp(opt.solver_Lap,'direct'); LF = chol(K_Lap,'lower'); clear K_Lap i_KF j_KF s_KF;
else; LF = ichol(K_Lap, struct('type','ict','droptol',1e-3,'diagcomp',0.1)); clear i_KF j_KF s_KF; end
i_xi = reshape(connect',n_nodes*nelx*nely,1);
N_T = N_matrix(posgp4).*W4/4;
%% Loop over steps
psi = alpha0*ones(n,1); psi(passive_node) = -alpha0; psi(active_node) = alpha0; psi_vec(:,1)=psi;
[~,chi] = compute_volume (psi,connect); chi0_step = chi; chi_vec(:,1) = chi';
% Initialize variables
iter = 1; J_vec = []; vol_vec = []; lambda_vec = 0; lambda = 0; fhandle6 = [];
[fhandle2,ohandle2] = plot_isosurface([],[],0,psi,coord,connect,1,opt);
for i_step = 1:nsteps
	[t_ref] = set_reference_volume(i_step,Vol0,Vol,nsteps,k);
	% Main loop by steps
	Tol_chi = 1;
	Tol_lambda = 1;
	iter_step = 1;
	while (((Tol_chi>1e-1 || Tol_lambda>1e-1) && iter_step<iter_max_step) || iter_step<=iter_min_step)
		% FE-analysis
		[K] = assmebly_stiff_mat (chi,KE,KE_cut,beta,m,iK,jK,n_unkn,nelx,nely);
		U(free_dofs,:) = K(free_dofs,free_dofs) \ (F(free_dofs,:) - K(free_dofs,fixed_dofs)*U(fixed_dofs,:));
		if iter == 1; U_vec(:,:,1)=U; J_ref = full(abs(sum(sum(F.*U,1),2))); end; J = full(sum(sum(F.*U,1),2))/J_ref;
		% Calculate sensitivities
		Energy = zeros(n_gauss,nelx*nely);
		id = chi==1|chi==0; int_chi = interp_property (m,m-1,beta,chi(id));
		u_e = reshape(U(edofMat(id,:)',1),n_nodes*n_unkn,[]); w_e = u_e;
		for i=1:n_gauss; Energy(i,id) = sum(w_e.*(KE_i(:,:,i)*u_e),1); end
		Energy(:,id) = int_chi.*Energy(:,id);
		id = ~id; int_chi = interp_property (m,m-1,beta,chi(id));
		u_e = reshape(U(edofMat(id,:)',1),n_nodes*n_unkn,[]); w_e = u_e;
		Energy(:,id) = repmat(int_chi.*sum(w_e.*(K_cut*u_e),1),n_gauss,1);
		if iter == 1; xi_shift = min(0,min(Energy(:))); xi_norm = max(range(Energy(:)),max(Energy(:))); end
		% Apply Laplacian regularization
		xi_int = N_T*(Energy-xi_shift*chi)/xi_norm;
		if strcmp(opt.solver_Lap,'direct')
			xi = LF'\(LF\accumarray(i_xi,xi_int(:),[n 1]));
		else
			[xi,flag] = minres(K_Lap,accumarray(i_xi,xi_int(:),[n 1]),1e-6,500,LF,LF'); assert(flag == 0);
		end
		% Compute topology
		[lambda,chi_n,psi,vol] = find_volume (xi,connect,active_node,passive_node,t_ref,lambda,alpha0);
		lambda_vec = [lambda_vec,lambda];
		Tol_lambda = (lambda_vec(iter+1)-lambda_vec(iter))/lambda_vec(iter+1);
		% Plot topology
		[fhandle2,ohandle2] = plot_isosurface(fhandle2,ohandle2,iter,psi,coord,connect,J,opt);
		% Update variables
		Tol_chi = sqrt(sum((chi-chi_n).^2))/sqrt(sum(chi0_step.^2));
		chi = chi_n;
		fprintf(' Step:%5i It.:%5i Obj.:%11.4f Vol.:%7.3f \n',i_step,iter_step,J,vol);
		iter_step = iter_step+1; iter = iter+1;
		drawnow;
	end
	chi0_step = chi;
	if J<10
		[fhandle6,J_vec,vol_vec] = plot_volume_iter(fhandle6,i_step,J_vec,J,vol_vec,vol,opt.Plot_vol_step,6,'Cost function Step','#step','+-b');
		psi_vec(:,i_step+1)=psi; chi_vec(:,i_step+1)=chi'; U_vec(:,:,i_step+1)=U;
	end
	if iter_step >= iter_max_step; warning('UNVARTOP_2D_compliance:Max_iter_step','Maximum number of in-step iterations achieved.'); break; end
	if iter > iter_max; warning('UNVARTOP_2D_compliance:Max_iter','Maximum number of iterations achieved.'); break; end
end
%% Animation
Topology_evolution(coord,connect,[Vol0,vol_vec],psi_vec,chi_vec,U_vec);
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [posgp,W] = gauss_points (n_gauss)
if n_gauss==1; s = 0; w = 2; else; s = sqrt(3)/3*[-1 1]; w = [1 1]; end
[s,t] = meshgrid(s,s); posgp = [s(:) t(:)]';
W=w'*w; W=W(:)';
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [N] = N_matrix(posgp)
N = 0.25*(1+[-1 1 1 -1]'*posgp(1,:)).*(1+[-1 -1 1 1]'*posgp(2,:));
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [BE,Det_Jacobian,cart_deriv] = B_matrix(posgp,n_unkn,n_nodes)
dshape = 0.25*[-1 -1;1 -1;1 1;-1 1]'.* flip(1+[-1 -1;1 -1;1 1;-1 1]'.*posgp,1);
Jacobian_mat = dshape*[0 0;1 0;1 1;0 1];
Det_Jacobian = det(Jacobian_mat);
cart_deriv = Jacobian_mat\dshape;
BE = zeros(3,n_unkn*n_nodes); BE([1 3],1:n_unkn:end) = cart_deriv; BE([3 2],2:n_unkn:end) = cart_deriv;
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [DE] = D_matrix_stress(E,nu) %Planestress
DE = E/(1-nu^2)*[1 nu 0; nu 1 0;0 0 (1-nu)/2];
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [coeff] = interp_property (m,n,beta,chi)
coeff = chi + (1-chi).*beta;
coeff = double(m==n).*coeff.^m + double(m~=n).*m*coeff.^n*(1-beta);
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [K] = assmebly_stiff_mat (chi,KE,KE_cut,beta,m,iK,jK,n_unkn,nelx,nely)
sK = interp_property(m,m,beta,chi).*KE(:); sK(:,chi~=1&chi~=0) = interp_property(m,m,beta,chi(chi~=1&chi~=0)).*KE_cut(:);
K = sparse(iK,jK,sK,n_unkn*(1+nelx)*(1+nely),n_unkn*(1+nelx)*(1+nely)); K = (K+K')/2; clear sK;
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [lambda,chi,psi,vol,Tol_constr] = find_volume (xi,connect,active_node,passive_node,t_ref,lambda,alpha0)
l1 = min(xi); c1 = t_ref; l2 = max(xi); c2 = t_ref-1; Tol_constr = 1; iter=1;
if lambda>l1 && lambda<l2
	[chi,psi,vol,l1,l2,c1,c2,Tol_constr] = compute_volume_lambda(xi,connect,active_node,passive_node,t_ref,lambda,l1,l2,c1,c2,alpha0);
end
while (abs(Tol_constr)>1e-4) && iter<1000
	lambda = 0.5*(l1+l2);
	[chi,psi,vol,l1,l2,c1,c2,Tol_constr] = compute_volume_lambda(xi,connect,active_node,passive_node,t_ref,lambda,l1,l2,c1,c2,alpha0);
	iter = iter + 1;
end
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [chi,psi,vol,l1,l2,c1,c2,Tol_constr] = compute_volume_lambda(xi,connect,active_node,passive_node,t_ref,lambda,l1,l2,c1,c2,alpha0)
psi = xi - lambda; psi(passive_node) = -alpha0; psi(active_node) = alpha0;
[vol,chi] = compute_volume (psi,connect);
Tol_constr = -(vol-t_ref);
if Tol_constr > 0, l1 = lambda; c1 = Tol_constr; else; l2 = lambda; c2 = Tol_constr; end
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [volume,chi] = compute_volume (psi,connect)
P = [-1 -1;1 -1;1 1;-1 1]; dvol = 1/4;
s = [-0.9324695142031521 -0.6612093864662645 -0.2386191860831969 0.2386191860831969 0.6612093864662645 0.9324695142031521]; [s,t] = meshgrid(s,s);
w = [ 0.1713244923791704  0.3607615730481386  0.4679139345726910 0.4679139345726910 0.3607615730481386 0.1713244923791704]; W=w'*w; W=W(:)';
psi_n = psi(connect); chi = sum((sign(psi_n)+1),2)'/8; id = chi~=1&chi~=0;
phi_x = psi_n(id,:)*((1+P(:,1)*s(:)').*(1+P(:,2)*t(:)')/4);
chi(1,id) = (W*(phi_x>0)'+ 0.5*W*(phi_x==0)')*dvol;
volume = 1 - sum(chi) / size(connect,1);
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [fig_handle,obj_handle] = plot_isosurface(fig_handle,obj_handle,iter,psi,coord,connect,J,opt)
if opt.Plot_top_iso
	if iter==0; fig_handle = figure(2); set(fig_handle,'Name','Topology'); caxis([-1 1]); colormap(flip(gray(2)));
		axis equal tight; xlabel('x'); ylabel('y'); title(['J = ',num2str(J)]);
		obj_handle = patch('Vertices',coord,'Faces',connect,'FaceVertexCData',psi,'EdgeColor',opt.EdgeColor,'FaceColor','interp');
	else
		set(0, 'CurrentFigure', fig_handle); set(get(gca,'Title'),'String',['J = ',num2str(J)]);
		set(obj_handle,'FaceVertexCData',psi);
	end
end
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [fig_handle,J_vec,vol_vec] = plot_volume_iter(fig_handle,iter,J_vec,J,vol_vec,vol,opt_plot,fig_num,fig_name,xlabel_name,linestyle)
J_vec = [J_vec,J]; vol_vec = [vol_vec,vol];
if opt_plot
	if iter==1; fig_handle = figure(fig_num); set(fig_handle,'Name',fig_name);
		subplot(2,1,1); plot(J_vec,linestyle);   ylabel('$\mathcal{J}_\chi$','Interp','Latex'); xlabel(xlabel_name); grid; grid minor;
		subplot(2,1,2); plot(vol_vec,linestyle); ylabel('|\Omega^-|');                          xlabel(xlabel_name); grid; grid minor;
	else; set(0, 'CurrentFigure', fig_handle);
		subplot(2,1,1); set(findobj(gca,'Type','line'),'Xdata',1:numel(J_vec),'YData',J_vec);
		subplot(2,1,2); set(findobj(gca,'Type','line'),'Xdata',1:numel(vol_vec),'YData',vol_vec);
	end
end
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
function [vol] = set_reference_volume(iter,Vol0,Volf,nsteps,k)
if k==0; vol=Vol0+(Volf-Vol0)/nsteps*iter;
else; C1=(Vol0-Volf)/(1-exp(k)); C2=Vol0-C1; vol=C1*exp(k*iter/nsteps)+C2; end