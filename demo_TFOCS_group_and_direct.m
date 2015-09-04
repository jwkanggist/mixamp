%% -------------------------------------------------------------------------
% Demo for the 2D sparse separation via the TFOCS method available at
% http://cvxr.com
% where we are handling a spars mixture of 
%         1) direct sparsity and
%         2) group sparsity.
%-------------------------------------------------------------------------
% This testbench is for the "Split L1 analysis" approach 
%  X= [ Xa; Xb], Ta=[ W | 0],  Tb= group_norm()
%    We solve 
%
%    X_hat = arg min   lambda2*group_norm(X)  + mu/2*|| X - X0 ||_2^2 
%   s.t.
%                       || Y - [A|A] X || <= epsilon
%                       || Ta X ||_1 <= 0.00000
%-------------------------------------------------------------------------
% Jaewook Kang @ GIST-CSNL
% Final update JWKANG 2015, May. (jwkkang@gist.ac.kr)
%% ---------------------------------------------------------------------------
clc
clear all
close all

%Handle random seed
if verLessThan('matlab','7.14')
  defaultStream = RandStream.getDefaultStream;
else
  defaultStream = RandStream.getGlobalStream;
end;

if 1
    savedState = defaultStream.State;
    save random_state.mat savedState;
else
    load random_state.mat
end

defaultStream.State = savedState;

% put key subdirectories in path if not already there
path(path, './images');
path(path,'./etc_functions');
path(path, './solvers');
path(path,genpath(pwd));

disp('%----------------------------------------------------------------------------------------------------------------------%');
disp('% 2D Sparse signal separation via the TFOCS method (available at http://cvxr.com)');
disp('% Sparse Mixture : Direct sparsity + Group Sparsity');
disp('% Copyright@Jaewook Kang with CSNL lab in GIST, Republic of Korea.');
disp('%');
disp('% Written by Jaewook Kang, Phd student in GIST-DIC, jwkkang@gist.ac.kr');
disp('% final update 2015 May');
disp('%');
disp('% The authors gratefully acknowledge the support from Electronic Warfare Research Center ')
disp('% at Gwangju Institute of Science and Technology (GIST), originally funded by Defense Acquisition Program Administration.');
disp('%----------------------------------------------------------------------------------------------------------------------%');
%% Problem dimension setting 
alpha=0.7;% undersampling ratio M/N %0.5 / 0.05 / 0.8
maxiter=1000;
iter_tol =5e-4;
%%  signal generation 
% Xb image loading
Xb = double(imread('qrcode128.png'));
[N,N]=size(Xb);blocksize=N/16;% group sparsity blocksize

load Xa_q005_128.mat 
%% 2D compressive measurement gernation by Gaussian matrix
P = rand(N,N);P = double(P<alpha);% % Sampling matrix generation
M=nnz(P); nz_index=find(vec(P));
A=orth(randn(N)/sqrt(N));  % Standard Gaussian measurement matrices
Y=P.*(A*(Xa+Xb)*A.');      % generating the noiseless measurement
 
%% For TFOCS solving using combining L1 analysis
Xoriginal=[vec(Xa) ;vec(Xb)];
A_kron= kron(A,A); A_kron=A_kron(nz_index, : );
A_TFOCS=[A_kron A_kron];
Yvec=nonzeros(vec(Y));
X0=zeros(2*N^2,1);% initial guess
clear A_kron 

%----------------------------------------------------------------------%
X_ref=Xoriginal+randn(2*N^2,1)*0.001;% 
norm_x_ref      = norm(X_ref);
norm_x_orig     = norm(Xoriginal);
er_ref          = @(x) norm(x-X_ref)/norm_x_ref;
er_signal       = @(x) norm(x-Xoriginal)/norm_x_orig;
resid           = @(x) norm(A_kron*x-Yvec)/norm(Yvec);  
%% Analysis matrix 
% For the direct sparsity (only for Xa part)
ZM=spalloc(N^2,N^2,0);
W1=[speye(N^2) ZM]; 
% For the group sparsity (only for Xb part)
W2=[ZM speye(N^2)]; 
clear ZM
%% control parameters 
% % alpha=0.7, q=0.05, cameraman128
lambda1=0.5;lambda2=1.2 ;
% lambda1=0.8; lambda2=2.5;
epsilon= 1e-10;
mu =5e-10*norm( group_norm(W2*X_ref,blocksize) ,Inf);
% ----------------------------------------------------------------------%
er              = er_ref;  % error with reference solution (from IPM)
obj_ref         = lambda2*group_norm(W2*X_ref,blocksize)+ mu/2*norm(X_ref-X0).^2;
opts            = [];
% opts.normW2     = normW^2;
opts.errFcn     = { @(f,dual,primal) er(primal), ...
                    @(f,dual,primal) obj_ref - f }; 
opts.maxIts     = maxiter;
opts.tol        = iter_tol;
z0  = [];   % we don't have a good guess for the dual
%-------------------------------------------------------------------------%
normA2 = linop_normest( A_TFOCS ).^2;
normW12 = linop_normest( W1 ).^2;
proxScale1 = sqrt( normW12 / normA2 );
prox       = { prox_l2( epsilon ), ...
               proj_linf( proxScale1 * lambda1 )};
%-------------------------------------------------------------------------
%% TFOCS solving
tstart=tic;
Xvec = tfocs_SCD( prox_group_l2_splitL1(blocksize,lambda2,W2), { A_TFOCS, -Yvec ;W1 0}, prox, mu, X0, z0, opts );
telapsed_TFOCS=toc(tstart);

Xa_matout=reshape(Xvec(1:N^2),N,N);
Xb_matout=reshape(Xvec(N^2+1:end),N,N);

MSE_Xa=norm(Xvec(1:N^2)     - Xa(:)   )^2 / N^2;
MSE_Xb=norm(Xvec(N^2+1:end) - Xb(:)   )^2 / N^2;
PSNR_Xa=10*log10(255^2/MSE_Xa);
PSNR_Xb=10*log10(255^2/MSE_Xb);

%% Recovery result Display
disp('<Recovery result>')
disp(sprintf('PSNR of the recovered Xa = %2.4f dB',PSNR_Xa)); 
disp(sprintf('PSNR of the recovered Xb = %2.4f dB',PSNR_Xb)); 
disp(sprintf('Running Time of TFOCS  = %2.7f sec',telapsed_TFOCS)); 
disp('%------------------------------------------------------------------------------------------%');


figure(1); clf;
colormap(gray)
subplot(2,3,1);imagesc(Xb+Xa,[0 255]);title('Original Xb+Xa');
subplot(2,3,2);imagesc(inv(A)*Y*inv(A.'),[0 255]);title('A^{-1}*Y*(A^T)^{-1}');
subplot(2,3,3);imagesc(P,[0 1]);title('Sampling matrix P');
subplot(2,3,4);imagesc(Xa_matout, [0 200]);title('Recovery of Xa');
subplot(2,3,5);imagesc(Xb_matout,[0 255]);title('Recovery of Xb');
% subplot(2,3,6);semilogy(1:outer_maxiter,theta(1:outer_maxiter));xlabel('Iterations','fontsize',13);title('\theta','fontsize',12);
box on
