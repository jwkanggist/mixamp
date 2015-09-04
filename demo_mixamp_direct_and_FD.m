% -------------------------------------------------------------------------
% Demo  of the MixAMP framework for 2D sparsity separation 
% where we are handling a sparse mixture of 
%         1) the direct sparsity and
%         2) the finite-difference sparsity.
% Jaewook Kang @ GIST-CSNL
% Final update JWKANG 2015, May. (jwkkang@gist.ac.kr)
%---------------------------------------------------------------------------
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
path(path, './solvers');
path(path, './solvers/denoiser');
path(path, './solvers/denoiser/splitBregmanROF_mex');
path(path,genpath(pwd));

disp('%----------------------------------------------------------------------------------------------------------------------%');
disp('% 2D Sparse signal separation via the MixAMP iteration');
disp('% Sparse Mixture : Direct sparsity + Finite difference Sparsity');
disp('% Copyright@Jaewook Kang with CSNL lab in GIST, Republic of Korea.');
disp('%');
disp('% Written by Jaewook Kang, Phd student in GIST-DIC, jwkkang@gist.ac.kr');
disp('% final update 2015 May.');
disp('%');
disp('% The authors gratefully acknowledge the support from Electronic Warfare Research Center ')
disp('% at Gwangju Institute of Science and Technology (GIST), originally funded by Defense Acquisition Program Administration.');
disp('%----------------------------------------------------------------------------------------------------------------------%');
%%  Problem dimension setting 
alpha=0.6;% undersampling ratio M/N %0.5 / 0.05 / 0.8
Delta=1e-10; % noise variance
% Delta=eps;
%% AMP denoiser parameter setting
maxiter=2000;
iter_tol =5e-4;
% damping_factor=0.7;
%  TVconst=1.2;% % for alpha=0.5, q=0.1 cameramen128

damping_factor=0.6;
 TVconst=1.1;

%%  signal generation 
% Xb image loading
Xb = double(imread('cameraman128.tif'));
Xb = Xb (:,:,1);[N,N]=size(Xb);
load Xa_q01_128.mat 
%% compressive measurement gernation  by Gaussian matrix
P = rand(N,N);P = double(P<alpha);M=nnz(P);
A=orth(randn(N)/sqrt(N));% Standard Gaussian measurement matrices
Y=P.*(A*(Xa+Xb)*A.');% generating the noiseless measurement 

%%  MixAMP solving
tstart=tic;
[est_Xa,est_Xb,theta,stop_iter]=...
    solve_MixAMP_direct_and_FD(A,P,Y,TVconst,damping_factor, maxiter, iter_tol,Delta); 
telapsed_MixAMP=toc(tstart);


%% Normalized MSE calculation 
MSE_Xa_MixAMP=norm(est_Xa(:)-Xa(:))^2 / norm(Xa(:))^2;
MSE_Xb_MixAMP=norm(est_Xb(:)-Xb(:))^2 / norm(Xb(:))^2;
MSE_Xa=norm(est_Xa(:)-Xa(:))^2 / N^2;
MSE_Xb=norm(est_Xb(:)-Xb(:))^2 / N^2;

PSNR_Xb=10*log10(255^2/MSE_Xb);
PSNR_Xa=10*log10(255^2/MSE_Xa);
%% Recovery result Display
disp('<Recovery result>')
disp(sprintf('MixAMP was terminated at t=%d',stop_iter-1))
disp(sprintf('Nomalized MSE wrt to Xa  = %8.7f',MSE_Xa_MixAMP));
disp(sprintf('Nomalized MSE wrt to Xb  = %8.7f',MSE_Xb_MixAMP)); 
disp(sprintf('PSNR of the recovered Xa = %2.4f dB',PSNR_Xa));
disp(sprintf('PSNR of the recovered Xb = %2.4f dB',PSNR_Xb));
disp(sprintf('Running Time of MixAMP  = %8.7f sec',telapsed_MixAMP)); 
disp('%------------------------------------------------------------------------------------------%');


figure(1); clf;
colormap(gray)
subplot(2,3,1);imagesc(Xb+Xa,[0 255]);title('Original Xb+Xa');
subplot(2,3,2);imagesc(inv(A)*Y*inv(A.'),[0 255]);title('A^{-1}*Y*(A^T)^{-1}');
subplot(2,3,3);imagesc(P,[0 1]);title('Sampling matrix P');
subplot(2,3,4);imagesc(est_Xa,[0 200]);title('Recovery of Xa');
subplot(2,3,5);imagesc(est_Xb,[0 255]);title('Recovery of Xb');
subplot(2,3,6);semilogy(1:stop_iter,theta(1:stop_iter));xlabel('Iterations','fontsize',13);title('\theta','fontsize',12);
box on
