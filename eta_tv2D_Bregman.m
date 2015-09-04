%------------------------------------------------------------------------
% function: eta_tv2D_Bregman
% AMP denoiser for finite-difference sparse signal
% input parameters : 1) [N x 1] object vector :u_in
%                    2) mu : control parameter for  || Y - A X A^T||_2^2
%                    3) tol : tolerance for this denoising optimization  
% output           : [N x 1] u_out
% 2015 Mar, written by Jaewook Kang
%-----------------------------------------------------------------------
function u_out = eta_tv2D_Bregman(u_in,mu,tol)
N=length(u_in);
u_out=zeros(N);


%% without knowledge of X0

u_out= SplitBregmanROF(u_in,mu,tol);% by Tom Goldstein' code


end   

