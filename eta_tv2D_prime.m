%------------------------------------------------------------------------
% function:  eta_tv2D_prime.m
% Derivative of AMP denoiser for finite difference sparse signal
% input parameters : 1) [N x 1] object vector :u_in
% 2015 Mar, written by Jaewook Kang
%-----------------------------------------------------------------------
function u_out = eta_tv2D_prime (u_in)

N=length(u_in);

Dx_h=diff(u_in);
Dx_v=diff(u_in.')';


u_out = ( nnz( Dx_h  )+nnz( Dx_v))/N^2;
% u_out =  length(nonzeros( Dx_h(:,1)  ))+length(nonzeros( Dx_v(1,:)));
end
