%------------------------------------------------------------------------
% function:  eta_direct2D_prime.m
% Derivative of AMP denoiser for direct sparse signal
% input parameters : 1) [N x 1] object vector :u_in
%                    2) scalar threshold: th
% 2015 Mar, written by Jaewook Kang
%-----------------------------------------------------------------------
function u_out = eta_direct2D_prime (u_in, th)

N=length(u_in);
u_out=zeros(N);

set=find(abs(u_in) > th);

if ~isempty(set)
    u_out(set)=1;
end

end
