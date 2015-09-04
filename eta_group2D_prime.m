%------------------------------------------------------------------------
% function:   eta_group2D_prime 
% Derivative of AMP denoiser for group sparse signal
% input parameters : 1) [N x 1] object vector :u_in
%                    2) scalar threshold: th
%                    3) block size : blocksize

% 2015 Mar, written by Jaewook Kang
%-----------------------------------------------------------------------
function u_out = eta_group2D_prime(u_in, th, blocksize)

N=length(u_in);
u_out=0;
B=im2col(u_in, [blocksize,blocksize], 'distinct');
blockwise_norm=sqrt(sum(B.^2));

set =find(blockwise_norm > th);

if ~isempty(set)
    u_out=sum( blocksize - (blocksize - 1)*th./blockwise_norm(set));

end
