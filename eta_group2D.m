%------------------------------------------------------------------------
% function: eta_group2D
% AMP denoiser for group sparse signal
% input parameters : 1) [N x 1] object vector :u_in
%                    2) scalar threshold: th
%                    3) block size : blocksize 1by1 2by2  4by4 8by8
% output           : [N x 1] u_out
% 2015 Mar, written by Jaewook Kang
%-----------------------------------------------------------------------
function u_out = eta_group2D (u_in, th, blocksize)

N=length(u_in);
u_out=zeros(N);

B=im2col(u_in, [blocksize,blocksize], 'distinct');


% soft block thresholding
% blockwise_norm=sqrt(sum(B.^2));
blockwise_norm=sum(abs(B));


set1 =find(blockwise_norm < th);
set2 =find(blockwise_norm > th);

B(:,set1)=zeros(blocksize^2,length(set1));
B(:,set2)=repmat(max(1- th./blockwise_norm(set2),0),[blocksize^2 1]).*B(:,set2); 
u_out=col2im(B, [blocksize,blocksize], [N,N], 'distinct');

end
