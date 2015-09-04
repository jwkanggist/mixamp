%------------------------------------------------------------------------
% function: eta_direct2D
% AMP denoiser for direct sparse signal
% input parameters : 1) [N x 1] object vector :u_in
%                    2) scalar threshold: th
% output           : [N x 1] u_out
% 2015 Mar, written by Jaewook Kang
%-----------------------------------------------------------------------
function u_out = eta_direct2D (u_in, th)

N=length(u_in);
u_out=zeros(N);

set1=find(u_in > th);
set2=find(u_in < -th);

if ~isempty(set1)
    u_out(set1)=u_in(set1) - th;
end

if ~isempty(set2)
    u_out(set2)=u_in(set2) +th;
end

end
