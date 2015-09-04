function [mu_a,mu_b,theta,stop_iter]=...
    solve_MixAMP_direct_and_FD(A,P,Y,TVconst,damping_factor, maxiter, iter_tol,Delta)

N=length(A);
M=nnz(P);
%----------------------------------------------------------------------%
mu_a=zeros(N,N);prev_mu_a=zeros(N,N);
mu_b=zeros(N,N);prev_mu_b=zeros(N,N);
r=zeros(N,N)   ;prev_r=zeros(N,N);

theta=zeros(maxiter+1,1);
% initialization
theta(1)=norm(Y,2)^2/M;
prev_r=Y;

%------------------- AMP iteration start-----------------------------%
for t=2:maxiter+1  

    rho_a=A'*prev_r*A+prev_mu_a;
    rho_b=A'*prev_r*A+prev_mu_b;

%   for mixture 2D sparsity
   theta(t)=  Delta ...
              +   damping_factor  * theta(t-1)...
              +(1-damping_factor) * theta(t-1)*( sum(sum(eta_direct2D_prime(rho_a,theta(t-1)) ) )...
                                                        +eta_tv2D_prime(rho_b) )/M;



%     theta(t)=norm(prev_r,2)^2/M;
    
    %signal message update
   mu_a=eta_direct2D      (rho_a, theta(t-1));
   mu_b=eta_tv2D_Bregman  (rho_b,TVconst,0.01);
%    mu_b=eta_tv2D_Bregman  (rho_b,theta(t-1),0.01);

% % residual message update for direct 2D sparsity 
   r        =     damping_factor  * prev_r...
             + (1-damping_factor) * (Y - P.* (A*(mu_a+mu_b)*A') ...
             +  prev_r       *  ( sum(   sum( eta_direct2D_prime(rho_a,theta(t-1)) )  )...
                                             +eta_tv2D_prime(rho_b))/M);
         
%     % Termination condition
%     if  (norm(prev_mu_a - mu_a,'fro' ) / norm(mu_a,'fro') ...
%        + norm(prev_mu_b - mu_b,'fro' ) / norm(mu_b,'fro') ) / 2   < iter_tol 
%         break;
%     end

    if  sqrt( norm(prev_mu_a - mu_a,'fro' )^2  + norm(prev_mu_b - mu_b,'fro' )^2)/...
            sqrt(norm(mu_a,'fro')^2+norm(mu_b,'fro')^2) < iter_tol 
        break;
    end
    
    prev_mu_a = mu_a;
    prev_mu_b = mu_b;
    prev_r    = r;
end
stop_iter=t;

end
