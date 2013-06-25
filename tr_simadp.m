function [Res_bi,Res_nv,t] = tr_simadp(T_tot,T_step,tol,step_tol,max_iter)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Adaptive Transient simulation using adaptive step control
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
global plotbi plotnv T_ X X_pre_t Res_bi_mos
global Delta_T numNodes G F tr_ok

fprintf('**************************************************\n');
fprintf('   Adaptive TRAN simulation Ver.1 starting...\n   ');
% attach gmin at each node to improve convergence
gmin = 1e-12;
setup=0;
%num_t_pts = ceil(T_tot/T_step);
%Delta_T = T_step;
t=[];
%     save time points
% for i=1:num_t_pts
%     t(i) = T_step * (i-1);
% end
%     total iterations of tr
tot_iter = 0;
%  intialize Res_bi and Res_nv
Res_bi = [];%zeros(num_t_pts,size(plotbi,1));
Res_nv = [];%zeros(num_t_pts,size(plotnv,1));
t_tr = cputime;

% Put your codes here
len = length(X);

%num_iter_tr = zeros(num_t_pts,1);
T_ = 0;
n = 1;
%for n=1:num_t_pts
while(T_ <= T_tot)
    X_pre = X;
    sc = 1;
    
    iter=0;
    while (sc > tol)
        evaluate(numNodes,setup);
        for i=1:numNodes
            G(i,i) = G(i,i) + gmin;
        end
        
        X = G\F;
        for i=1:len
            if (X(i) - X_pre(i) > step_tol)
                X(i) = X_pre(i) + step_tol;
            elseif (X(i) - X_pre(i) < -1 * step_tol)
                X(i) = X_pre(i) - step_tol;
            end
        end
        
        tot_iter = tot_iter + 1;
        Temp = X - X_pre;
        X_pre = X;
        sc = norm(Temp,inf);
        iter = iter + 1;
        if(iter > max_iter)
            fprintf(' Error: can not converge at time step %e\n', T_);
            fprintf('   Transient simulation failed\n');
            fprintf('**************************************************\n');
            return
        end
    end
    X_pre_t = X;
    t(n) = T_;
    
    if(mod(n,4)==0)
        fprintf('.');
    end
    %fprintf('.');
    if(mod(n,100)==0)
        fprintf('\n   ');
    end

    tflag=0;
    if (T_ < T_tot && T_ + 2*T_step > T_tot)
        tflag=1;
    end
    
    if (iter <= 2 && tflag==0)
        T_= T_ + 2 * T_step;
        Delta_T = 2 * T_step;
    else
        T_= T_ + T_step;
        Delta_T = T_step;
    end

    
    for j=1:size(plotnv,1)
        Res_nv(n,j) = X_pre_t(plotnv(j));
    end
    
    for j=1:size(plotbi,1)
        Res_bi(n,j) = Res_bi_mos(plotbi(j,1));
    end
    %Res_bi_mos;
    n = n + 1;
end
tr_ok=1;
t_tr = cputime - t_tr;

fprintf('\n     finished!\n');
fprintf('   (%d) steps for TRAN analysis \n',n);
fprintf('   (%d) N-R iterations for TRAN analysis \n',tot_iter);
fprintf('   CPU time for TRAN analysis is %.4f(s) \n',t_tr);
fprintf('**************************************************\n');

end
