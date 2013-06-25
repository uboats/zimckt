function [Res_bi,Res_nv,t] = tr_sim(T_tot,T_step,tol,step_tol,max_iter)

%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Transient simulation using fixed time step
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
global plotbi plotnv T_ X X_pre_t X_pre_2t Res_bi_mos I_pre I_eqp
global G F numNodes Delta_T num_C tr_n tr_ok

fprintf('**************************************************\n');
fprintf('   Standard TRAN simulation starting...\n   ');

% attach gmin at each node to improve convergence
gmin = 1e-12;
num_t_pts = ceil(T_tot/T_step);
Delta_T = T_step;
%     save time points
for i=1:num_t_pts
    t(i) = T_step * (i-1);
end
t(num_t_pts)=T_tot;

%     total iterations of tr
tot_iter = 0;
%  intialize Res_bi and Res_nv
Res_bi = zeros(num_t_pts,size(plotbi,1));
Res_nv = zeros(num_t_pts,size(plotnv,1));

tr_ok=0;
setup=0;
% Put your codes here
len = length(X);
I_pre = zeros(1,num_C);
I_eqp = zeros(1,num_C);
X_pre_t = X;
X_pre_2t = X;
T_ = 0;
tr_n = 0;

t_tr = cputime;
%num_iter_tr = zeros(num_t_pts,1);
%T_=0;n=1;
for n=1:num_t_pts
    tr_n = n;
    X_pre = X;
    sc = 1;
    T_ = t(n);
    iter = 0;
    while (sc > tol)
        evaluate(numNodes,setup);
        for i=1:numNodes
            G(i,i) = G(i,i) + gmin;
        end
        
        X = G\F;
        for i=1:len
            if (X(i) - X_pre(i) > step_tol)
                X(i) = X_pre(i) + step_tol;
            elseif (X(i) - X_pre(i)< -1 * step_tol)
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
    X_pre_2t = X_pre_t;
    X_pre_t = X;
    if(mod(n,5) == 0)
        fprintf('.');
    end
    %fprintf('.');
    if(mod(n,200) == 0)
        fprintf('\n   ');
    end
    %size(plotnv, 1)
    %size(Res_nv)
    %size(X_pre_t)
    for j=1:size(plotnv,1)
        %plotnv(j)
        Res_nv(n,j) = X_pre_t(plotnv(j));
    end
    
    for j=1:size(plotbi,1)
        Res_bi(n,j) = Res_bi_mos(plotbi(j,1));
    end
    %Res_bi_mos;
end
tr_ok=1;
t_tr = cputime - t_tr;

fprintf('\n     finished!\n');
fprintf('   (%d) steps for TRAN analysis \n',n);
fprintf('   (%d) N-R iterations for TRAN analysis \n',tot_iter);
fprintf('   CPU time for TRAN analysis is %.4f(s) \n',t_tr);

fprintf('**************************************************\n');

end
