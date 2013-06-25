function [Res_bi,Res_nv,t] = tr_simadp2(T_tot,T_step,tol,step_tol,max_iter)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Adaptive Transient simulation using adaptive step control
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
global plotbi plotnv T_ X X_pre_t Res_bi_mos
global Delta_T numNodes G F tr_ok

fprintf('**************************************************\n');
fprintf('   Adaptive TRAN simulation Ver.2 starting...\n   ');
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
            fprintf(' Error: can not converge at time step %f\n', T_);
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
    
    % check the variation of srcs in the following 2~3 steps
    % if srcs keep unchanged, then increase time step
    % otherwise, keep the initial time step
    [option, activity] = checksrc(T_, T_step);
    
    tflag=0;
    if (T_ < T_tot && T_ + 2*T_step > T_tot)
        tflag=1;
    end
    
    if (activity==1 && tflag==0)
        if(option == 1)
            T_= T_ + 2 * T_step;
            Delta_T = 2 * T_step;
        elseif(option == 2)
            T_= T_ + 4 * T_step;
            Delta_T = 4 * T_step;
        else
            %disp('  checksrc error: conflicting activity and option');
            T_= T_ + T_step;
            Delta_T = T_step;
            %tr_ok=0;
            %return
        end
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


function [option, activity] = checksrc(T_, T_step)
global LINELM V_ TYPE_
global PWL_ V_TYPE_ V_POINTS_ V_VALUE_

numLINE = size(LINELM,1);
option = 10; % 1: 2steps; 2: 3steps.
activity = 1; % 1: active; 2:deactive

for i = 1:numLINE
    switch LINELM(i, TYPE_)
        case V_
            
            if(LINELM(i,V_TYPE_) == PWL_)
                pwl_times=[];
                pwl_vals=[];
                slope=0;
                numpts=LINELM(i,V_POINTS_)+1;
                startv=LINELM(i,V_VALUE_);
                pwl_times = zeros(numpts,1);
                pwl_vals = zeros(numpts,1);
                pwl_vals(1) = startv;
                for j = 1: numpts-1,
                    pwl_times(j+1)= LINELM(i,V_POINTS_ + 2*j-1) ;
                    pwl_vals(j+1)= LINELM(i,V_POINTS_ + 2*j);
                end
                
                find=0;
                j=1;
                while (T_ >= pwl_times(j)) && (T_<=pwl_times(numpts) && (j<numpts))
                    if(T_ == pwl_times(j)),
                        v = pwl_vals(j);
                        find=1;
                        break
                    end
                    j=j+1;
                end
                
                if (j<=numpts) && (find==0) && (T_<= pwl_times(numpts))
                    slope = (pwl_vals(j) - pwl_vals(j-1))/(pwl_times(j) - pwl_times(j-1));
                    v = pwl_vals(j-1) + slope * (T_- pwl_times(j-1));
                elseif (find==0)
                    v = pwl_vals(numpts);
                end
                
                % evaluate future T + 1 * t_step
                find=0;
                j=1;
                while (T_+T_step >= pwl_times(j)) && (T_+T_step<=pwl_times(numpts) && (j<numpts))
                    if(T_+T_step == pwl_times(j)),
                        v1 = pwl_vals(j);
                        find=1;
                        break
                    end
                    j=j+1;
                end
                
                if (j<=numpts) && (find==0) && (T_+T_step<= pwl_times(numpts))
                    slope = (pwl_vals(j) - pwl_vals(j-1))/(pwl_times(j) - pwl_times(j-1));
                    v1 = pwl_vals(j-1) + slope * (T_+T_step- pwl_times(j-1));
                elseif (find==0)
                    v1 = pwl_vals(numpts);
                end
                
                if(v1 == v)
                    
                    % evaluate future T + 2 * t_step
                    find=0;
                    j=1;
                    while (T_+2*T_step >= pwl_times(j)) && (T_+2*T_step<=pwl_times(numpts) && (j<numpts))
                        if(T_+2*T_step == pwl_times(j)),
                            v2 = pwl_vals(j);
                            find=1;
                            break
                        end
                        j=j+1;
                    end
                    
                    if (j<=numpts) && (find==0) && (T_+2*T_step<= pwl_times(numpts))
                        slope = (pwl_vals(j) - pwl_vals(j-1))/(pwl_times(j) - pwl_times(j-1));
                        v2 = pwl_vals(j-1) + slope * (T_+2*T_step- pwl_times(j-1));
                    elseif (find==0)
                        v2 = pwl_vals(numpts);
                    end
                    
                    if(v2 == v1)
                        option = min(option, 2);
                        activity = activity & 1;
                    else
                        option = min(option, 1);
                        activity = activity & 1;
                    end
                    
                else
                    option = min(option, 0);
                    activity = activity & 0;
                end
 
            end
            
    end
end

end
