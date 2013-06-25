function [Res_dc] = dc_sim(tol,step_tol,max_iter,gmin_init)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   DC simulation
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
global X G F numNodes dc_ok 
global DC_SWEEP info METHOD_

if info(METHOD_) ~= DC_SWEEP
    fprintf('**************************************************\n');
    fprintf('   DC simulation starting...\n   ');
end
num_iter_dc = 0;

dc_ok = 0;
t_dc = cputime;

setup = 0;
% Put your codes here
len = length(X);
X_pre = X;
Gmin = gmin_init;
sc = 1;

%info(METHOD_)

while Gmin > 1e-20
    X_pre = X;
    evaluate(numNodes,setup);
    for i=1:numNodes
        G(i,i) = G(i,i) + Gmin;
    end
    X = G\F;
    for i=1:len
        if X(i) - X_pre(i) > step_tol
            X(i) = X_pre(i) + step_tol;
            
        elseif X(i)-X_pre(i) < -1 * step_tol
            X(i) = X_pre(i) - step_tol;
        end
    end
    Temp = X - X_pre;
    sc = norm(Temp,inf);
    while sc > tol
        X_pre = X;
        evaluate(numNodes,0);
        for i=1:numNodes
            G(i,i) = G(i,i) + Gmin;
        end

        X = G\F;

        num_iter_dc = num_iter_dc + 1;
        if num_iter_dc >= max_iter
            break;
        end
        for i=1:len
            if X(i) - X_pre(i) > step_tol
                X(i) = X_pre(i) + step_tol;
                
            elseif X(i)-X_pre(i) < -1 * step_tol
                X(i) = X_pre(i) - step_tol;
                
            end
        end
        Temp = X - X_pre;
        sc=norm(Temp,inf);
        
    end
    if num_iter_dc >= max_iter
        break;
    end
    Gmin = Gmin/10;
    
    if info(METHOD_) ~= DC_SWEEP
        fprintf('.');
    end
end


t_dc = cputime - t_dc;
Res_dc = X;
%num_iter_dc = max_iter+1;
if num_iter_dc <= max_iter
    dc_ok = 1;
else
    fprintf('\n   Warning:\n    DC analysis does not converge after %d iterations\n',num_iter_dc);
    fprintf('**************************************************\n');
    return
end;


if dc_ok == 1 && info(METHOD_) ~= DC_SWEEP
    fprintf('\n     finished!\n');
    fprintf('   (%d) N-R iterations for DC analysis \n',num_iter_dc);
    fprintf('   CPU time for DC analysis is %.4f(s) \n',t_dc);
    fprintf('**************************************************\n');
elseif dc_ok == 1 && info(METHOD_) == DC_SWEEP
    fprintf('.');
end


end
