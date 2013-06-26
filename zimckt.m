function zimckt(ckt,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% zimckt: top-level entry of program zimckt
%%
%% - ckt : input circuit file
%% - type: transient simulation type (default 0)
%%         if no TRAN analysis is request,
%%         this option will be ignored
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global_def;
quitall = 0;
oncerun = 0;

while(~quitall)
    clc;
    format long
    
    %type;
    tr_type=0;
    PI = pi;
    
    %type = 0;
    if(nargin == 0)
        Logo;
        
        cmdline = 0;
        while(~cmdline)
            cmdinput = input('cmd>> ','s');
            if(isempty(cmdinput))
                cmdline = 0;
                fprintf('  Type ''help'' to see available inst\n\n');
            elseif(strcmp(cmdinput, 'help') || strcmp(cmdinput, 'HELP'))
                ckthelp;
                cmdline = 0;
            elseif(strcmp(cmdinput, 'quit') || strcmp(cmdinput, 'exit'))
                cmdline = 1;
                quitall = 1;
            elseif(strcmp(cmdinput, 'clc'))
                clc;
                Logo;
            elseif(strcmp(cmdinput, 'clear'))
                clear;
                global_def;
                PI = pi;
                tr_type = 0;
                quitall = 0;
                oncerun = 0;
                cmdline = 0;
            elseif(strcmp(cmdinput, 'run'))
                cmdline = 1;
                cc = input(' - circuit file: ', 's');
                if(isempty(cc))
                    disp('  ** empty file **');
                    cmdline = 0;
                    continue;
                else
                    ckt = cc;
                end
                
            elseif(strcmp(cmdinput, 'ls'))
                eval(cmdinput);
                cmdline = 0;
            else
                fprintf(' Unknown cmd\n');
                fprintf(' Type ''help'' to see available inst\n');
                cmdline = 0;
            end
        end
        
        if(quitall)
            return;
        end
        
    elseif(nargin > 2)
        fprintf('\n Usage: zimckt(''cktname'',type)\n');
        fprintf('  *''type'' is used for transient analysis only:\n    standard       -- 0\n    adaptive ver.1 -- 1\n    adaptive ver.2 -- 2\n');
        fprintf('  *for other analysis (DC, DC SWEEP and AC), \n   you only need to enter "zimckt(''filename'')"\n');
        fprintf('\n ----  Quit  ----\n\n');
        quitall = 1;
        return;
        
    elseif(nargin == 2)
        oncerun = 1;
        
        if(max(size(ckt))>200)
            fprintf(' Check your input argument:\n');
            fprintf('   - correct: zimckt(''file'')\n');
            fprintf('   - wrong:   zimckt(file)\n');
            fprintf('\n ----  Quit  ----\n\n');
            quitall = 1;
            return;
        end
        
        if (type == 0)
            tr_type = 0;
        elseif (type == 1)
            tr_type = 1;
        elseif (type == 2)
            tr_type = 2;
        else
            fprintf('\n Unknown transient option\n');
            fprintf(' Please enter correct transient analysis type:\n    standard       -- 0\n    adaptive ver.1 -- 1\n    adaptive ver.2 -- 2\n');
            fprintf('\n ----  Quit  ----\n\n');
            quitall = 1;
            return;
        end
        
    elseif(nargin == 1)
        oncerun = 1;
        if(max(size(ckt))>200)
            fprintf(' Check your input argument:\n');
            fprintf('   - correct: zimckt(''file'')\n');
            fprintf('   - wrong:   zimckt(file)\n');
            fprintf('\n ----  Quit  ----\n\n');
            quitall = 1;
            return;
        else
            tr_type = 0;
        end
        
    end
    
    if(quitall)
        return;
    end
    
    % init macro
    parser_init;
    
    circuit = ckt;
    
    % Constant settings for nonlinear solver
    %Total tol
    t = [];
    Res_nv = [];
    Res_bi = [];
    LINELM = [];
    NLINELM = [];
    % tolerance for NR iteration
    tol = 1e-6;
    % set Step bounding for nonlinear solver
    step_tol = 0.1;
    % maximum iterations for NR
    max_iter = 500;
    
    % parse the input ckt;
    % reading the netlists: only for getting info of linear elements
    fprintf('\n');
    [parse_ok,LINELM,NLINELM,info,nodes,names,printnv,plotnv,plotbi] = parser(circuit);
    
    if(parse_ok == 0)
        fprintf('\n Parser failed\n');
        if(oncerun)
            fprintf('\n ----  Quit  ----\n\n');
            quitall = 1;
            return;
        else
            cc = input(' Press any key to continue ...');
            continue;
        end
    end
    
    numNodes = length(nodes);
    
    if(info(SIM_) ~= RD_)
        % time step for transient simulation
        Delta_T = info(TSTEP_);
        Res_bi_mos = zeros(size(NLINELM,1),1);
        setup = 1;
        evaluate(numNodes,setup);
        
        % rhs vector and the unknow vector
        F = zeros(dim,1);
        X = zeros(dim,1);
        M = zeros(dim,1);
        % inductors
        X_pre_t = zeros(dim,1);
        X_pre_2t = zeros(dim,1);
        % Jacobian matrix for solving delta X or X
        G = spalloc(dim,dim,0);
        % attach gmin at each node to improve DC convergence
        gmin_init = 1e-10;
    end
    
    fprintf(' >  system dimension is %d-by-%d\n\n', dim, dim);
    %info(SIM_) = DC_;
    if(info(SIM_) == DC_)
        S_type = DC_;
        %sweep=1;
        if(info(METHOD_) == DC_SWEEP)
            s_start = info(SWEEP_START);
            s_end = info(SWEEP_END);
            s_h = info(SWEEP_STEP);
            sweep = zeros((s_end-s_start)/s_h+1,1);
            for k=1:(s_end-s_start)/s_h+1
                sweep(k) = s_start + (k-1)*s_h;
            end
            var = info(SWEEP_ELEM);
            Res_bi = zeros(max(size(sweep)),size(plotbi,1));%! save complex number
            Res_nv = zeros(max(size(sweep)),size(plotnv,1));%
            
            fprintf('**************************************************\n');
            fprintf('   DC SWEEP simulation starts ...\n   ');
            t_dc = cputime;
            for k=1:(s_end-s_start)/s_h+1
                LINELM(var,V_VALUE_) = sweep(k);
                % DC solution should be stored into the Res_dc vector
                Res_dc = dc_sim(tol,step_tol,max_iter,gmin_init);
                
                if(dc_ok == 0)
                    quitall = 1;
                    break;
                end
                
                for l=1:size(plotnv,1)
                    Res_nv(k,l) = Res_dc(plotnv(l));
                end
                for l=1:size(plotbi,1)
                    Res_bi(k,l) = Res_bi_mos(plotbi(l,1));
                end
                
                if(mod(k,30)==0)
                    fprintf('\n   ');
                end
            end
            
            %dc_ok = 0;
            if(dc_ok==1)
                t_dc = cputime - t_dc;
                fprintf('\n     finished!\n');
                fprintf('   CPU time for %d DC SWEEP is %.4f(s) \n', max(size(sweep)), t_dc);
                fprintf('\n**************************************************\n');
                
                output_sol(circuit, nodes, names, plotbi, plotnv, printnv, Res_nv, Res_bi, sweep);
                
                fprintf('**************************************************\n');
                fprintf(' Mission Accomplished!\n');
                %fprintf('\n ----  Quit  ----\n\n');
            else
                %fprintf('\n**************************************************\n');
                fprintf('\n Mission Failed!\n');
                fprintf('\n ----  Quit  ----\n\n');
                quitall = 1;
                return;
            end
            
        else
            % DC solution should be stored into the Res_dc vector
            t_dc = cputime;
            Res_dc = dc_sim(tol,step_tol,max_iter,gmin_init);
            %dc_ok = 0;
            if(dc_ok==1)
                fprintf(' The DC simulation result is\n');
                for i=1:numNodes
                    fprintf('  node (%d): %f(V)\n', nodes(i), Res_dc(i));
                end
                fprintf('**************************************************\n');
                fprintf(' Mission Accomplished!\n');
                %fprintf('\n ----  Quit  ----\n\n');
            else
                %fprintf('**************************************************\n');
                fprintf('\n Mission Failed!\n');
                fprintf('\n ----  Quit  ----\n\n');
                quitall = 1;
                return;
            end
            
        end
        
    elseif(info(SIM_) == TRAN_)
        disp('**************************************************');
        disp('  Do transient analysis, enter TR type');
        disp('    1: adaptive method');
        disp('    2: adaptive method');
        disp('    default: type ''Enter''');
        cc = input('  [default: normal]: ');
        if(isempty(cc))
            tr_type = 0;
        else
            tr_type = cc;
        end
        
        S_type = DC_;
        Res_dc = dc_sim(tol,step_tol,max_iter,gmin_init);
        
        if(dc_ok==0)
            fprintf('\n Mission Failed!\n');
            fprintf('\n ----  Quit  ----\n\n');
            quitall = 1;
            return;
        end
        
        S_type = TRAN_;
        T_tot = info(TSTOP_);
        T_step = info(TSTEP_);
        T_step_init = T_step;
        % standard or adaptive transient simulation
        max_iter_tr = 100;
        % Set initial solution for transient analysis
        X_pre_t = Res_dc;
        % The interested Mosfet currents Ids and node voltages should be stored into Res_bi and Res_nv vectors
        
        if (tr_type == 0)
            % Do standard TRAN analysis using traditional method
            [tr_ok,Res_bi,Res_nv,t] = tr_sim(T_tot,T_step,tol,step_tol,max_iter_tr);
        elseif (tr_type == 1)
            % Do adaptive TRAN analysis using adaptive method
            [tr_ok,Res_bi,Res_nv,t] = tr_simadp(T_tot,T_step_init,tol,step_tol,max_iter_tr);
        elseif (tr_type == 2)
            % Do adaptive TRAN analysis using adaptive method
            [tr_ok,Res_bi,Res_nv,t] = tr_simadp2(T_tot,T_step_init,tol,step_tol,max_iter_tr);
        end
        
        % plot and print the solution
        if(tr_ok==1)
            % DC solution should be stored into the Res_dc vector
            output_sol(circuit, nodes, names, plotbi, plotnv, printnv, Res_nv, Res_bi, t);
            
            fprintf('**************************************************\n');
            fprintf(' Mission Accomplished!\n');
            %fprintf('\n ----  Quit  ----\n\n');
        else
            fprintf('**************************************************\n');
            fprintf('\n Mission Failed!\n');
            fprintf('\n ----  Quit  ----\n\n');
            quitall = 1;
            return;
        end
        
    elseif(info(SIM_) == AC_)
        % find dc operating point
        S_type = DC_;
        % DC solution should be stored into the Res_dc vector
        dc_point = dc_sim(tol,step_tol,max_iter,gmin_init);
        
        if(dc_ok==0)
            fprintf('\n Mission Failed!\n');
            fprintf('\n ----  Quit  ----\n\n');
            quitall = 1;
            return;
        end
        
        S_type = AC_;
        
        f_h = info(AC_PPD_);
        fstart = info(AC_FSTART_);
        fend = info(AC_FSTOP_);
        ac_type = info(METHOD_);
        
        [Res_nv, Res_bi, freq] = ac_sim(dc_point,ac_type,f_h,fstart,fend);
        
        if(ac_ok==1)
            output_sol(circuit, nodes, names, plotbi, plotnv, printnv, Res_nv, Res_bi, freq);
            
            fprintf('**************************************************\n');
            fprintf(' Mission Accomplished!\n');
            %fprintf('\n ----  Quit  ----\n\n');
        else
            %fprintf('**************************************************\n');
            fprintf('\n Mission Failed!\n');
            fprintf('\n ----  Quit  ----\n\n');
            quitall = 1;
            return;
        end
        
    elseif(info(SIM_) == RD_)
        S_type = RD_;
        
        if(info(METHOD_) == PRIMA_)
            orders = 10;
            exps = 0;
            port_ind1 = 1;
            port_ind2 = 1;
            plotf = 1;
            
            %% PRIMA
            [rd_ok] = prima(circuit, names, nodes, orders, exps, port_ind1, port_ind2, plotf);
            
            if(rd_ok)
                fprintf('\n**************************************************\n');
                fprintf(' Mission Accomplished!\n');
                %fprintf('\n ----  Quit  ----\n\n');
            else
                fprintf('\n Mission Failed!\n');
                fprintf('\n ----  Quit  ----\n\n');
                quitall = 1;
                return;
            end
            
        else
            fprintf('\n Unknown MOR method\n');
            %fprintf(' ----  Quit  ----\n\n');
            oncerun = 1;
            quitall = 1;
        end
    else
        fprintf('\n Unknown Sim type\n');
        %fprintf(' ----  Quit  ----\n\n');
        oncerun = 1;
        quitall = 1;
    end
    
    if(oncerun)
        quitall = 1;
        break;
    else
        oncerun = 0;
        quitall = 0;
        cmdinput = input('\n Press any key to continue ...\n', 's');
    end
    
end

fprintf(' ----  Quit  ----\n\n');
return
end
%% end of function zimckt


function Logo()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Logo: print simple help information
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('===============================================\n');
fprintf('     Zimckt: Circuit Simulator for MATLAB      \n');
fprintf('             by xueqian 2013 \n');
fprintf('===============================================\n');
fprintf('\n   Analysis:');
fprintf('\n     * AC');
fprintf('\n     * DC & DC SWEEP');
fprintf('\n     * TRANS (w/ adaptive step control)');
fprintf('\n     * MOR: PRIMA');
fprintf('\n   Sources types available:');
fprintf('\n     * AC, DC, PWL, SIN and PULSE\n');
fprintf('\n---------------------------------------------\n');
fprintf('\n   More details: please read tutorial');
fprintf('\n   Type ''help'' to see instructions\n');
fprintf('\n---------------------------------------------\n\n');
end
%% end of function Logo