function output_sol(circuit, nodes, names, plotbi, plotnv, printnv, Res_nv, Res_bi, t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% output_sol: print and plot the solutions
%%
%% - circuit: input circuit file
%% - nodes  : node dict map
%% - names  : device name
%% - plotbi : plot branch current
%% - plotnv : plot node voltage
%% - printnv: print node voltage
%% - Res_nv : node voltage result
%% - Res_bi : branch current result
%% - t      : time points
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global S_type DC_ TRAN_ AC_ DEC_ LIN_ DC_SWEEP RD_ PRIMA_
global info METHOD_ SWEEP_DEV

fprintf('\n');

nv=[];
tmpstr = keys(nodes);
tmpval = cell2mat(values(nodes));
for m=1:size(plotnv,1)
    %nv = combp(tmpstr(plotnv(m)),nv,m);
    for ll=1:length(tmpval)
        if(tmpval(ll)==plotnv(m))
            nv = [nv; tmpstr(ll)];
        end
    end
end
bi=[];
for m=1:size(plotbi,1)
    %bi = combp(tmpstr(plotbi(m)),bi,m);
    for ll=1:length(tmpval)
        if(tmpval(ll)==plotbi(m))
            bi = [bi; tmpstr(ll)];
        end
    end
    %bi = [bi; tmpstr(plotbi(m))];
end

if(S_type == TRAN_)
    % Plot the final Mosfet currents as well as node voltages
    if(isempty(plotbi) == 0 && isempty(plotnv) == 0)
        ext = '.binv.fig';
        resultfile = [circuit ext];
        subplot(2,1,1); plot(t,Res_nv);
        xlabel('Time points (s)');ylabel('Amplitude (V)');
        legend(nv)
        
        subplot(2,1,2); plot(t,Res_bi);
        xlabel('Time points (s)');ylabel('Amplitude (Amp)');
        legend(bi)
        saveas(gcf,resultfile,'fig');
        fprintf(' The figure is saved in "%s" \n',resultfile);
    elseif(isempty(plotnv) == 0)
        ext = '.nv.fig';
        resultfile = [circuit ext];
        plot(t,Res_nv);
        xlabel('Time points (s)');ylabel('Amplitude (V)');
        legend(nv)
        saveas(gcf,resultfile,'fig');
        fprintf(' The figure is saved in "%s" \n',resultfile);
    elseif(isempty(plotbi) == 0)
        ext = '.bi.fig';
        resultfile = [circuit ext];
        plot(t,Res_bi);
        xlabel('Time points (s)');ylabel('Amplitude (Amp)');
        legend(bi)
        saveas(gcf,resultfile,'fig');
        fprintf(' The figure is saved in "%s" \n',resultfile);
    end
    
    % Print the required results into files using matlab format
    if(isempty(printnv) == 0)
        ext = '.nv.m';
        resultfile = [circuit ext];
        %size(Res_nv)
        fip = fopen(resultfile, 'w');
        for i=1:size(printnv,1)
            for ll=1:length(tmpval)
                if(tmpval(ll)==printnv(i))
                    fprintf(fip, '%% node: %s\n', char(tmpstr(ll)));
                end
            end            
            fprintf(fip, 'n%d = [\n', i);
            for j=1:size(Res_nv,1)
                fprintf(fip, '%6.4f;\n', Res_nv(j,i));
            end
            fprintf(fip, '];\n\n\n');
        end
        
        fprintf(fip, 't = [\n');
        for i=1:max(size(t))
            fprintf(fip, '%e;\n', t(i));
        end
        fprintf(fip, '];\n\n');
        
        fclose(fip);
        fprintf(' The data is saved in "%s" \n',resultfile);
    end
    
elseif(S_type == AC_)
    % each figure contains two subfigures: amplitude and phase
    nv_amp = zeros(max(size(Res_nv)),size(plotnv,1));
    nv_phase = zeros(max(size(Res_nv)),size(plotnv,1));
    
    bi_amp = zeros(max(size(Res_bi)),size(plotbi,1));
    bi_phase = zeros(max(size(Res_bi)),size(plotbi,1));
    
    for m=1:size(plotnv,1)
        for n=1:max(size(Res_nv))
            nv_amp(n,m) = norm(Res_nv(n,m));
            nv_phase(n,m) = phase(Res_nv(n,m));
        end
    end
    
    for m=1:size(plotbi,1)
        for n=1:max(size(Res_bi))
            bi_amp(n,m) = norm(Res_bi(n,m));
            bi_phase(n,m) = phase(Res_bi(n,m));
        end
    end
    
    if(isempty(plotnv) == 0)
        ext = '.nv.fig';
        resultfile = [circuit ext];
        
        if(info(METHOD_) == LIN_)
            subplot(2,1,1);loglog(t,nv_amp);
            xlabel('Frequence in dec log (Hz)');ylabel('Voltage Gain');
            legend(nv)
            subplot(2,1,2);semilogx(t,nv_phase);
            xlabel('Frequence in dec log (Hz)');ylabel('Phase (degree)');
            legend(nv)
        elseif(info(METHOD_) == DEC_)
            subplot(2,1,1);loglog(t,nv_amp);
            xlabel('Frequence in dec log (Hz)');ylabel('Voltage Gain');
            legend(nv)
            subplot(2,1,2);semilogx(t,nv_phase);
            xlabel('Frequence in dec log (Hz)');ylabel('Phase (degree)');
            legend(nv)
        end
        saveas(gcf,resultfile,'fig');
        fprintf(' The figure is saved in "%s" \n',resultfile);
    end
   
% currently not evaluate current AC
%     if(isempty(plotbi) == 0)
%         ext='.bi.fig';
%         resultfile=[circuit ext];
%         
%         if(info(METHOD_)==LIN_)
%             subplot(2,1,1);loglog(t,bi_amp);
%             xlabel('Frequence in linear (Hz)');ylabel('Current Gain');
%             legend(bi)
%             subplot(2,1,2);semilogx(t,bi_phase);
%             xlabel('Frequence in linear (Hz)');ylabel('Phase (degree)');
%             legend(bi)
%         elseif(info(METHOD_)==DEC_)
%             subplot(2,1,1);loglog(t,bi_amp);
%             xlabel('Frequence in dec log (Hz)');ylabel('Current Gain');
%             legend(bi)
%             subplot(2,1,2);semilogx(t,bi_phase);
%             xlabel('Frequence in dec log (Hz)');ylabel('Phase (degree)');
%             legend(bi)
%         end
%         %saveas(gcf,'node_voltage','fig');
%         saveas(gcf,resultfile,'fig');
%         fprintf(' The figure is saved in "%s" \n',resultfile);
%     end
    
elseif(S_type == DC_)
    
    if(info(METHOD_) == DC_SWEEP)
        if(isempty(plotbi) == 0 && isempty(plotnv) == 0)
            ext = '.binv.fig';
            resultfile = [circuit ext];
            subplot(2,1,1); plot(t,Res_nv);
            xlabel('Time points (s)');ylabel('Amplitude (V)');
            legend(nv)
            
            subplot(2,1,2); plot(t,Res_bi);
            xlabel('Time points (s)');ylabel('Amplitude (Amp)');
            legend(bi)
            saveas(gcf,resultfile,'fig');
            fprintf(' The figure is saved in "%s" \n',resultfile);
        elseif(isempty(plotnv) == 0)
            ext = '.nv.fig';
            resultfile = [circuit ext];
            plot(t,Res_nv);
            xtest = 'Sweep variable ';
            xtest = [xtest names(info(SWEEP_DEV),:)];
            xlabel(xtest);ylabel('Amplitude (V)');
            legend(nv)
            saveas(gcf,resultfile,'fig');
            fprintf(' The figure is saved in "%s" \n',resultfile);
        elseif(isempty(plotbi) == 0)
            ext = '.bi.fig';
            resultfile = [circuit ext];
            plot(t,Res_bi);
            xlabel('Time points (s)');ylabel('Amplitude (Amp)');
            legend(bi)
            saveas(gcf,resultfile,'fig');
            fprintf(' The figure is saved in "%s" \n',resultfile);
        end
    end
elseif(S_type == RD_)
    
    if(info(METHOD_) == PRIMA_)
        if(isempty(plotbi) == 0 && isempty(plotnv) == 0)
            ext = '.binv.fig';
            resultfile = [circuit ext];
            subplot(2,1,1); plot(t,Res_nv);
            xlabel('Time points (s)');ylabel('Amplitude (V)');
            legend(nv)
            
            subplot(2,1,2); plot(t,Res_bi);
            xlabel('Time points (s)');ylabel('Amplitude (Amp)');
            legend(bi)
            saveas(gcf,resultfile,'fig');
            fprintf(' The figure is saved in "%s" \n',resultfile);
        elseif(isempty(plotnv) == 0)
            ext = '.nv.fig';
            resultfile = [circuit ext];
            plot(t,Res_nv);
            xtest = 'Sweep variable ';
            xtest = [xtest names(info(SWEEP_DEV),:)];
            xlabel(xtest);ylabel('Amplitude (V)');
            legend(nv)
            saveas(gcf,resultfile,'fig');
            fprintf(' The figure is saved in "%s" \n',resultfile);
        elseif(isempty(plotbi) == 0)
            ext = '.bi.fig';
            resultfile = [circuit ext];
            plot(t,Res_bi);
            xlabel('Time points (s)');ylabel('Amplitude (Amp)');
            legend(bi)
            saveas(gcf,resultfile,'fig');
            fprintf(' The figure is saved in "%s" \n',resultfile);
        end
    end
end

if(isempty(plotbi) == 1 && isempty(plotnv) == 1)
    fprintf(' No plot or print request\n');
end

fprintf('\n');

end
%% end of function output_sol

function [nv] = combp(A, nv, number_elem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% combp: print space extension
%%
%% - A          : temp vector includes a tracking node name
%% - nv         : vector stores tracking node names 
%% - number_elem: current total number of devices in nv
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_a = size(A,2);
s_e = size(nv,2);
if(s_a > s_e)
    kk = num2str(zeros(size(nv,1),s_a-s_e));
    nv = [kk,nv];
elseif(s_a < s_e)
    kk = num2str(zeros(1,s_e-s_a));
    A = [kk,A];
end
nv(number_elem,:) = A;
end
%% end of function combp