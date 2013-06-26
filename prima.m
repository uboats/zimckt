function [rd_ok] =  prima(circuit, names, nodes, orders, exps, port_ind1, port_ind2, plotf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% prima: apply PRIMA model order reduction
%%
%% - circuit     : input file
%% - orders      : moment matching order for each expansion: one for each
%%                 expansion point in Hz
%% - exps        : expansion points (assuming each is real)
%% - port_ind1/2 : input/output ports
%% - plotf       : if 1, plot full model 
%%                 otherwise, only plot the reduced model
%% - res_r, res_f: contain the port TF data of the reduced order
%%                 model and the full model for later plotting.
%%
%% modified by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global G1_ C1_ B1_ GR_ CR_ BR_ L1_ LR_
fprintf('**************************************************\n');
fprintf('   MOR using PRIMA starts ...\n');

numNodes = length(nodes);
t_rd = cputime;

G1_ = [];
C1_ = [];
B1_ = [];
[stamp_ok, G1_, C1_, B1_, L1_, port_ind, port_name] = stamp_prima(numNodes,1,names);

if(~stamp_ok)
	rd_ok = 0;
	return
end

j_ = sqrt(-1);

% *************
[dim, numports] = size(B1_);

PROJ_PRIMA = [];

for k = 1:length(exps),
    if rem(orders(k), numports) == 0
        border = orders(k)/numports;
    else
        border = ceil(orders(k)/numports);
    end
    
    Gi = G1_ + 2*pi*exps(k)*C1_;
    M0 = Gi\B1_;
    PROJ = orth(M0);
    for i = 1:border-1,
        % Mk = G\C*Xk-1
        M = Gi\(C1_*PROJ(:,(numports*(i-1)+1):(numports*i)));
        [PROJ, r] = qr([PROJ M], 0);
    end
    
    if border * numports > orders(k),
        % extract the first q columns
        %[u, s,v ] = svd(PROJ, 0);
        %PROJ_PRIMA = [PROJ_PRIMA u(:,1:orders(k))];
        dim1 = (border-1)*numports;
        [u,s,v] = svd(PROJ(:, dim1 + 1: border*numports), 0);
        PROJ_PRIMA = [PROJ_PRIMA PROJ(:, 1:dim1) u(:,1:orders(k) - dim1)];
    else
        PROJ_PRIMA = [PROJ_PRIMA PROJ];
    end
end
[PROJ_PRIMA, r] = qr([PROJ_PRIMA], 0);

% calculate reduced system
GR_ = PROJ_PRIMA.'*G1_*PROJ_PRIMA;
CR_ = PROJ_PRIMA.'*C1_*PROJ_PRIMA;
BR_ = PROJ_PRIMA.'*B1_;
LR_ = PROJ_PRIMA.'*L1_;

fprintf('     - Original system size: %d\n', size(PROJ_PRIMA,1));
fprintf('     - Reduced prima-model size: %d\n', size(PROJ_PRIMA,2));
fprintf('\n');

% plot results
w = 2*pi*j_*logspace(4,12,300);
x_axis = logspace(4,12,300);

if(length(port_ind)<1)
    fprintf(' * No port available\n');
    rd_ok = 0;
    return
end

% choose in/out ports
fprintf(' Plot response curve through i/o ports\n');
fprintf(' * %d Ports :[ ', length(port_ind));
for i = 1:length(port_ind)
	%fprintf('%d ', port_ind(i));
    fprintf('%s ', lower(port_name(i,:)));
end
fprintf(']\n');
inp = input(' - choose an input port (?-th): ');
outp = input(' - choose output ports (e.g. [1,3]): ');
if(isempty(inp))
    fprintf('   - default in-port: %d-th\n', 1);
	inport = 1;
else
    if(length(inp)>1)
        fprintf(' * Single input port is restricted, default choose 1st\n');
        inport = 1;
    else
        inport = inp;
    end
end
if(isempty(outp))
	if(length(port_ind)>1)
        fprintf('   - default out-port: %d-th\n', 2);
		outport = 2;
    else
        fprintf('   - default out-port: %d-th\n', 1);
		outport = 1;
	end
else
    if(length(outp)>length(port_ind))
        fprintf('   * exceeds #ports, default choose 1st\n');
        outport = 1;
    else
        outport = outp;
    end
end

nfreqs = length(w);
for i = 1:nfreqs,
    if plotf,
        x = (G1_ + w(i) * C1_)\B1_(:, inport);
        for j = 1: length(outport),
            res_f(j,i) = x(port_ind(outport(j)));
        end
    end
    x = (GR_ + w(i) * CR_)\BR_(:, inport);
    y = LR_' * x;
    for j = 1: length(outport),
        res_r(j, i) = y(outport(j));
    end
end

t_rd = cputime - t_rd;
fprintf('\n     finished!\n');
fprintf('   CPU time for MOR PRIMA is %.4f(s) \n', t_rd);
fprintf('\n**************************************************\n');
% x_axis, res_f, res_r, port_ind2, port_ind1
fprintf('\n');
for j = 1: length(outport),
        figure;
		ext = '.mor.fig';
		circuit = [circuit '.'];
        resultfile = [circuit num2str(outport(j)) ext];
        %loglog(x_axis, abs(res_f(j,:)),'b', x_axis, abs(res_r(j,:)), 'r-.');
        loglog(x_axis, abs(res_r(j,:)), 'b');
        if(plotf)
            hold on
            loglog(x_axis, abs(res_f(j,:)), 'r-.');
            legend('Full Model', 'PRIMA');
        else
            legend('PRIMA');
        end
        str = sprintf('T%d%d Magnitude', outport(j), inport);
        title(str);
		saveas(gcf,resultfile,'fig');
        fprintf(' The %d-th fig is saved in "%s" \n', j, resultfile);
end

rd_ok = 1;

end
%% end of function prima


function [stamp_ok, G, C, B, L, port_ind, port_name] = stamp_prima(numNodes, grd_vsrc, names)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% stamp_prima: Voltage source/inductor stamps follow 
%%              what in the PRIMA paper to ensure the 
%%              passivity of MOR. The inputs can be a 
%%              combination of vsources and isources.
%%
%% - numNodes: number of circuit nodes
%% - port_ind: port information. For vsources: it gives the
%%             corresponding current unknown index; for csources: 
%%             it specifies the current node index (node voltage).
%% - grd_vsrc: if 1, a vsource with dc value being 0 will be treated
%%             as a ground vsource. It will be considered in the G matrix 
%%             but will not be included as part of B matrix(ports).
%%
%% modified by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global LINELM TYPE_ C_ L_ R_ Y_ V_ I_ BRANCH_CUR_IND_
global R_N1_ R_N2_ Y_N1_ Y_N2_ C_N1_ C_N2_ V_N1_ V_N2_ L_N1_ L_N2_
global I_N1_ I_N2_ R_VALUE_ Y_VALUE_ C_VALUE_ L_VALUE_ V_VALUE_ ...
    K_VALUE_ K_ELEM1_ K_ELEM2_ K_ L_CUR_


G = spalloc(numNodes, numNodes, 6*numNodes);
C = spalloc(numNodes, numNodes, 20*numNodes);
%B = zeros(numNodes, 1);
B = [];
port_ind = [];
port_name = [];

stamp_ok = 1;

numElm = size(LINELM,1);
for i = 1:numElm,
    switch LINELM(i, TYPE_)
        case R_ %resistance
            n1 = LINELM(i, R_N1_);
            n2 = LINELM(i, R_N2_);
            g =1/ LINELM(i, R_VALUE_);
            if n1>0 && n2>0,
                G(n1, n1) = G(n1, n1) + g;
                G(n1, n2) = G(n1, n2) - g;
                G(n2, n1) = G(n2, n1) - g;
                G(n2, n2) = G(n2, n2) + g;
            elseif n1>0,
                G(n1, n1) = G(n1, n1) + g;
            elseif n2>0,
                G(n2, n2) = G(n2, n2) + g;
            end
        case Y_ %conductance
            n1 = LINELM(i, Y_N1_);
            n2 = LINELM(i, Y_N2_);
            g  = LINELM(i, Y_VALUE_);
            if n1>0 && n2>0,
                G(n1, n1) = G(n1, n1) + g;
                G(n1, n2) = G(n1, n2) - g;
                G(n2, n1) = G(n2, n1) - g;
                G(n2, n2) = G(n2, n2) + g;
            elseif n1>0,
                G(n1, n1) = G(n1, n1) + g;
            elseif n2>0,
                G(n2, n2) = G(n2, n2) + g;
            end
            
        case L_ %inductance
            n1 = LINELM(i, L_N1_);
            n2 = LINELM(i, L_N2_);
            % need to add inductance current as one variable
            len = length(G);
            C(len+1, len+1) = LINELM(i, L_VALUE_);
            if n1 > 0,
                G(n1, len+1) = 1;
                G(len+1, n1) = -1;
            end
            if n2 > 0,
                G(n2, len+1) = -1;
                G(len+1, n2) = 1;
            end
            % make the dimension agree
            [tmp, num_cols] = size(B);
            if tmp,
                B(len+1, num_cols) = 0;
            end
            % record the branch current
            LINELM(i, L_CUR_) = len + 1;
            
        case K_ %mutual inductance
            k = LINELM(i, K_VALUE_);
            ind1 = LINELM(i, K_ELEM1_);
            ind2 = LINELM(i, K_ELEM2_);
            m = k * sqrt( LINELM(ind1, L_VALUE_) * LINELM(ind2, L_VALUE_));
            cur1 = LINELM(ind1, L_CUR_);
            cur2 = LINELM(ind2, L_CUR_);
            C(cur1, cur2) = m;
            C(cur2, cur1) = m;
            
        case C_ %capacitance
            n1 = LINELM(i, C_N1_);
            n2 = LINELM(i, C_N2_);
            if n1 > 0,
                C(n1, n1) = C(n1, n1) + LINELM(i, C_VALUE_);
                if n2 > 0,
                    C(n1, n2) = C(n1, n2) - LINELM(i, C_VALUE_);
                end
            end
            if n2 > 0,
                C(n2, n2) = C(n2, n2) + LINELM(i, C_VALUE_);
                if n1 >0,
                    C(n2, n1) = C(n2, n1) - LINELM(i,C_VALUE_);
                end
            end
            
        case V_ % independent voltage source
            
            n1 = LINELM(i, V_N1_);
            n2 = LINELM(i, V_N2_);
            len = length(G);
            
            [tmp, num_cols] = size(B);
            add2b = 0;
            
            if (~grd_vsrc)
                add2b = 1;
            else
                if (abs(LINELM(i, V_VALUE_)) < 1e-20)
                    add2b = 0;
                else
                    add2b = 1;
                end
            end
            
            % add one column
            if (add2b)
                B(1, num_cols+1) = 0;
            end
            % need to add current as one variable
            % save the location when voltage equation is inserted.
            if n1 > 0 && n2 > 0,
                fprintf(' *Error stamp device Vsrc %d\n', i);
				stamp_ok = 0;
				continue
            end
            LINELM(i, BRANCH_CUR_IND_) = len + 1;
            if n1 > 0,
                G(len+1, n1) = -1;
                G(n1, len+1) = 1;
                % form b vector
                if add2b
                    B(len+1, num_cols+1) = -1;
                    port_ind = [port_ind len + 1];
                    port_name =[port_name; names(i,:)];
                else
                    if (size(B,2))
                        % only make the dimension equal
                        B(len+1, 1) = 0;
                    end
                end
                
            end
            if n2 > 0,
                G(len+1, n2) = 1;
                G(n2,   len+1) = -1;
                % form b vector
                if add2b
                    B(len+1, num_cols+1) = 1;
                    port_ind = [port_ind len + 1];
                    port_name =[port_name; names(i,:)];
                else
                    % only make the dimension equal
                    B(len+1, 1) = 0;
                end
            end
            % only makes the dimensions agree
            C(len+1, len+1) = 0;
            
        case I_
            [tmp, num_cols] = size(B);
            % add one column
            len = length(G);
            B(1, num_cols+1) = 0;
            num_cols = num_cols + 1;
            
            n1 = LINELM(i, I_N1_);
            n2 = LINELM(i, I_N2_);
            if n1 > 0 && n2 > 0,
                fprintf(' *Error stamp device Isrc %d\n', i);
				stamp_ok = 0;
				continue
            end
            
            if n1 >0,
                B(n1, num_cols) = 1;
                port_ind = [port_ind n1];
            end
            if n2 > 0,
                B(n2, num_cols) = -1;
                port_ind = [port_ind n2];
            end
            
        otherwise
            LINELM(i, TYPE_)
            fprintf(' *Error stamp unknown device %d\n', i);
			stamp_ok = 0;
			continue
    end
end

len = length(G);
if (length(B) < len)
    B(len, 1) = 0.0;
end
L = B;

end
%% end of function stamp_prima
