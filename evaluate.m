function evaluate(numNodes,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% evaluate: evaluate device and stamp device values
%%           all evaluation is to solve delta X in N-R solver
%%
%% - numNodes: circuit size
%% - setup   : if 1, create system matrices
%%             otherwise, stamp elements to matrices
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global LINELM NLINELM BE_ TR_ tt
global F G C_idx C_num num_C dim

% clear G matrix and Capacitor idx
G = spalloc(numNodes,numNodes,0);
F = zeros(dim,1);
% number of capacitors
C_num=0;
C_idx=1;

tt = BE_;

if (setup==1)
    dim=numNodes;
end
%   port_ind = [];

%   linear elements
numLINE = size(LINELM,1);
ldev(numLINE,setup);

%   Nonlinear elements
numNLINE = size(NLINELM,1);
nldev(numNLINE,setup);

if(setup==1)
    num_C = C_num;
end

end
%% end of function evaluate


function ldev(numLINE,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% ldev: stamp linear devices
%%
%% - numLINE: number of linear devices
%% - setup  : if 1, create system matrices
%%            otherwise, stamp elements to matrices
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global LINELM TYPE_ C_ L_ R_ V_ I_ PI 
global R_N1_ R_N2_ C_N1_ C_N2_ L_N1_ L_N2_ L_IC_
global R_VALUE_ C_VALUE_ L_VALUE_
global V_N1_ V_N2_ V_VALUE_ V_TYPE_ V_POINTS_ 
global I_N1_ I_N2_ I_VALUE_ I_TYPE_ I_POINTS_
global PWL_ SIN_ L_CUR_ T_ Delta_T BE_ TR_ tt
global F G TRAN_ DC_ S_type C_idx C_num dim
global AC_ freq_t vin_amp cur_amp DA_
global V_SINFREQ_ V_SINPHASE_ V_DCAMP_ V_ACAMP_
global I_SINFREQ_ I_SINPHASE_ I_DCAMP_ I_ACAMP_
global PULSE_ PLS_L_ PLS_H_ PLS_D_ PLS_R_ PLS_W_ PLS_F_ PLS_P_

for i = 1:numLINE
    switch LINELM(i, TYPE_)
        case R_ %resistance
            if(setup==0)
                n1 = LINELM(i, R_N1_);
                n2 = LINELM(i, R_N2_);
                g = 1/LINELM(i, R_VALUE_);
                Stamp_Res(n1,n2,g);
            end
            
        case L_ %inductance
            n1 = LINELM(i, L_N1_);
            n2 = LINELM(i, L_N2_);
            % need to add inductance current as one variable
            len = length(G);
            if(setup==1)
                dim=dim+1;
                LINELM(i,L_IC_) = dim;
            else
                if(S_type == AC_)
                    val = LINELM(i, L_VALUE_);
                    omiga = 2*PI*freq_t;
                    val = omiga*val*1i;
                    val = 1/val;
                else
                    val = Delta_T/(LINELM(i, L_VALUE_));
                end
                L_CUR_idx = LINELM(i, L_IC_);
                Stamp_Ind(n1,n2,val,L_CUR_idx)
                LINELM(i, L_CUR_) = len + 1;
            end
            
        case C_ %capacitance
            if(setup==1)
                C_num=C_num+1;
            else
                if(S_type == TRAN_) && (setup==0)
                    n1 = LINELM(i, C_N1_);
                    n2 = LINELM(i, C_N2_);
                    
                    if(tt == BE_)
                        val = (LINELM(i, C_VALUE_))/Delta_T;
                    elseif(tt == TR_)
                        val = 2*(LINELM(i, C_VALUE_))/Delta_T;
                    end
                    Stamp_Cap(n1,n2,val,C_idx);
                    C_idx=C_idx+1;
                elseif(S_type == AC_)
                    n1 = LINELM(i, C_N1_);
                    n2 = LINELM(i, C_N2_);
                    
                    val = LINELM(i, C_VALUE_);
                    omiga = 2*PI*freq_t;
                    c_val = 1/(omiga*val*1i);
                    c_val = 1/c_val;
                    Stamp_Cap(n1,n2,c_val,C_idx);
                    C_idx=C_idx+1;
                end
            end
            
        case I_
            if(setup==0)
                n1 = LINELM(i, I_N1_);
                n2 = LINELM(i, I_N2_);
                
                if((S_type == TRAN_) && (LINELM(i,I_TYPE_) == PWL_))
                    pwl_times=[];
                    pwl_vals=[];
                    slope=0;
                    numpts=LINELM(i,I_POINTS_);
                    starti=LINELM(i,I_VALUE_);
                    pwl_times = zeros(numpts,1);
                    pwl_vals = zeros(numpts,1);
                    pwl_vals(1) = starti;
                    for j = 1: numpts,
                        pwl_times(j) = LINELM(i,I_POINTS_ + 2*j-1) ;
                        pwl_vals(j) = LINELM(i,I_POINTS_ + 2*j);
                    end
                    find=0;
                    j=1;
                    while (T_ >= pwl_times(j)) && (T_<=pwl_times(numpts) && (j<numpts))
                        if(T_ == pwl_times(j)),
                            cur = pwl_vals(j);
                            find=1;
                            break
                        end
                        j=j+1;
                    end
                    if (j<=numpts) && (find==0) && (T_<= pwl_times(numpts))
                        slope = (pwl_vals(j) - pwl_vals(j-1))/(pwl_times(j) - pwl_times(j-1));
                        cur = pwl_vals(j-1) + slope * (T_- pwl_times(j-1));
                    elseif (find==0)
                        cur = pwl_vals(numpts);
                    end
                    % SIN sources
                elseif((S_type == TRAN_) && (LINELM(i,I_TYPE_) == SIN_))
                    cur_dc = LINELM(i,I_DCAMP_);
                    cur_ac = LINELM(i,I_ACAMP_);
                    freq = LINELM(i,I_SINFREQ_);
                    phase = LINELM(i,I_SINPHASE_);
                    cur = cur_dc + cur_ac * sin(freq*T_ - phase);
                
				elseif((S_type == TRAN_) && (LINELM(i,I_TYPE_) == PULSE_))
					i_l = LINELM(i,PLS_L_);
					i_h = LINELM(i,PLS_H_);
					t_d = LINELM(i,PLS_D_);
					t_r = LINELM(i,PLS_R_);
					t_w = LINELM(i,PLS_W_);
					t_f = LINELM(i,PLS_F_);
					t_p = LINELM(i,PLS_P_);

					t_ = T_ - t_d;
					n_p = floor(t_/t_p);
					t_ = t_ - n_p * t_p;

					if(t_ <= t_r)
						cur = (i_h-i_l)/t_r*t_;
					elseif(t_ <= t_r+t_w)
						cur = i_h;
					elseif(t_ <= t_r+t_w+t_f)
						cur = i_h - (i_h-i_l)/t_f*(t_-t_r-t_w);
					else
						cur = i_l;
					end

                    % DC sources
                elseif((S_type == TRAN_) && (LINELM(i,I_TYPE_) == DC_ || LINELM(i,I_TYPE_) == DA_))
                    cur = LINELM(i,I_VALUE_);
                elseif(S_type == DC_ || LINELM(i,I_TYPE_) == DA_)
                    cur = LINELM(i,I_VALUE_);
                elseif(S_type == AC_ && LINELM(i,I_TYPE_) == DC_)
                    cur = 0;
                elseif(S_type == AC_ && LINELM(i,I_TYPE_) == AC_)
                    cur = LINELM(i,I_VALUE_);
                    cur_amp = max(cur_amp, i);
                elseif(S_type == AC_ && LINELM(i,I_TYPE_) == DA_)
                    cur = LINELM(i,I_VALUE_);
                    cur_amp = max(cur_amp, v);
                end
                
                if(n1>0)
                    F(n1) = F(n1)-cur;
                end
                if(n2>0)
                    F(n2) = F(n2)+cur;
                end
            end
            
        case V_ % independent voltage source
            
            n1 = LINELM(i, V_N1_);
            n2 = LINELM(i, V_N2_);
            len = length(G);
            if(setup==1)
                dim=dim+1;
            else
                if(n1 > 0)
                    G(n1, len+1) = -1;
                    G(len+1, n1) = 1;
                    F(n1)=F(n1);%-X(len+1);
                    F(len+1)=F(len+1);%-X(n1);
                end
                if(n2 > 0)
                    G(n2, len+1) = 1;
                    G(len+1, n2) = -1;
                    F(n2)=F(n2);%+X(len+1);
                    F(len+1)=F(len+1);%+X(n2);
                end
                
                if((S_type == TRAN_) && (LINELM(i,V_TYPE_) == PWL_))
                    pwl_times=[];
                    pwl_vals=[];
                    slope=0;
                    numpts=LINELM(i,V_POINTS_);
                    startv=LINELM(i,V_VALUE_);
                    pwl_times = zeros( numpts,1);
                    pwl_vals = zeros( numpts,1);
                    pwl_vals(1) = startv;
                    for j = 1: numpts,
                        pwl_times(j) = LINELM(i,V_POINTS_ + 2*j-1) ;
                        pwl_vals(j) = LINELM(i,V_POINTS_ + 2*j);
                    end
                    find=0;
                    j=1;
                    while (T_ >= pwl_times(j)) && (T_<=pwl_times(numpts) && (j<numpts))
                        if(T_ == pwl_times(j)),
                            v = pwl_vals(j);
                            F(len+1) = v;
                            find=1;
                            break
                        end
                        j=j+1;
                    end
                    if (j<=numpts) && (find==0) && (T_<= pwl_times(numpts))
                        slope = (pwl_vals(j) - pwl_vals(j-1))/(pwl_times(j) - pwl_times(j-1));
                        v = pwl_vals(j-1) + slope * (T_- pwl_times(j-1));
                        F(len+1) = v;
                        if(n2>0)
                            % F(len+1)=F(len+1)+ X(n2);
                            find=1;
                        end
                        find=1;
                    elseif (find==0)
                        v = pwl_vals(numpts);
                        F(len+1) = v;
                        if(n2>0)
                            % F(len+1)=F(len+1)+ X(n2);
                            find=1;
                        end
                    end
                % SIN sources
                elseif((S_type == TRAN_) && (LINELM(i,V_TYPE_) == SIN_))
                    v_dc = LINELM(i,V_DCAMP_);
                    v_ac = LINELM(i,V_ACAMP_);
                    freq = LINELM(i,V_SINFREQ_);
                    phase = LINELM(i,V_SINPHASE_);
                    v = v_dc + v_ac * sin(freq*T_ - phase);
                    F(len+1) = v;
                elseif((S_type == TRAN_) && (LINELM(i,V_TYPE_) == PULSE_))
					v_l = LINELM(i,PLS_L_);
					v_h = LINELM(i,PLS_H_);
					t_d = LINELM(i,PLS_D_);
					t_r = LINELM(i,PLS_R_);
					t_w = LINELM(i,PLS_W_);
					t_f = LINELM(i,PLS_F_);
					t_p = LINELM(i,PLS_P_);

					t_ = T_ - t_d;
					n_p = floor(t_/t_p);
					t_ = t_ - n_p * t_p;

					if(t_ <= t_r)
						v = (v_h-v_l)/t_r*t_;
					elseif(t_ <= t_r+t_w)
						v = v_h;
					elseif(t_ <= t_r+t_w+t_f)
						v = v_h - (v_h-v_l)/t_f*(t_-t_r-t_w);
					else
						v = v_l;
					end

					F(len+1) = v;
                    
                % DC sources
                elseif((S_type == TRAN_) && (LINELM(i,V_TYPE_) == DC_ || LINELM(i,V_TYPE_) == DA_))
                    v = LINELM(i,V_VALUE_);
                    F(len+1) = v;
                elseif(S_type == DC_ || LINELM(i,V_TYPE_) == DA_)
                    v = LINELM(i,V_VALUE_);
                    F(len+1) = v;
                elseif(S_type == AC_ && LINELM(i,V_TYPE_) == DC_)
                    F(len+1) = 0;
                elseif(S_type == AC_ && LINELM(i,V_TYPE_) == AC_)
                    v = LINELM(i,V_VALUE_);
                    F(len+1) = v;
                    vin_amp = max(vin_amp, v);
                elseif(S_type == AC_ && LINELM(i,V_TYPE_) == DA_)
                    v = LINELM(i,V_VALUE_);
                    F(len+1) = v;
                    vin_amp = max(vin_amp, v);
                else
                    v = LINELM(i,V_VALUE_);
                    F(len+1) = v;
                end
                
            end
    end
end
end
%% end of function ldev


function nldev(numNLINE,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% evaluate: stamp nonlinear device
%%
%% - numNLINE: number of nonlinear devices
%% - setup   : if 1, create system matrices
%%             otherwise, stamp elements to matrices
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global NLINELM TYPE_ Delta_T BE_ TR_ tt M_TYPE_
global M_ND_ M_NG_ M_NS_ M_W_ M_L_ M_VT_ M_MU_ M_COX_ M_LAMBDA_ M_CJ0_ M_
global TRAN_ S_type AC_ Res_bi_mos C_idx C_num

for i = 1:numNLINE,
    switch NLINELM(i, TYPE_)
        case M_    %mosfets
            if(setup==1)
                C_num=C_num+4;
            end
            if(setup==0)
                nd = NLINELM(i, M_ND_);
                ng = NLINELM(i, M_NG_);
                ns = NLINELM(i, M_NS_);
                W=NLINELM(i, M_W_);
                L=NLINELM(i, M_L_);
                Vt=NLINELM(i,M_VT_);
                Mu=NLINELM(i,M_MU_);
                Cox=NLINELM(i,M_COX_);
                Lambda=NLINELM(i,M_LAMBDA_);
                CJ0=NLINELM(i,M_CJ0_);
                type=NLINELM(i, M_TYPE_);
                ids=Stamp_Mos(nd,ng,ns,Vt,Mu,Cox,W,L,Lambda,type);
                Res_bi_mos(i)=ids;
                if(S_type == TRAN_)
                    % stamp 4 capacitors
                    Cgs=0.5*Cox*W*L;
                    Cgd=Cgs;
                    Cd=CJ0;
                    Cs=Cd;
                    ratio = 1;
                    if(tt==BE_)
                        ratio=1;
                    elseif(tt==TR_)
                        ratio=2;
                    end
                    % Cgs
                    val = ratio*(Cgs)/Delta_T;
                    Stamp_Cap(ng,ns,val,C_idx);
                    C_idx=C_idx+1;
                    % Cgd
                    val = ratio*(Cgd)/Delta_T;
                    Stamp_Cap(ng,nd,val,C_idx);
                    C_idx=C_idx+1;
                    % Cd: grounded cap
                    val = ratio*(Cs)/Delta_T;
                    Stamp_Cap(ns,-1,val,C_idx);
                    C_idx=C_idx+1;
                    % Cs:grounded cap
                    val = ratio*(Cd)/Delta_T;
                    Stamp_Cap(nd,-1,val,C_idx)
                    C_idx=C_idx+1;
                    
                elseif(S_type == AC_)
                    
                end
            end
    end
end
end
%% end of function nldev


function Stamp_Res(n1,n2,g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Stamp_Res: stamp resistor
%%
%% - n1: node 1
%% - n2: node 2
%% - g : value
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global G
% Put your codes here
if (n1>0) && (n2>0)
    G(n1,n1)=G(n1,n1)+g;
    G(n1,n2)=G(n1,n2)-g;
    G(n2,n1)=G(n2,n1)-g;
    G(n2,n2)=G(n2,n2)+g;
elseif (n1>0) && (n2<0)
    G(n1,n1)=G(n1,n1)+g;
elseif (n1<0) && (n2>0)
    G(n2,n2)=G(n2,n2)+g;
elseif (n1<0) && (n2<0)
    disp(' Error: unknow res connection')
    return;
else
    disp(' Error: unknow nodes')
    return;
end
end
%% end of function Stamp_Res


function Stamp_Cap(n1,n2,val,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Stamp_Cap: stamp capacitor (BE approximation)
%%
%% - n1  : node 1
%% - n2  : node 2
%% - val : value
%% - i   : device index
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TR_ BE_ G F X_pre_t I_pre I_eqp tt tr_n
global S_type TRAN_
%tt = BE_;

if (n1>0 && n2>0)
    G(n1, n1) = G(n1, n1) + val;
    G(n1, n2) = G(n1, n2) - val;
    G(n2, n1) = G(n2, n1) - val;
    G(n2, n2) = G(n2, n2) + val;
elseif n1>0
    G(n1, n1) = G(n1, n1) + val;
elseif n2>0
    G(n2, n2) = G(n2, n2) + val;
end

if(S_type == TRAN_)
    if(tt == TR_)
        if (n1>0 && n2>0)
            I_eq = val * (X_pre_t(n1) - X_pre_t(n2)) + I_pre(i);
            if(tr_n==1)
                I_eqp(i)=I_eq;
            end
            F(n1) = F(n1) + I_eq;
            F(n2) = F(n2) - I_eq;
            I_pre(i) = val * (X_pre_t(n1) - X_pre_t(n2)) - I_eqp(i);
            I_eqp(i) = I_eq;
        elseif n1>0
            I_eq = val * X_pre_t(n1) + I_pre(i);
            if(tr_n==1)
                I_eqp(i)=I_eq;
            end
            F(n1) = F(n1) + I_eq;
            I_pre(i) = val * X_pre_t(n1) - I_eqp(i);
            I_eqp(i) = I_eq;
        elseif n2>0
            I_eq = val * X_pre_t(n2) + I_pre(i);
            if(tr_n==1)
                I_eqp(i)=I_eq;
            end
            F(n2) = F(n2) - I_eq;
        end
    elseif(tt == BE_)
        if (n1>0 && n2>0)
            I_eq = val * (X_pre_t(n1) - X_pre_t(n2));
            F(n1) = F(n1) + I_eq;
            F(n2) = F(n2) - I_eq;
        elseif n1>0
            I_eq = val * X_pre_t(n1);
            F(n1) = F(n1) + I_eq;
        elseif n2>0
            I_eq = val * X_pre_t(n2);
            F(n2) = F(n2) - I_eq;
        end
    end
end
end
%% end of function Stamp_Cap


function Stamp_Ind(n1,n2,val,L_CUR_idx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Stamp_Ind: stamp inductor (BE approximation)
%%
%% - n1       : node 1
%% - n2       : node 2
%% - val      : value
%% - L_CUR_idx: pseudo node index
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global G F X_pre_t
global  S_type TRAN_ DC_ AC_
% Put your codes here
%gshort = 1e+6;

if(S_type == DC_)
    if (n1>0 && n2>0)
        G(n1, L_CUR_idx) = 1;
        G(n2, L_CUR_idx) = -1;
        G(L_CUR_idx, n1) = -1;
        G(L_CUR_idx, n2) = 1;
    elseif n1>0
        G(n1, L_CUR_idx) = 1;
        G(L_CUR_idx, n1) = -1;
    elseif n2>0
        G(n2, L_CUR_idx) = -1;
        G(L_CUR_idx, n2) = 1;
    end
    
elseif(S_type == TRAN_)
    if (n1>0 && n2>0)
        G(n1, L_CUR_idx) = 1;
        G(n2, L_CUR_idx) = -1;
        G(L_CUR_idx, n1) = -1;
        G(L_CUR_idx, n2) = 1;
        G(L_CUR_idx, L_CUR_idx) = 1/val;
        
        F(L_CUR_idx) = X_pre_t(L_CUR_idx)/val;
    elseif n1>0
        G(n1, L_CUR_idx) = 1;
        G(L_CUR_idx, n1) = -1;
        G(L_CUR_idx, L_CUR_idx) = 1/val;
        
        F(L_CUR_idx) = X_pre_t(L_CUR_idx)/val;
    elseif n2>0
        G(n2, L_CUR_idx) = -1;
        G(L_CUR_idx, n2) = 1;
        G(L_CUR_idx, L_CUR_idx) = 1/val;
        
        F(L_CUR_idx) = X_pre_t(L_CUR_idx)/val;
    end
    
elseif(S_type == AC_)
    if (n1>0 && n2>0)
        G(n1, L_CUR_idx) = 1;
        G(n2, L_CUR_idx) = -1;
        G(L_CUR_idx, n1) = -1;
        G(L_CUR_idx, n2) = 1;
        G(L_CUR_idx, L_CUR_idx) = 1/val;
        
        %   %F(L_CUR_idx) = X_pre_t(L_CUR_idx)/val;
    elseif n1>0
        G(n1, L_CUR_idx) = 1;
        G(L_CUR_idx, n1) = -1;
        G(L_CUR_idx, L_CUR_idx) = 1/val;
        
        %   %F(L_CUR_idx) = X_pre_t(L_CUR_idx)/val;
    elseif n2>0
        G(n2, L_CUR_idx) = -1;
        G(L_CUR_idx, n2) = 1;
        G(L_CUR_idx, L_CUR_idx) = 1/val;
        
        %   %F(L_CUR_idx) = X_pre_t(L_CUR_idx)/val;
    end
end
end
%% end of function Stamp_Ind


% stamp mosfet
function [ids]=Stamp_Mos(nd,ng,ns,Vt,Mu,Cox,W,L,Lambda,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Stamp_Mos: stamp mosfet
%%
%% - nd    : drain
%% - ng    : gate
%% - ns    : source
%% - Vt    : threshold
%% - Mu    : transconductance coefficient
%% - Cox   : coupling capacitance
%% - W     : channel width
%% - L     : channel length
%% - Lambda: channel-length modulation
%% - type  : 1 nmos, 0 pmos
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOS stamping function:
% type =1 for NMOS;type=0 for PMOS
% solver for delta X
% 6 situations are considered here:
%    1    2   3   4    5     6
% Ng >0   >0  >0  <0   <0    OTHERS
% Ns >0   >0  <0  >0   >0    OTHERS
% Nd >0   <0  >0  >0   <0    OTHERS
global G F X
if  nd>0 && ng>0 && ns>0
    Vgs=X(ng)-X(ns);
    Vds=X(nd)-X(ns);
elseif nd<0 && ng>0 && ns>0
    Vgs=X(ng)-X(ns);
    Vds=-X(ns);
elseif nd>0 && ng>0 && ns<0
    Vgs=X(ng);
    Vds=X(nd);
elseif nd>0 && ng<0 && ns>0
    Vgs=-X(ns);
    Vds=X(nd)-X(ns);
elseif nd<0 && ng<0 && ns>0
    Vgs=-X(ns);
    Vds=-X(ns);
else
    disp(' Error: unknown MOSFET DC condition. Quit!')
    return;
end
%PMOS fac=-1; NMOS fac=1;
fac = (-1)^(1-type);
Vgs = fac*(Vgs);
Vds = fac*(Vds);
gm = MOSVI_Gm(Vgs,Vds,Vt,Mu,Cox,W,L,Lambda,type);
gds = MOSVI_Gds(Vgs,Vds,Vt,Mu,Cox,W,L,Lambda,type);
ids = MOSVI_Ids(Vgs,Vds,Vt,Mu,Cox,W,L,Lambda,type);
ids = fac*ids;
%Stamp voltage-control current source Vgs,gm
if nd>0 && ns>0 && ng>0,
    G(nd, ng) = G(nd, ng) + gm;
    G(nd, ns) = G(nd, ns) - gm;
    G(ns, ng) = G(ns, ng) - gm;
    G(ns, ns) = G(ns, ns) + gm;
    %         %F(nd)=F(nd)-((X(ng)-X(ns))*gm+gds*(X(nd)-X(ns)));
    %         %F(ns)=F(ns)+((X(ng)-X(ns))*gm+gds*(X(nd)-X(ns)));
    F(nd)=F(nd)+((X(ng)-X(ns))*gm+gds*(X(nd)-X(ns)));
    F(ns)=F(ns)-((X(ng)-X(ns))*gm+gds*(X(nd)-X(ns)));
    F(nd)=F(nd)-ids;
    F(ns)=F(ns)+ids;
    %         %F(nd)=F(nd)+ids;
    %         %F(ns)=F(ns)-ids;
elseif nd<0 && ns>0 && ng>0,
    G(ns, ns) = G(ns, ns) + gm;
    G(ns, ng) = G(ns, ng) - gm;
    %         %F(ns)=F(ns)+((X(ng)-X(ns))*gm+gds*(-X(ns)));
    F(ns)=F(ns)-((X(ng)-X(ns))*gm+gds*(-X(ns)));
    
    F(ns)=F(ns)+ids;
    %        %F(ns)=F(ns)-ids;
elseif nd>0 && ns<0 && ng>0,
    G(nd, ng) = G(nd, ng) + gm;
    %        %F(nd)=F(nd)-((X(ng))*gm+gds*(X(nd)));
    F(nd)=F(nd)+((X(ng))*gm+gds*(X(nd)));
    F(nd)=F(nd)-ids;
    %        %F(nd)=F(nd)+ids;
    
elseif nd>0 && ns>0 && ng<0,
    G(nd, ns) = G(nd, ns) - gm;
    G(ns, ns) = G(ns, ns) + gm;
    %        %F(nd)=F(nd)-((-X(ns))*gm+gds*(X(nd)-X(ns)));
    %        %F(ns)=F(ns)+((-X(ns))*gm+gds*(X(nd)-X(ns)));
    F(nd)=F(nd)+((-X(ns))*gm+gds*(X(nd)-X(ns)));
    F(ns)=F(ns)-((-X(ns))*gm+gds*(X(nd)-X(ns)));
    F(nd)=F(nd)-ids;
    F(ns)=F(ns)+ids;
    %        %F(nd)=F(nd)+ids;
    %        %F(ns)=F(ns)-ids;
elseif nd<0 && ns>0 && ng<0,
    G(ns, ns) = G(ns, ns) + gm;
    %        %F(ns)=F(ns)+((-X(ns))*gm+gds*(-X(ns)));
    F(ns)=F(ns)-((-X(ns))*gm+gds*(-X(ns)));
    
    F(ns)=F(ns)+ids;
    %        %F(ns)=F(ns)-ids;
end
%Stamp Gds
if nd>0 && ns>0,
    G(nd, nd) = G(nd, nd) + gds;
    G(nd, ns) = G(nd, ns) - gds;
    G(ns, nd) = G(ns, nd) - gds;
    G(ns, ns) = G(ns, ns) + gds;
elseif nd>0,
    G(nd, nd) = G(nd, nd) + gds;
elseif ns>0,
    G(ns, ns) = G(ns, ns) + gds;
end
end
%% end of function Stamp_Mos


% Mosfet function: Ids vs. Vds for NMOS and Isd vs. Vsd
function Ids= MOSVI_Ids(Vgs,Vds,Vt,Mu,Cox,W,L,Lambda,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% MOSVI_Ids: evaluate Ids
%%
%% - Vgs   : gate-2-source voltage
%% - Vds   : drain-2-source voltage
%% - Vt    : threshold
%% - Mu    : transconductance coefficient
%% - Cox   : coupling capacitance
%% - W     : channel width
%% - L     : channel length
%% - Lambda: channel-length modulation
%% - type  : 1 nmos, 0 pmos
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if type==1
    %NMOS
    if (Vgs>Vt) && (Vds>(Vgs-Vt)),
        Ids=0.5*Mu*Cox*W/L*(Vgs-Vt)^2*(1+Lambda*Vds);
    elseif (Vgs>Vt) && (Vds<=Vgs-Vt),
        Ids=Mu*Cox*W/L*((Vgs-Vt)*Vds-0.5*Vds^2)*(1+Lambda*Vds);
    elseif Vgs<=Vt,
        Ids=0;
    end
elseif type==0
    %PMOS
    Vsg=Vgs;
    Vsd=Vds;
    if (Vsg>-Vt) && (Vsd>(Vsg+Vt)),
        Isd=0.5*Mu*Cox*W/L*(Vsg+Vt)^2*(1+Lambda*Vds);
    elseif (Vsg>-Vt) && (Vsd<=Vsg+Vt),
        Isd=Mu*Cox*W/L*((Vsg+Vt)*Vsd-0.5*Vsd^2)*(1+Lambda*Vds);
    elseif Vsg<=-Vt,
        Isd=0;
    end
    Ids=Isd;
end
end
%% end of function MOSVI_Ids


% Mosfet function: Gds vs. Vds for NMOS and Gsd vs. dVsd
function Gds = MOSVI_Gds(Vgs,Vds,Vt,Mu,Cox,W,L,Lambda,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% MOSVI_Gds: evaluate equivalent conductance
%%
%% - Vgs   : gate-2-source voltage
%% - Vds   : drain-2-source voltage
%% - Vt    : threshold
%% - Mu    : transconductance coefficient
%% - Cox   : coupling capacitance
%% - W     : channel width
%% - L     : channel length
%% - Lambda: channel-length modulation
%% - type  : 1 nmos, 0 pmos
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if type==1
    %NMOS
    if (Vgs>Vt) && (Vds>(Vgs-Vt)),
        Gds=0.5*Mu*Cox*W/L*(Vgs-Vt)^2*(Lambda);
    elseif (Vgs>Vt) && (Vds<=Vgs-Vt),
        Gds=Mu*Cox*W/L*(((Vgs-Vt)*Vds-0.5*Vds^2)*(Lambda)+...
            ((Vgs-Vt)-Vds)*(1+Lambda*Vds));
    else
        Gds=0;
    end
elseif type==0
    %PMOS
    Vsg=Vgs;
    Vsd=Vds;
    if (Vsg>-Vt) && (Vsd>(Vsg+Vt)),
        Gds=0.5*Mu*Cox*W/L*(Vsg+Vt)^2*(Lambda);
    elseif (Vsg>-Vt) && (Vsd<=Vsg+Vt),
        Gds=Mu*Cox*W/L*(((Vsg+Vt)*Vsd-0.5*Vsd^2)*(Lambda)+...
            ((Vsg+Vt)-Vsd)*(1+Lambda*Vds));
    else
        Gds=0;
    end
end
end
%% end of function MOSVI_Gds


% Mosfet function: Ids vs. Vds for NMOS and Isd vs. Vsd
function Gm = MOSVI_Gm(Vgs,Vds,Vt,Mu,Cox,W,L,Lambda,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% MOSVI_Gm: evaluate trans-conductance
%%
%% - Vgs   : gate-2-source voltage
%% - Vds   : drain-2-source voltage
%% - Vt    : threshold
%% - Mu    : transconductance coefficient
%% - Cox   : coupling capacitance
%% - W     : channel width
%% - L     : channel length
%% - Lambda: channel-length modulation
%% - type  : 1 nmos, 0 pmos
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if type==1
    %NMOS
    if (Vgs>Vt) && (Vds>(Vgs-Vt)),
        Gm=Mu*Cox*W/L*(Vgs-Vt)*(1+Lambda*Vds);
    elseif (Vgs>Vt) && (Vds<=Vgs-Vt),
        Gm=Mu*Cox*W/L*(Vds)*(1+Lambda*Vds);
    else
        Gm=0;
    end
elseif type==0
    %PMOS
    Vsg=Vgs;
    Vsd=Vds;
    if (Vsg>-Vt) && (Vsd>(Vsg+Vt)),
        Isd=Mu*Cox*W/L*(Vsg+Vt)*(1+Lambda*Vds);
    elseif (Vsg>-Vt) && (Vsd<=Vsg+Vt),
        Isd=Mu*Cox*W/L*(Vsd)*(1+Lambda*Vds);
    else
        Isd=0;
    end
    Gm=Isd;
end
end
%% end of function MOSVI_Gm