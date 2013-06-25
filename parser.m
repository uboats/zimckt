function [parse_ok,LINELEM,NLNELEM,INFO,NODES,NAMES,PRINTNV,PLOTNV,PLOTBI] = parser(ckt)


%cparse_init;
parser_init;

% EA's additions;
NLNELEM=[]; LINELEM=[];
INFO=[]; NODES=[];
NLNNAME=[]; LINNAME=[];
PRINTNV=[]; 
PLOTNV=[]; PLOTBI=[];

% load the input ckt
[parse_ok,ELEM,INFO,NODES,NAMES,PRINTNV,PLOTNV,PLOTBI_INIT] = loadckt(ckt);

if(parse_ok == 0)
    %fprintf(' Parser failed\n ---- Quit ----\n');
    return;
end

% analyze the elements
BOGUS_ = 0;
LINEAR_ = 1;
NONLINEAR_ = 2;
countLIN = 0;
countNLN = 0;
for i=1:size(ELEM,1),
    if (ELEM(i,TYPE_) == R_),
        countLIN = countLIN + 1;
        LINELEM(countLIN,TYPE_)      = R_;
        LINELEM(countLIN,R_VALUE_)   = ELEM(i,VALUE_);
        LINELEM(countLIN,R_N1_)      = ELEM(i,N1_);
        LINELEM(countLIN,R_N2_)      = ELEM(i,N2_);
        LINNAME(countLIN,:) = NAMES(i,:);
        ELEMNUM(i,1:2) = [countLIN,LINEAR_];
        
    elseif (ELEM(i,TYPE_) == C_),
        countLIN = countLIN + 1;
        LINELEM(countLIN,TYPE_)      = C_;
        LINELEM(countLIN,C_VALUE_)   = ELEM(i,VALUE_);
        LINELEM(countLIN,C_IC_)      = ELEM(i,IC_);
        LINELEM(countLIN,C_N1_)      = ELEM(i,N1_);
        LINELEM(countLIN,C_N2_)      = ELEM(i,N2_);
        LINNAME(countLIN,:) = NAMES(i,:);
        ELEMNUM(i,1:2) = [countLIN,LINEAR_];
        
    elseif (ELEM(i,TYPE_) == L_),
        countLIN = countLIN + 1;
        LINELEM(countLIN,TYPE_)      = L_;
        LINELEM(countLIN,L_VALUE_)   = ELEM(i,VALUE_);
        LINELEM(countLIN,L_IC_)      = ELEM(i,IC_);
        LINELEM(countLIN,L_N1_)      = ELEM(i,N1_);
        LINELEM(countLIN,L_N2_)      = ELEM(i,N2_);
        LINNAME(countLIN,:) = NAMES(i,:);
        ELEMNUM(i,1:2) = [countLIN,LINEAR_];
        
    elseif (ELEM(i,TYPE_) == V_),
        countLIN = countLIN + 1;
        LINELEM(countLIN,TYPE_)      = V_;
        LINELEM(countLIN,V_VALUE_)   = ELEM(i,VALUE_);
        LINELEM(countLIN,V_N1_)      = ELEM(i,N1_);
        LINELEM(countLIN,V_N2_)      = ELEM(i,N2_);
        LINELEM(countLIN,V_TYPE_)    = ELEM(i,DAC_);
        
        if (LINELEM(countLIN,V_TYPE_) == AC_)
            LINELEM(countLIN,V_PHASE_) = ELEM(i,PHASE_);
        elseif (LINELEM(countLIN,V_TYPE_) == PWL_)
            LINELEM(countLIN,V_POINTS_)  = ELEM(i,PWL_START_V_);
            for j = 1:2*ELEM(i,PWL_START_V_)
                LINELEM(countLIN,V_POINTS_ + j) = ELEM(i,PWL_START_V_ + j);
            end
        elseif (LINELEM(countLIN,V_TYPE_) == SIN_);
            LINELEM(countLIN,V_DCAMP_) = ELEM(i,V_DCAMP_);
            LINELEM(countLIN,V_ACAMP_) = ELEM(i,V_ACAMP_);
            LINELEM(countLIN,V_SINFREQ_) = ELEM(i,V_SINFREQ_);
            LINELEM(countLIN,V_SINPHASE_) = ELEM(i,V_SINPHASE_);
        elseif(LINELEM(countLIN,V_TYPE_) == PULSE_)
            %LINELEM(countLIN,I_TYPE_)
            LINELEM(countLIN,PLS_L_) = ELEM(i,PLS_L_);
            LINELEM(countLIN,PLS_H_) = ELEM(i,PLS_H_);
            LINELEM(countLIN,PLS_D_) = ELEM(i,PLS_D_);
            LINELEM(countLIN,PLS_R_) = ELEM(i,PLS_R_);
            LINELEM(countLIN,PLS_W_) = ELEM(i,PLS_W_);
            LINELEM(countLIN,PLS_F_) = ELEM(i,PLS_F_);
            LINELEM(countLIN,PLS_P_) = ELEM(i,PLS_P_);
        end
        
        LINNAME(countLIN,:) = NAMES(i,:);
        ELEMNUM(i,1:2) = [countLIN,LINEAR_];
        
        if(INFO(METHOD_) == DC_SWEEP)
            if(i == INFO(SWEEP_DEV) )
                INFO(SWEEP_ELEM) = countLIN;
            end
        end
        
    elseif (ELEM(i,TYPE_) == I_),
        countLIN = countLIN + 1;
        LINELEM(countLIN,TYPE_)      = I_;
        LINELEM(countLIN,I_VALUE_)   = ELEM(i,VALUE_);
        LINELEM(countLIN,I_N1_)      = ELEM(i,N1_);
        LINELEM(countLIN,I_N2_)      = ELEM(i,N2_);
        LINELEM(countLIN,I_TYPE_)    = ELEM(i,DAC_);
        
        if (LINELEM(countLIN,I_TYPE_) == AC_)
            LINELEM(countLIN,I_PHASE_) = ELEM(i,PHASE_);
        elseif (LINELEM(countLIN,I_TYPE_) == PWL_)
            LINELEM(countLIN,I_POINTS_)  = ELEM(i,PWL_START_I_);
            for j = 1:2*ELEM(i,PWL_START_I_)
                LINELEM(countLIN,I_POINTS_ + j) = ELEM(i,PWL_START_I_ + j);
            end
        elseif (LINELEM(countLIN,I_TYPE_) == SIN_);
            LINELEM(countLIN,I_DCAMP_) = ELEM(i,I_DCAMP_);
            LINELEM(countLIN,I_ACAMP_) = ELEM(i,I_ACAMP_);
            LINELEM(countLIN,I_SINFREQ_) = ELEM(i,I_SINFREQ_);
            LINELEM(countLIN,I_SINPHASE_) = ELEM(i,I_SINPHASE_);
        elseif(LINELEM(countLIN,I_TYPE_) == PULSE_)
            %LINELEM(countLIN,I_TYPE_)
            LINELEM(countLIN,PLS_L_) = ELEM(i,PLS_L_);
            LINELEM(countLIN,PLS_H_) = ELEM(i,PLS_H_);
            LINELEM(countLIN,PLS_D_) = ELEM(i,PLS_D_);
            LINELEM(countLIN,PLS_R_) = ELEM(i,PLS_R_);
            LINELEM(countLIN,PLS_W_) = ELEM(i,PLS_W_);
            LINELEM(countLIN,PLS_F_) = ELEM(i,PLS_F_);
            LINELEM(countLIN,PLS_P_) = ELEM(i,PLS_P_);
        end
        
        LINNAME(countLIN,:) = NAMES(i,:);
        ELEMNUM(i,1:2) = [countLIN,LINEAR_];
        
        if(INFO(METHOD_) == DC_SWEEP)
            if(i == INFO(SWEEP_DEV) )
                INFO(SWEEP_ELEM) = countLIN;
            end
        end
        
    elseif (ELEM(i,TYPE_) == M_),
        countNLN = countNLN + 1;
        NLNELEM(countNLN,TYPE_)      = M_;
        NLNELEM(countNLN,M_W_)       = ELEM(i,MOS_W_);
        NLNELEM(countNLN,M_L_)       = ELEM(i,MOS_L_);
        NLNELEM(countNLN,M_ND_)      = ELEM(i,MOS_ND_);
        NLNELEM(countNLN,M_NG_)      = ELEM(i,MOS_NG_);
        NLNELEM(countNLN,M_NS_)      = ELEM(i,MOS_NS_);
        NLNELEM(countNLN,M_TYPE_)    = ELEM(i,MOS_TYPE_);
        NLNELEM(countNLN,M_MID_)     = ELEM(i,MOS_MID_);
        NLNNAME(countNLN,:) = NAMES(i,:);
        ELEMNUM(i,1:2) = [countNLN,NONLINEAR_];
        
    elseif (ELEM(i,TYPE_) == '.'),
        MDLELEM(ELEM(i,MOD_ID_),MOD_ID_)       = ELEM(i,MOD_ID_);
        MDLELEM(ELEM(i,MOD_ID_),MOD_VT_)       = ELEM(i,MOD_VT_);
        MDLELEM(ELEM(i,MOD_ID_),MOD_MU_)       = ELEM(i,MOD_MU_);
        MDLELEM(ELEM(i,MOD_ID_),MOD_COX_)      = ELEM(i,MOD_COX_);
        MDLELEM(ELEM(i,MOD_ID_),MOD_LAMBDA_)   = ELEM(i,MOD_LAMBDA_);
        MDLELEM(ELEM(i,MOD_ID_),MOD_CJ0_)      = ELEM(i,MOD_CJ0_);
        ELEMNUM(i,1:2) = [0,BOGUS_];
    end
        
end

% append the model parameters to the MOSFET line in NLNELEM
for i = 1:size(NLNELEM,1),
    index = NLNELEM(i,M_MID_);
    for j = M_MID_:(M_MID_+size(MDLELEM,2)-MOD_ID_-1),
        NLNELEM(i,j) = MDLELEM(index,j-M_MID_+MOD_ID_+1);
    end
end

for i = 1:size(PLOTBI_INIT,1),
    PLOTBI(i,1:2) = ELEMNUM(PLOTBI_INIT(i),1:2);
end

% correct
for i = 1:size(PLOTBI_INIT,1),
    PLOTBI(i,1) = ELEMNUM(PLOTBI_INIT(i),1);
    PLOTBI(i,2) = NLNELEM(ELEMNUM(PLOTBI_INIT(i),1),M_TYPE_);
    PLOTBI(i,3) = ELEM(PLOTBI_INIT(i),MOS_ND_);
    PLOTBI(i,4) = ELEM(PLOTBI_INIT(i),MOS_NG_);
    PLOTBI(i,5) = ELEM(PLOTBI_INIT(i),MOS_NS_);
    PLOTBI(i,6) = PLOTBI_INIT(i);
end


%PLOTBI
