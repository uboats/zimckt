%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% parser_init - contains short cut constants which should be used
%%               for accessing the data contained in the output of
%%               parser (see: help parser)
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global R_ M_ N_ C_ L_ V_ I_
global TYPE_
global R_N1_ R_N2_ R_VALUE_
global C_N1_ C_N2_ C_VALUE_ C_IC_
global L_N1_ L_N2_ L_VALUE_ L_IC_
global DC_ AC_ PWL_ SIN_ TRAN_ DA_
global V_N1_ V_N2_ V_VALUE_ V_TYPE_ V_POINTS_ V_PHASE_ V_DCAMP_ V_ACAMP_ V_SINFREQ_ V_SINPHASE_
global I_N1_ I_N2_ I_VALUE_ I_TYPE_ I_POINTS_ I_PHASE_ I_DCAMP_ I_ACAMP_ I_SINFREQ_ I_SINPHASE_
% global E_VALUE_ E_N1_ E_N2_ E_CN1_ E_CN2_
% global F_VALUE_ F_N1_ F_N2_ F_SOURCE_IND_ F_SOURCE_MAT_
% global G_VALUE_ G_N1_ G_N2_ G_CN1_ G_CN2_
% global H_VALUE_ H_N1_ H_N2_ H_SOURCE_IND_ H_SOURCE_MAT_
global M_TYPE_ M_ND_ M_NG_ M_NS_ M_W_ M_L_ M_VT_ M_MU_ M_COX_ M_LAMBDA_ M_CJ0_
global NMOS_ PMOS_
global METHOD_ TSTEP_ TSTOP_ ORDER_ AC_PPD_ AC_FSTART_ AC_FSTOP_
global FE_ BE_ TR_ AWE_ PRIMA_ DEC_ LIN_ DC_SWEEP
global LSTYPE_ L_CIR_ SWEEP_START SWEEP_END SWEEP_STEP SWEEP_DEV SWEEP_ELEM
global PULSE_ PLS_L_ PLS_H_ PLS_D_ PLS_R_ PLS_W_ PLS_F_ PLS_P_

global MAG_ VALUE_ PHASE_ IC_ DAC_ N1_ N2_ CN1_ CN2_ FNUM_
global MOS_MID_ MOS_W_ MOS_L_ MOS_TYPE_
global MOD_ID_ MOD_VT_ MOD_MU_ MOD_COX_ MOD_LAMBDA_ MOD_CJ0_
global M_MID_

% For the ELEM Vector
%! TYPE_  = 1; %! is being defined now in parser_init.m, MB
MAG_   = 2;
VALUE_ = 2;
PHASE_ = 3;
IC_    = 4;
DAC_   = 5;
N1_    = 6;
N2_    = 7;
CN1_   = 8;
CN2_   = 9;
FNUM_  = 10;

% for PWL independent sources
PWL_START_V_ = 8;
PWL_START_I_ = 8;

% for MOSFET
MOS_MID_    = 2;
MOS_W_      = 3;
MOS_L_      = 4;
MOS_ND_     = 5;
MOS_NG_     = 6;
MOS_NS_     = 7;
MOS_TYPE_   = 8;

M_MID_      = 8; %! from parser_init.m to avoid confusion; MB

% for Model
MOD_ID_     = 2;
MOD_VT_     = 3;
MOD_MU_     = 4;
MOD_COX_    = 5;
MOD_LAMBDA_ = 6;
MOD_CJ0_    = 7;

% element identifiers
R_ = abs('R');
M_ = abs('M');
%Y_ = abs('Y');
N_ = abs('N');
C_ = abs('C');
L_ = abs('L');
%S_ = abs('S');
%K_ = abs('K');
%W_ = abs('W');
V_ = abs('V');
I_ = abs('I');
% G_ = abs('G');
% E_ = abs('E');
% F_ = abs('F');
% H_ = abs('H');

% element type;
TYPE_       = 1; %! type ID for the element
% for resistors
R_VALUE_    = 2;
R_N1_       = 3;
R_N2_       = 4;
% for capacitors
C_VALUE_    = 2;
C_N1_       = 3;
C_N2_       = 4;
C_IC_       = 5; %! initial capacitor voltage
% for inductors
L_VALUE_    = 2;
L_N1_       = 3;
L_N2_       = 4;
L_IC_       = 5; %! initial inductor current

% independent sources
DC_         = 0;
AC_         = 1;
PWL_       = 2;
TRAN_     = 3;
DA_         = 4; %! if src not specific type, then available for both dc and ac.
SIN_        = 30;
PULSE_   = 31;
RD_     = 100;
% for voltage sources
V_VALUE_    = 2;
V_N1_       = 3;
V_N2_       = 4;
V_TYPE_     = 5; %! DC_ AC_ PWL_ SIN_
V_POINTS_   = 6;
V_PHASE_    = 6;
V_DCAMP_    = 31;
V_ACAMP_    = 32;
V_SINFREQ_  = 33;
V_SINPHASE_ = 34;
% for current sources
I_VALUE_    = 2;
I_N1_       = 3;
I_N2_       = 4;
I_TYPE_     = 5; %! DC_ AC_ PWL_
I_POINTS_   = 6;
I_PHASE_    = 6;
I_DCAMP_    = 31;
I_ACAMP_    = 32;
I_SINFREQ_  = 33;
I_SINPHASE_ = 34;

% pulse
PLS_L_  = 31; % low
PLS_H_ = 32; % high
PLS_D_  = 33; % delay
PLS_R_  = 34; % rising
PLS_W_  = 35; % width
PLS_F_  = 36; % falling
PLS_P_ = 37; % period


% % for voltage controlled voltage sources
% E_VALUE_    = 2;
% E_N1_       = 3;
% E_N2_       = 4;
% E_CN1_      = 5;
% E_CN2_      = 6;
% % for current controlled current sources
% F_VALUE_    = 2;
% F_N1_       = 3;
% F_N2_       = 4;
% F_SOURCE_IND_   = 5; %! source (controlling) element index
% F_SOURCE_MAT_   = 6; %! matrix, where to find source element data
% %! 1 = LINELEM, 2 = NLNELEM
% % for voltage controlled current sources
% G_VALUE_    = 2;
% G_N1_       = 3;
% G_N2_       = 4;
% G_CN1_      = 5;
% G_CN2_      = 6;
% % for current controlled voltage sources
% H_VALUE_    = 2;
% H_N1_       = 3;
% H_N2_       = 4;
% H_SOURCE_IND_   = 5; %! source (controlling) element index
% H_SOURCE_MAT_   = 6; %! matrix, where to find source element data
%! 1 = LINELEM, 2 = NLNELEM
% for MOSFET
M_TYPE_     = 2; %! type of MOSFET; NMOS_ or PMOS_
M_ND_       = 3; %! drain node
M_NG_       = 4; %! gate node
M_NS_       = 5; %! source node
M_W_        = 6; %! width
M_L_        = 7; %! length
% and those from the .MODEL card
M_VT_       = 8; %! threshold voltage
M_MU_       = 9; %! mobility
M_COX_      = 10; %! gate oxide capacitance per area
M_LAMBDA_   = 11; %! channel-width modulator
M_CJ0_      = 12; %! junction capacitance

NMOS_       = 1; %! type ID for NMOS (-> M_TYPE_)
PMOS_       = 0; %! type ID for PMOS (-> M_TYPE_)
% Elements of the INFO vector
METHOD_    = 1; %! integration method
TSTEP_     = 2; %! internal timestep
TSTOP_     = 3; %! end time
ORDER_     = 4; %! order of model order reduction
AC_PPD_    = 5; %! evaluation points per decade
AC_FSTART_ = 6; %! lower evaluation frequency bound
AC_FSTOP_  = 7; %! upper evaluation frequency bound
LSTYPE_	   = 8; %! inductor/susceptor circuit indicator
SIM_       = 9; %! simulation
SWEEP_START = 10;
SWEEP_END   = 11;
SWEEP_STEP  = 12;
SWEEP_DEV   = 13; %! dev idx in NAMES
SWEEP_ELEM  = 14; %! dev idx in LINELEM

% values for INFO(METHOD_) -> 0, if no transient analysis of any sort
FE_ = 1; %! Forward Euler
BE_	= 2; %! Backward Euler
TR_	= 3; %! Trapezoidal Rule
AWE_    = 4; %! Asymptotic Waveform Evaluation
PRIMA_  = 5; %! PRIMA (ask Altan what the acronym means)
DEC_ = 6;%! AC analysis freq decade sweep
LIN_ = 7;%! AC analysis freq linear sweep
DC_SWEEP = 8;
% values for INFO(LSTYPE_)
L_CIR_  = 1;
% S_CIR_  = 2;
% HB_     = 30;
% HBBASETONE_ = 31;
% HBFREQNUM_  = 32;
