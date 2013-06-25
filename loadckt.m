%function [ELEM, NAMES] = loadckt(ckt)
function [parse_ok,ELEM,INFO,NODES,NAMES,PRINTNV,PLOTNV,PLOTBI_INIT] = loadckt(ckt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%cparse_init;
parser_init;

parse_ok = 1;

ELEM=[];INFO=[];NODES=[];NAMES=[];
PRINTNV=[];PLOTNV=[];PLOTBI_INIT=[];
nodelist=[];

%global number_nunk
number_nunk = 0;
number_elem = 0;
number_name = 0;
ELEM = zeros(1,1);
INFO = zeros(1,10);
%NAMES = zeros(1,1);

%ckt = 'a.txt';
fip = fopen(ckt,'r');
if(fip == -1)
    fprintf(' Error:  Can not find the file "%s"\n', ckt);
    parse_ok = 0;
    return;
else
    fprintf(' The input ckt is "%s" \n\n', ckt);
    disp(' Load and analyze the input circuit')
end

plnv=0;
plbi=0;
prnv=0;
prbi=0;
num_row = 0;
num_option = 0;

while ~feof(fip)
    tline = deblank(fgetl(fip));
    num_row = num_row + 1;
    if ~isempty(tline)
        tletter = upper(tline(1));
        
        switch tletter
            case '*'
                
            case 'R'
                %disp('find resistor');
                number_elem = number_elem + 1;
                number_name = number_name + 1;
                rrow = regexp(tline,'\s+','split');
                
                if(max(size(rrow))~=4)
                    fprintf(' *Error at row %d: Rxx node1 node2 val\n', num_row);
                    parse_ok=0;
                    continue
                end
                
                dname = upper(char(rrow(1)));
                %upper(char(rrow(1)))
                nname1 = str2num(char(rrow(2)));
                nname2 = str2num(char(rrow(3)));
                val = str2double(char(rrow(4)));
                
                [idx1,number_nunk, nodelist] = addnode(nname1, nodelist, number_nunk);
                [idx2,number_nunk, nodelist] = addnode(nname2, nodelist, number_nunk);
                
                A=[];
                A(TYPE_) = abs(dname(1));
                A(R_VALUE_) = val;
                A(N1_) = idx1; A(N2_) = idx2;
                
                ELEM = combELEM(A, ELEM, number_elem);
                NAMES = combNAMES(dname, NAMES, number_name);
                
            case 'L'
                %disp('find inductor');
                number_elem = number_elem + 1;
                number_name = number_name + 1;
                rrow = regexp(tline,'\s+','split');
                
                if(max(size(rrow))~=4)
                    fprintf(' *Error at row %d: Lxx node1 node2 val\n', num_row);
                    parse_ok=0;
                    continue
                    %return
                end
                
                dname = upper(char(rrow(1)));
                nname1 = str2num(char(rrow(2)));
                nname2 = str2num(char(rrow(3)));
                val = str2double(char(rrow(4)));
                
                [idx1,number_nunk, nodelist] = addnode(nname1, nodelist, number_nunk);
                [idx2,number_nunk, nodelist] = addnode(nname2, nodelist, number_nunk);
                
                A=[];
                A(TYPE_) = abs(dname(1));
                A(L_VALUE_) = val;
                A(N1_) = idx1; A(N2_) = idx2;
                
                ELEM = combELEM(A, ELEM, number_elem);
                NAMES = combNAMES(dname, NAMES, number_name);
                
            case 'C'
                %disp('find capacitor');
                number_elem = number_elem + 1;
                number_name = number_name + 1;
                rrow = regexp(tline,'\s+','split');
                
                if(max(size(rrow))~=4)
                    fprintf(' *Error at row %d: Cxx node1 node2 val\n', num_row);
                    parse_ok = 0;
                    continue
                    %return
                end
                
                dname = upper(char(rrow(1)));
                nname1 = str2num(char(rrow(2)));
                nname2 = str2num(char(rrow(3)));
                val = str2double(char(rrow(4)));
                
                [idx1,number_nunk, nodelist] = addnode(nname1, nodelist, number_nunk);
                [idx2,number_nunk, nodelist] = addnode(nname2, nodelist, number_nunk);
                
                A=[];
                A(TYPE_) = abs(dname(1));
                A(C_VALUE_) = val;
                A(IC_) = 0/0;
                A(N1_) = idx1; A(N2_) = idx2;
                
                ELEM = combELEM(A, ELEM, number_elem);
                NAMES = combNAMES(dname, NAMES, number_name);
                
            case 'V'
                number_elem = number_elem + 1;
                number_name = number_name + 1;
                rrow = regexp(tline,'\s+','split');
                
                if(max(size(rrow))<4)
                    fprintf(' *Error at row %d: Vxx node1 node2 type val ...\n', num_row);
                    parse_ok=0;
                    continue
                    %return
                    %elseif(max(size(rrow))==4)
                end
                
                dname = upper(char(rrow(1)));
                nname1 = str2num(char(rrow(2)));
                nname2 = str2num(char(rrow(3)));
                
                vtype=[];vvtype=[];
                if(max(size(rrow))>4)
                    vtype = char(rrow(4));
                    vvtype = upper(vtype);
                end
                %val = str2double(char(rrow(5)));
                
                [idx1,number_nunk, nodelist] = addnode(nname1, nodelist, number_nunk);
                [idx2,number_nunk, nodelist] = addnode(nname2, nodelist, number_nunk);
                
                A=[];
                A(TYPE_) = abs(dname(1));
                A(N1_) = idx1; A(N2_) = idx2;
                
                if(strcmp(vvtype,'DC') && max(size(rrow))>4)
                    if(max(size(rrow))~=5)
                        fprintf(' *Error at row %d: incorrect DC syntax\n', num_row);
                        parse_ok = 0;
                        continue
                    end
                    A(V_TYPE_) = DC_;
                    A(V_VALUE_) = str2double(char(rrow(5)));
                    
                elseif(strcmp(vvtype,'AC') && max(size(rrow))>4)
                    if(max(size(rrow))~=5)
                        fprintf(' *Error at row %d: incorrect AC syntax\n', num_row);
                        parse_ok = 0;
                        continue
                    end
                    A(V_TYPE_) = AC_;
                    A(V_VALUE_) = str2double(char(rrow(5)));
                    
                elseif(strcmp(vvtype,'PWL') && max(size(rrow))>4)
                    A(V_TYPE_) = PWL_;
                    A(V_VALUE_) = str2double(char(rrow(6)));
                    
                    pwlpair = size(rrow,2) - 4;
                    
                    if(mod(pwlpair,2) ~= 0)
                        fprintf(' *Error at row %d: incorrect PWL input\n', num_row);
                        parse_ok=0;
                        continue
                    else
                        pwlcnt = pwlpair/2;
                        A(PWL_START_V_) = pwlcnt;
                    end
                    
                    cnt = PWL_START_V_;
                    for i=5:size(rrow,2)
                        cnt = cnt + 1;
                        A(cnt) = str2double(char(rrow(i)));
					end
					for i=PWL_START_V_+1:2:cnt-3
						if(A(i)>=A(i+2))
                        fprintf(' *Error at row %d: incorrect time points squence\n', num_row);
							parse_ok=0;
							break
						end
					end
                    
                elseif(strcmp(vvtype,'SIN') && max(size(rrow))>4)
                    if(max(size(rrow))~=8)
                        fprintf(' *Error at row %d: incorrect SIN syntax\n', num_row);
                        parse_ok = 0;
                        continue
                    end
                    
                    A(V_TYPE_) = SIN_;
                    A(V_VALUE_) = str2double(char(rrow(5)));
                    A(V_DCAMP_) = str2double(char(rrow(5)));
                    A(V_ACAMP_) = str2double(char(rrow(6)));
                    A(V_SINFREQ_) = str2double(char(rrow(7)));
                    A(V_SINPHASE_) = str2double(char(rrow(8)));
                    
                elseif(strcmp(vvtype,'PULSE') && max(size(rrow))>4)
                    if(max(size(rrow))~=11)
                        fprintf(' *Error at row %d: incorrect PULSE syntax\n', num_row);
                        parse_ok = 0;
                        continue
                    end
                    
                    A(V_TYPE_) = PULSE_;
                    A(V_VALUE_) = str2double(char(rrow(5)));
                    A(PLS_L_) = str2double(char(rrow(5)));
                    A(PLS_H_) = str2double(char(rrow(6)));
                    A(PLS_D_) = str2double(char(rrow(7)));
                    A(PLS_R_) = str2double(char(rrow(8)));
                    A(PLS_W_) = str2double(char(rrow(9)));
                    A(PLS_F_) = str2double(char(rrow(10)));
                    A(PLS_P_) = str2double(char(rrow(11)));
                    
                elseif(max(size(rrow))==4)
                    A(V_TYPE_) = DA_;
                    A(V_VALUE_) = str2double(char(rrow(4)));
                    
                else
                    fprintf(' *Error at row %d: unknown voltage type\n', num_row);
                    parse_ok = 0;
                    continue
                end
                
                ELEM = combELEM(A, ELEM, number_elem);
                NAMES = combNAMES(dname, NAMES, number_name);
                
            case 'M'
                number_elem = number_elem + 1;
                number_name = number_name + 1;
                rrow = regexp(tline,'\s+','split');
                
                if(max(size(rrow))~=8)
                    fprintf(' *Error at row %d: Mxx n_d n_g n_s type ch_w ch_l mod\n', num_row);
                    parse_ok=0;
                    continue
                    %return
                end
                
                dname = upper(char(rrow(1)));
                nname1 = str2num(char(rrow(2)));
                nname2 = str2num(char(rrow(3)));
                nname3 = str2num(char(rrow(4)));
                mtype = upper(char(rrow(5)));
                ch_w = str2double(char(rrow(6)));
                ch_l = str2double(char(rrow(7)));
                midx = str2double(char(rrow(8)));
                
                [idx1,number_nunk, nodelist] = addnode(nname1, nodelist, number_nunk);
                [idx2,number_nunk, nodelist] = addnode(nname2, nodelist, number_nunk);
                [idx3,number_nunk, nodelist] = addnode(nname3, nodelist, number_nunk);
                
                A=[];
                A(TYPE_) = abs(dname(1));
                A(MOS_MID_) = midx;
                A(MOS_ND_) = idx1; A(MOS_NG_) = idx2; A(MOS_NS_) = idx3;
                A(MOS_W_) = ch_w; A(MOS_L_) = ch_l;
                if(strcmp(mtype,'N'))
                    A(MOS_TYPE_) = 1;
                elseif(strcmp(mtype,'P'))
                    A(MOS_TYPE_) = 0;
				else
					fprintf(' *Error: unknown MOS type at row %d\n', num_row);
					parse_ok = 0;
					continue
                end
                
                ELEM = combELEM(A, ELEM, number_elem);
                NAMES = combNAMES(dname, NAMES, number_name);
            
            case 'I'
                number_elem = number_elem + 1;
                number_name = number_name + 1;
                rrow = regexp(tline,'\s+','split');
                
                if(max(size(rrow))<4)
                    fprintf(' *Error at row %d: Ixx node1 node2 type val ...\n', num_row);
                    parse_ok=0;
                    continue
                    %return
                    %elseif(max(size(rrow))==4)
                end
                
                dname = upper(char(rrow(1)));
                nname1 = str2num(char(rrow(2)));
                nname2 = str2num(char(rrow(3)));
                
                itype=[];iitype=[];
                if(max(size(rrow))>4)
                    itype = char(rrow(4));
                    iitype = upper(itype);
                end
                %val = str2double(char(rrow(5)));
                
                [idx1,number_nunk, nodelist] = addnode(nname1, nodelist, number_nunk);
                [idx2,number_nunk, nodelist] = addnode(nname2, nodelist, number_nunk);
                
                A=[];
                A(TYPE_) = abs(dname(1));
                A(N1_) = idx1; A(N2_) = idx2;
                
                if(strcmp(iitype,'DC') && max(size(rrow))>4)
                    if(max(size(rrow))~=5)
                        fprintf(' *Error at row %d: incorrect DC syntax\n', num_row);
                        parse_ok = 0;
                        continue
                    end
                    A(I_TYPE_) = DC_;
                    A(I_VALUE_) = str2double(char(rrow(5)));
                    
                elseif(strcmp(iitype,'AC') && max(size(rrow))>4)
                    if(max(size(rrow))~=5)
                        fprintf(' *Error at row %d: incorrect AC syntax\n', num_row);
                        parse_ok = 0;
                        continue
                    end
                    A(I_TYPE_) = AC_;
                    A(I_VALUE_) = str2double(char(rrow(5)));
                    
                elseif(strcmp(iitype,'PWL') && max(size(rrow))>4)
                    A(I_TYPE_) = PWL_;
                    A(I_VALUE_) = str2double(char(rrow(6)));
                    
                    pwlpair = size(rrow,2) - 4;
                    
                    if(mod(pwlpair,2) ~= 0)
                        fprintf(' *Error at row %d: incorrect PWL input\n', num_row);
                        parse_ok=0;
                        continue
                    else
                        pwlcnt = pwlpair/2;
                        A(PWL_START_I_) = pwlcnt;
                    end
                    
                    cnt = PWL_START_I_;
                    for i=5:size(rrow,2)
                        cnt = cnt + 1;
                        A(cnt) = str2double(char(rrow(i)));
					end
					for i=PWL_START_I_+1:2:cnt-3
						if(A(i)>=A(i+2))
                        fprintf(' *Error at row %d: incorrect time points squence\n', num_row);
							parse_ok=0;
							break
						end
					end

                    
                elseif(strcmp(iitype,'SIN') && max(size(rrow))>4)
                    if(max(size(rrow))~=8)
                        fprintf(' *Error at row %d: incorrect SIN syntax\n', num_row);
                        parse_ok = 0;
                        continue
                    end
                    A(I_TYPE_) = SIN_;
                    A(I_VALUE_) = str2double(char(rrow(5)));
                    A(I_DCAMP_) = str2double(char(rrow(5)));
                    A(I_ACAMP_) = str2double(char(rrow(6)));
                    A(I_SINFREQ_) = str2double(char(rrow(7)));
                    A(I_SINPHASE_) = str2double(char(rrow(8)));
                    
                elseif(strcmp(iitype,'PULSE') && max(size(rrow))>4)
                    if(max(size(rrow))~=11)
                        fprintf(' *Error at row %d: incorrect PULSE syntax\n', num_row);
                        parse_ok = 0;
                        continue
                    end
                    
                    A(I_TYPE_) = PULSE_;
                    A(I_VALUE_) = str2double(char(rrow(5)));
                    A(PLS_L_) = str2double(char(rrow(5)));
                    A(PLS_H_) = str2double(char(rrow(6)));
                    A(PLS_D_) = str2double(char(rrow(7)));
                    A(PLS_R_) = str2double(char(rrow(8)));
                    A(PLS_W_) = str2double(char(rrow(9)));
                    A(PLS_F_) = str2double(char(rrow(10)));
                    A(PLS_P_) = str2double(char(rrow(11)));
                    
                elseif(max(size(rrow))==4)
                    A(I_TYPE_) = DA_;
                    A(I_VALUE_) = str2double(char(rrow(4)));
                    
                else
                    fprintf(' *Error at row %d: unknown current type\n', num_row);
                    parse_ok = 0;
                    continue
                end
                
                ELEM = combELEM(A, ELEM, number_elem);
                NAMES = combNAMES(dname, NAMES, number_name);
                
            case '.'
                tletter2 = upper(tline(2));
                switch tletter2
                    case 'M'
                        number_elem = number_elem + 1;
                        number_name = number_name + 1;
                        rrow = regexp(tline,'\s+','split');
                        
                        if(max(size(rrow))~=12)
                            fprintf(' *Error at row %d: .MODEL x VT x MU x COX x LAMBDA x CJ0 x\n', num_row);
                            parse_ok=0;
                            continue
                            %return
                        end
                        
                        dname = upper(char(rrow(1)));
                        mmtype = str2num(char(rrow(2)));
                        vt = str2double(char(rrow(4)));
                        mu = str2double(char(rrow(6)));
                        cox = str2double(char(rrow(8)));
                        lambda = str2double(char(rrow(10)));
                        cj0 = str2double(char(rrow(12)));
                        
                        A=[];
                        A(TYPE_) = abs(dname(1));
                        A(MOD_ID_) = mmtype;
                        A(MOD_VT_) = vt;
                        A(MOD_MU_) = mu;
                        A(MOD_COX_) = cox;
                        A(MOD_LAMBDA_) = lambda;
                        A(MOD_CJ0_) = cj0;
                        
                        ELEM = combELEM(A, ELEM, number_elem);
                        NAMES = combNAMES(dname, NAMES, number_name);
                        
                    case 'P'
                        rrow = regexp(tline,'\s+','split');
                        %print or plot and voltage or current
                        %pnname = char(pline(2));
                        ptype = upper(char(rrow(1)));
                        
                        if(strcmp(ptype(3),'R')) %print
                            if(strcmp(ptype(7:8),'NV'))
                                prnv = prnv+1;
                                nname = str2num(char(rrow(2)));
                                [parse_ok, node_idx, number_nunk, nodelist] = findnode(nname,nodelist,number_nunk,1,parse_ok);
                                if(parse_ok==0)
                                    continue
                                end
                                PRINTNV(prnv) = node_idx;
                            elseif(strcmp(ptype(7:8),'BI'))
                                %disp('printbi')
                                prbi = prbi + 1;
                                dname = char(rrow(2));
                                [parse_ok,didx] = findelem(dname, NAMES, number_name,parse_ok);
                                PRINTBI_INIT(prbi) = didx;
                            else
                                fprintf(' *Skip at row %d: unknown option `%s`\n', num_row, ptype);
                            end
                        elseif(strcmp(ptype(3),'L')) %plot
                            if(strcmp(ptype(6:7),'NV'))
                                %disp('plotnv')
                                plnv = plnv+1;
                                nname = str2num(char(rrow(2)));
                                [parse_ok, node_idx, number_nunk, nodelist] = findnode(nname,nodelist,number_nunk,1,parse_ok);
                                if(parse_ok==0)
                                    continue
                                end
                                PLOTNV(plnv) = node_idx;
                            elseif(strcmp(ptype(6:7),'BI'))
                                %disp('plotbi')
                                plbi = plbi + 1;
                                dname = char(rrow(2));
                                [parse_ok,didx] = findelem(dname, NAMES, number_name,parse_ok);
                                PLOTBI_INIT(plbi) = didx;
                            else
                                fprintf(' *Skip at row %d: unknown option `%s`\n', num_row, ptype);
                            end
                        end
                        
                    case 'T'
                        num_option = num_option + 1;
                        rrow = regexp(tline,'\s+','split');
                        if(max(size(rrow)) == 4)
                            dmod = upper(char(rrow(2)));
                            a = str2double(char(rrow(3)));
                            b = str2double(char(rrow(4)));
                            t_h = min(a,b);
                            tend = max(a,b);
                        elseif(max(size(rrow))== 3)
                            a = str2double(char(rrow(2)));
                            b = str2double(char(rrow(3)));
                            t_h = min(a,b);
                            tend = max(a,b);
                        else
                            fprintf(' *Error at row %d: .TRAN type step stop\n', num_row);
                            %fclose(fip);
                            parse_ok=0;
                            continue
                        end
                        
                        if(max(size(rrow))== 4)
                            if(strcmp(dmod, 'TR'))
                                INFO(METHOD_) = TR_;
                            elseif(strcmp(dmod, 'BE'))
                                INFO(METHOD_) = BE_;
                            else
                                fprintf(' *Error at row %d: .TRAN type(BE,TR) step stop\n', num_row);
                                %fclose(fip);
                                parse_ok=0;
                                continue
                            end
                        elseif(max(size(rrow))== 3)
                            INFO(METHOD_) = BE_;
                        end
                        
                        INFO(METHOD_) = BE_;
                        
                        INFO(TSTEP_) = t_h;
                        INFO(TSTOP_) = tend;
                        INFO(SIM_) = TRAN_;
                        
                    case 'D'
                        num_option = num_option + 1;
                        rrow = regexp(tline,'\s+','split');
                        
                        INFO(SIM_) = DC_;
                        
                        if(max(size(rrow)) == 6 && strcmp(char(rrow(2)),'SWEEP'))
                            var = upper(char(rrow(3)));
                            sweep_start = str2double(char(rrow(4)));
                            sweep_end = str2double(char(rrow(5)));
                            sweep__h = str2double(char(rrow(6)));
                            
                            INFO(METHOD_) = DC_SWEEP;
                            INFO(SWEEP_START) = sweep_start;
                            INFO(SWEEP_END) = sweep_end;
                            INFO(SWEEP_STEP) = sweep__h;
                            [parse_ok,didx] = findelem(var,NAMES,number_name,parse_ok);
                            INFO(SWEEP_DEV) = didx;
                        elseif(max(size(rrow)) > 6 && strcmp(char(rrow(2)),'SWEEP'))
                            fprintf(' *Warning: only one DC SWEEP variables allowed\n');
                            parse_ok = 0;
                            %fprintf('          reset to DC only.\n');
                        elseif(max(size(rrow))~=1)
                            fprintf(' *Error at row %d: .DC or .DC SWEEP val start stop step\n', num_row);
                            parse_ok = 0;
                        end
                        
                        
                    case 'A'
                        num_option = num_option + 1;
                        rrow = regexp(tline,'\s+','split');
                        if(max(size(rrow)) == 5)
                            dmod = upper(char(rrow(2)));
                            f_h = str2double(char(rrow(3)));
                            fstart = str2double(char(rrow(4)));
                            fend = str2double(char(rrow(5)));
                        elseif(max(size(rrow)) == 4)
                            f_h = str2double(char(rrow(2)));
                            fstart = str2double(char(rrow(3)));
                            fend = str2double(char(rrow(4)));
                        else
                            fprintf(' *Error at row %d: .AC type(LIN,DEC) points start stop\n', num_row);
                            %fclose(fip);
                            parse_ok=0;
                            continue
                        end
                        
                        if(max(size(rrow)) == 5)
                            if(strcmp(dmod, 'DEC'))
                                INFO(METHOD_) = DEC_;
                            elseif(strcmp(dmod, 'LIN'))
                                INFO(METHOD_) = LIN_;
                            else
                                fprintf(' *Error at row %d: unknown AC freq sweep option (LIN or DEC)\n', num_row);
                                %fclose(fip);
                                parse_ok=0;
                                continue
                            end
                        elseif(max(size(rrow)) == 4)
                            INFO(METHOD_) = LIN_;
                        end
                        
                        INFO(SIM_) = AC_;
                        INFO(AC_PPD_) = f_h;
                        INFO(AC_FSTART_) = fstart;
                        INFO(AC_FSTOP_) = fend;
                       
					case 'R'
						
                        num_option = num_option + 1;

                        rrow = regexp(tline,'\s+','split');
                        ptype = upper(char(rrow(2)));
                        
                        if(strcmp(ptype,'PRIMA'))
 	                       	INFO(METHOD_) = PRIMA_;
    	                    
                	        INFO(SIM_) = RED_;
						else
 	                       	fprintf(' *Error at row %d: unknown MOR option\n', num_row);
        	                parse_ok=0;
            	            continue
						end

                    otherwise
                        fprintf(' *Error at row %d: unknown option\n', num_row);
                        parse_ok=0;
                        continue
                end
            otherwise
                fprintf(' *Error at row %d: unknown device\n', num_row)
                %fclose(fip);
                parse_ok=0;
                continue
        end
    end
    
end

fclose(fip);

if(num_option == 0)
    fprintf(' *Error: simulation type is required (.DC, .AC, or .TRAN)!!\n');
    fprintf('         setup info can be found in cktsim(''help'')\n');
    parse_ok=0;
    return
elseif(num_option > 1)
    fprintf(' *Error: multiple simulation methods are set in the file!!\n');
    parse_ok=0;
    return
end

if(parse_ok == 1)
    kk = zeros(size(ELEM,1),3);
    ELEM = [ELEM,kk];
    NODES = nodelist;
    
    % flip some vectors
    INFO = INFO';
    NODES = NODES';
    PRINTNV = PRINTNV';
    PLOTNV = PLOTNV';
    PLOTBI_INIT = PLOTBI_INIT';
    
    fprintf(' >  ckt includes %d devices\n', number_elem);
    fprintf(' >  system has %d nodes\n', number_nunk);
else
    fprintf('\n');
    fprintf(' Please check the syntax of each device\n');
    fprintf(' Please check the sequence of the ckt list\n');
    fprintf('   - Place all the options at the end of the list\n');
end

end


function [nodeidx, number_nunk, nodelist] = addnode(nname, nodelist, number_nunk)

nflag=0;

if(number_nunk == 0 && nname ~= 0)
    number_nunk = number_nunk + 1;
    nodelist(number_nunk) = nname;
    nodeidx = number_nunk;
elseif(number_nunk ~= 0 && nname ~= 0)
    for i=1:number_nunk
        if(nname == nodelist(i))
            nflag = nflag+1;
            nodeidx = i;
        end
    end
    
    if(nflag > 1)
        fprintf(' Duplicated nodes already stored in nodelist\n');
    elseif(nflag == 0)
        if(nname ~= 0)
            number_nunk = number_nunk + 1;
            nodelist(number_nunk) = nname;
            nodeidx = number_nunk;
        else
            nodeidx = -1;
        end
    end
else
    nodeidx = -1;
end


end

function [parse_ok,nodeidx, number_nunk, nodelist] = findnode(nname, nodelist, number_nunk, option, parse_ok)
nflag = 0;
nodeidx = 1;

if(parse_ok == 0)
    return
end

for i=1:number_nunk
    if(nname == nodelist(i))
        nodeidx = i;
        nflag = 1;
        break;
    end
end

if(nflag == 0 && option==0)
    number_nunk = number_nunk + 1;
    nodelist(number_nunk) = nname;
    nodeidx = number_nunk;
elseif(nflag == 0 && option==1)
    fprintf(' Error: can not find node (%d) in option\n', nname);
    parse_ok = 0;
end

end


function [NAMES] = combNAMES(dname, NAMES, number_elem)
s_d = size(dname,2);
s_n = size(NAMES,2);
kk=[];
if(s_d > s_n)
    for i=1:(s_d-s_n)
        kk = [kk, ' '];
    end
    
    for i=1:size(NAMES,1)-1
        kk = [kk; kk(1,:)];
    end
    %kk = zeros(size(NAMES,1),s_d-s_n);
    NAMES = [NAMES,kk];
elseif(s_d < s_n)
    for i=1:(s_n-s_d)
        kk = [kk ' '];
    end
    %kk = zeros(1,s_n-s_d);
    dname = [dname,kk];
end
NAMES(number_elem, :) = dname;
end

function [ELEM] = combELEM(A, ELEM, number_elem)
s_a = size(A,2);
s_e = size(ELEM,2);
if(s_a > s_e)
    kk = zeros(size(ELEM,1),s_a-s_e);
    ELEM = [ELEM,kk];
elseif(s_a < s_e)
    kk = zeros(1,s_e-s_a);
    A = [A,kk];
end
ELEM(number_elem,:) = A;
end

function [parse_ok,didx] = findelem(dname, NAMES, number_name, parse_ok)
nflag = 0;
didx = 1;

if(parse_ok == 0)
    return
end

for i=1:number_name
    mname = deblank(NAMES(i,:));
    if(strcmp(dname, mname))
        didx = i;
        nflag = 1;
        break;
    end
end

if(nflag == 0)
    fprintf(' Error: can not find device (%s) in option\n', dname);
    parse_ok=0;
    return;
end
end
