function [ Res_nv, Res_bi, freq ] = ac_sim(dc_point,ac_type,f_h,fstart,fend)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% ac_sim: ac simulation
%%
%% - dc_point: dc solution
%% - ac_type : ac type (dec or lin)
%% - f_h     : frequence step size
%% - fstart  : start frequence
%% - fend    : stop frequence
%%
%% by xueqian 06/24/2012
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DEC_ LIN_ ac_ok numNodes freq_t
global plotbi plotnv G F X vin_amp

fprintf('**************************************************\n');
fprintf('   AC simulation starts ...\n   ');

ac_ok=0;
setup=0;

freq=[];
sweep =0;
vin_amp=0;

% set the dc operating point
X = dc_point;

if(ac_type == DEC_)
    dec=log10(fend/fstart);
    if(ceil(dec)==dec)
        sweep = dec * f_h + 1;
        freq = zeros(sweep,1);
        bot = fstart;
        for i=1:dec
            top = ceil(10^(log10(fstart)+i));
            for j=1:f_h
                freq(f_h * (i-1) + j) = bot + (top-bot)/f_h*(j-1);
            end
            bot = top;
        end
        freq(sweep) = fend;
    else
        sweep = ceil(dec) * f_h + 1;
        freq = zeros(sweep,1);
        bot = fstart;
        for i=1:floor(dec)
            top = ceil(10^(log10(fstart)+i));
            for j=1:f_h
                freq(f_h * (i-1) + j) = bot + (top-bot)/f_h*(j-1);
            end
            bot = top;
        end
        top = fend;
        for j=1:f_h
            freq(f_h * (ceil(dec)-1) + j) = bot + (top-bot)/f_h*(j-1);
        end
        freq(sweep)=fend;
    end
elseif(ac_type == LIN_)
    sweep = f_h+1;
    freq = zeros(sweep,1);
    for i=1:sweep-1
        freq(i) = fstart + (fend-fstart)/f_h*(i-1);
    end
    freq(sweep) = fend;
end

Res_bi = zeros(sweep,size(plotbi,1));%! save complex number
Res_nv = zeros(sweep,size(plotnv,1));%

t_ac = cputime;

% AC kernel
for k=1:sweep
    freq_t = freq(k);
    evaluate(numNodes,setup);
    
    D = G\F;
    
    for n=1:size(plotnv,1)
        if(vin_amp~=0)
            Res_nv(k,n) = D(plotnv(n))/vin_amp;
        else
            Res_nv(k,n) = D(plotnv(n));
        end
    end
    
    if(mod(k,4)==0)
        fprintf('.');
    end
    if(mod(k,100)==0)
        fprintf('\n   ');
    end
end

t_ac = cputime - t_ac;
ac_ok = 1;

fprintf('\n     finished!\n');
fprintf('   (%d) sweeps for AC analysis \n',sweep);
fprintf('   CPU time for AC analysis is %.4f(s) \n',t_ac);

fprintf('**************************************************\n');

end

