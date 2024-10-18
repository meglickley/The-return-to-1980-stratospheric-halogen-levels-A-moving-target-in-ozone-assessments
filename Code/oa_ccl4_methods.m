function [mf_out] = oa_ccl4_methods(lt, ccl4_mf, mf_yrs, Yend, oa)

ppt_per_kg = 10^(-5)*(1/25313.18968);
Gg2ppt = ppt_per_kg*(lt*(1-exp(-1/lt)));
if oa == 1
  
    %Top-down derived emissions, linear decline from 2005 to present. 
    ind = find(mf_yrs == Yend);

    if Yend == 2005 % If years are before 2015, linear interpolation to zero
        ind
        lt
        Emiss_2005 = (1/Gg2ppt)*(ccl4_mf(ind+1) - exp(-1/lt)*ccl4_mf(ind));
        Emiss_2015 = 0;
        x = interp1([2005, 2015],[Emiss_2005, Emiss_2015],[2005:2014]);
        length(x)
        ind2 = find(mf_yrs <= 2005);
        mf_out = nan(150,1);
        mf_out(ind2) = ccl4_mf(ind2);

        for ii = 1:length(x)
            mf_out(ind+ii) = exp(-1/lt)*mf_out(ind + ii - 1) + Gg2ppt*x(ii);
        end
    
        ind3 = ind + length(x);
        for ii = ind3 + 1:length(mf_out)
            mf_out(ii) = exp(-1/lt)*mf_out(ii-1);
        end
    elseif Yend == 2022

        ind2 = find(mf_yrs <= Yend);
        mf_out = nan(150,1);
        mf_out(ind2) = ccl4_mf(ind2);

        for ii = ind+1:length(mf_out)
            mf_out(ii) = exp(-1/lt)*mf_out(ii-1);
        end
    end
elseif oa == 2

    ind = find(mf_yrs == Yend);
    ind2 = find(mf_yrs == 2050);
    ind3 = find(mf_yrs == 2100);
    
    ind1_emiss_inf = find(mf_yrs == 2003);
    % Infer emissions from 2003 to end of available mf observations
    
    Emiss_5yrs = (1/Gg2ppt)*(ccl4_mf(ind1_emiss_inf+1:ind) - exp(-1/lt)*ccl4_mf(ind1_emiss_inf:ind-1));
    Emiss_5yrs
    bvals = Emiss_5yrs(2:6)./Emiss_5yrs(1:5)-1 
    % infer the b value just from 2003 to 2009 emissions values.
    b = mean(bvals)
    
    b = round(b,2)
    ind_obs = find(mf_yrs <= Yend);
    mf_out = nan(150,1);
    mf_out(ind_obs) = ccl4_mf(ind_obs);
    
    % Infer starting emissions for forward simulation from end of observational record. 
    Emiss_0 = Emiss_5yrs(end)*(1+b);
    Emiss_5yrs(end)
    Emiss_0
    for ii = ind:ind2
        mf_out(ii) = exp(-1/lt)*mf_out(ii-1) + Gg2ppt*Emiss_0;
        Emiss_0 = Emiss_0*(1+b);
    end
     Emiss_0   
    for ii = ind2+1:ind3
        mf_out(ii) = exp(-1/lt)*mf_out(ii-1);
    end
elseif oa == 3
    ind = find(mf_yrs == Yend);
    ind2 = find(mf_yrs == 2100);
    
    ind1_emiss_inf = find(mf_yrs == 2003);
    % Infer emissions from 2003 to end of available mf observations
    
    Emiss_5yrs = (1/Gg2ppt)*(ccl4_mf(ind1_emiss_inf+1:ind) - exp(-1/lt)*ccl4_mf(ind1_emiss_inf:ind-1));
    Emiss_5yrs
    if lt == 26
        b = -0.064; 
    else
        b = -0.07;
    end
    %bvals = Emiss_5yrs(2:6)./Emiss_5yrs(1:5)-1 
    % infer the b value just from 2003 to 2009 emissions values.
    %b = mean(bvals)
    
    %b = round(b,2)
    ind_obs = find(mf_yrs <= Yend);
    mf_out = nan(150,1);
    mf_out(ind_obs) = ccl4_mf(ind_obs);
    
    % Infer starting emissions for forward simulation from end of observational record. 
    Emiss_0 = Emiss_5yrs(end)*(1+b);
    Emiss_5yrs(end)
    Emiss_0
    
    Emiss_0 = Emiss_5yrs(end);
    mf_out(ind-1:ind)
    for ii = ind:ind2-1
        mf_out(ii+1) = exp(-1/lt)*mf_out(ii) + Gg2ppt*Emiss_0;
        Emiss_0 = Emiss_0*(1+b);
    end
elseif oa == 4
    ind = find(mf_yrs == Yend);
    ind2 = find(mf_yrs == 2100);
    Emiss_5yrs = (1/Gg2ppt)*(ccl4_mf(ind+1-21:ind) - exp(-1/lt)*ccl4_mf(ind-21:ind-1));
    size(Emiss_5yrs(2:end)./Emiss_5yrs(1:end-1)-1)
    b = -0.025; %mean(Emiss_5yrs(3:end)./Emiss_5yrs(2:end-1))-1

    ind_obs = find(mf_yrs <= Yend);
    mf_out = nan(150,1);
    mf_out(ind_obs) = ccl4_mf(ind_obs);
    
    Emiss_0 = Emiss_5yrs(end);
    for ii = ind:ind2-1
        mf_out(ii+1) = exp(-1/lt)*mf_out(ii) + Gg2ppt*Emiss_0;
        Emiss_0 = Emiss_0*(1+b);
    end
end



