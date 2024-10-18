function [bank_out, emiss_out, mf_out, bank_rf_out] = bank_emiss_mf(lt, Gg2ppt, mf_in, start_yr, end_yr, mf_obs_yr, bank_ref_yr, bank_ref, bank_rf, prod_in, prod_yrs ,inferRF, Bank2015, oa)

    Years = start_yr:end_yr;
    
    ind1  = find(Years == start_yr);
    ind2  = find(Years == bank_ref_yr);
    ind3  = find(Years == mf_obs_yr);
    ind4  = find(Years == end_yr); 
    ind2a = find(Years == 2015);


    % initiate output values
    mf_out    = nan(length(Years), 1);
    emiss_out = nan(length(Years), 1);
    bank_out  = nan(length(Years), 1);
    prod      = zeros(length(Years), 1);
    
    %initiate 'known' values
    bank_out(ind2)      = bank_ref;
    mf_out(ind1:ind3)   = mf_in(ind1:ind3);
    emiss_out(ind1:ind3-1) = (1/Gg2ppt)*(mf_in(ind1+1:ind3) - exp(-1/lt)*mf_in(ind1:ind3-1));

    ind_tmp1 = find(Years == prod_yrs(1));
    ind_tmp2 = find(Years == prod_yrs(end));

    prod(ind_tmp1:ind_tmp2) = prod_in;
    

    for ii = 1:ind2-1
        bank_out(ind2-ii) = bank_out(ind2-ii+1) - prod(ind2-ii) + emiss_out(ind2-ii);
    end
    
    for ii = ind2+1:ind3
        bank_out(ii) = bank_out(ii-1) + prod(ii-1) - emiss_out(ii-1);
    end
       
    bank_out(bank_out<0) = 0; % setting all negative bank values to zero

    if inferRF
        if oa == 1
            bank_tmp(ind3) = bank_out(ind3);
            rf_tmp = [0.001:0.001:0.1];
            for jj = 1:100
                for ii = ind3+1:ind2a
                    %bank_tmp(ii) = (1 - rf_tmp(jj))*bank_tmp(ii-1)+prod(ii-1);
                    bank_tmp(ii) = (1 - rf_tmp(jj))*(bank_tmp(ii-1)+prod(ii-1));
                end
                diff_banks(jj) = abs(Bank2015 - bank_tmp(ind2a));
            end
            
            ind_rf = find(diff_banks == min(diff_banks),1);
            bank_rf_out = rf_tmp(ind_rf)

            for ii = ind3+1:ind2a
                bank_out(ii) = bank_rf_out*(bank_out(ii-1)+prod(ii-1));
            end
            %bank_rf_out = 1 - (Bank2015/bank_out(ind3))^(1/10);
            %bank_rf_out
        elseif oa == 2
            for ii = 1:10
                rf_tmp(ii) = emiss_out(ind2-ii+1)/(bank_out(ind2-ii+1)+prod(ind2-ii+1));
            end
            
            cumsum(fliplr(rf_tmp))./[1:10];
            bank_rf_out = mean(rf_tmp)
        elseif oa == 3

            for ii = 1:7
                rf_tmp(ii) = emiss_out(ind2-2+ii)/(bank_out(ind2-2+ii)+prod(ind2-2+ii)); %2007 - 2013
            end
            cumsum(fliplr(rf_tmp))./[1:7];
            bank_rf_out = mean(rf_tmp)
        elseif oa == 4

            for ii = 1:7
                rf_tmp(ii) = emiss_out(ind2+ii+2)/(bank_out(ind2+ii+2)+prod(ind2+ii+2)); %2011 - 2017
            end
            cumsum(fliplr(rf_tmp))./[1:7];
            bank_rf_out = mean(rf_tmp)
        end

    else
        bank_rf_out = bank_rf;
    end

    if bank_out(ind3)<=0
        bank_out(ind3) = 0; 
        bank_rf_out = 1; 
    end

%     if bank_rf_out < 0
%         bank_rf_out = 1; 
%     end

    emiss_out(ind3) = (prod(ind3) + bank_out(ind3))*bank_rf_out;
    

    for ii = ind3+1:ind4
        bank_out(ii) = bank_out(ii - 1) + prod(ii-1) - emiss_out(ii-1);
        emiss_out(ii) = (prod(ii) + bank_out(ii))*bank_rf_out;
        mf_out(ii) = exp(-1/lt)*mf_out(ii-1) + Gg2ppt*emiss_out(ii-1);
    end

