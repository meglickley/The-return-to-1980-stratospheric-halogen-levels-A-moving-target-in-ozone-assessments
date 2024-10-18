% Kicking the Can in the Ozone Assessment
% cfc-11 analysis
clear all
close all

run_mf_update = 1; 
run_make_final_figures = 1; 
% Colors for lines for each scenario and update
LineCols = [0,      0.4470, 0.7410;
            0.8500, 0.3250, 0.0980;
            0.9290, 0.6940, 0.1250;
            0.4940, 0.1840, 0.5560;
            0.4660, 0.6740, 0.1880];

if run_mf_update
for gas_ii = 1:3
    if gas_ii == 1
        cfc11 = 1;
        cfc12 = 0;
        halon1301 = 0;
    elseif gas_ii == 2
        cfc11 = 0;
        cfc12 = 1;
        halon1301 = 0;
    elseif gas_ii == 3
         cfc11 = 0;
        cfc12 = 0;
        halon1301 = 1;
    end
    if cfc11
        load('cfc11_OA_data.mat')
        mole_str = 'cfc11';
        lt = [45, 45, 52, 52, 52];
        lt_2022 = 52; %54.3846;
        ppt_per_kg = 10^(-5)*[4.4243, 4.4237, 4.4237, 4.4237, 4.4237];
        Yend = [2006, 2010, 2014, 2018];
        
        OAproduction = cfc11_production; 
        OAmolefractions = cfc11_mf;
    
        Bank_RF = [0.033, 0.05, 0.05, 0.05];
        Bank_refyear = [2002, 2008, 2008, 2008];
        Bank_ref_val = [1653596, 1420017, 1420017, 1420017];
        Bank2015 = 1101137; % used only for the 2006 assessment
    elseif cfc12
        load('cfc12_OA_data.mat')
        mole_str = 'cfc12';
        lt_2022 = 102; %94.7364;
        % FIX THE LIFETIMES!
        lt = [100, 100, 102, 102, 102];
        ppt_per_kg = 10^(-5)*[5.0263, 5.0257, 5.0257, 5.0257, 5.0257];
        Yend = [2006, 2010, 2014, 2018];
        
        OAproduction = cfc_12_production; 
        OAmolefractions = cfc12_mf;
    
        Bank_RF = [0.0385, 0.15, 0.15, 0.15];
        Bank_refyear = [2002, 2008, 2008, 2008];
        Bank_ref_val = [710806, 393999, 393999, 393999];
        Bank2015 = 278471; % used only for the 2006 assessment
    elseif halon1301
        load('halon1301_OA_data.mat')
        mole_str = 'halon1301';
        lt_2022 = 72; 
        % FIX THE LIFETIMES!
        lt = [65, 65, 72, 72, 72];
        ppt_per_kg = 10^(-5)*[4.08112, 4.0808, 4.0808, 4.0808, 4.0808];
        Yend = [2006, 2010, 2014, 2018];
        
        OAproduction = halon1301_prod; 
        OAmolefractions = halon1301_mf;
    
        Bank_RF = [0.05, 0.04, 0.04, 0.04];
        Bank_refyear = [2002, 2008, 2008, 2008];
        Bank_ref_val = [41977, 46700, 46700, 46700];
        Bank2015 = 23538; % used only for the 2006 assessment
    end
    
    Bank_byscen = nan(151,5,4); 
    Emiss_byscen = nan(151,5,4); 
    MF_byscen = nan(151,5,4); 
    
    
    Prod(31:151,1:3) = OAproduction;
    Prod(:,4:5) = [Prod(:,3),Prod(:,3)];
    y1 = 1950;
    fighandle = figure;
    % Scenario 1: ozone assessment mf
    for scen_ii = 1:4
        if scen_ii == 1
            title_str= 'Scenario 1: Initial Assessment';
            MF_scen = OAmolefractions;  % sanity check against MF_byscen
            for oa = 1:4
    
                mf_obs_yr = Yend(oa);  % Edit this for each scenario
                LT = lt(oa);          % Edit this for each scenario
    
                Gg2ppt = ppt_per_kg(oa)*(LT*(1-exp(-1/LT)));
                mf_in = OAmolefractions(:, oa+1);
                
                start_yr =1950; 
                end_yr = 2100;
                
                bank_ref_yr = Bank_refyear(oa);
                bank_ref = Bank_ref_val(oa);
                bank_rf = Bank_RF(oa);
                
                prod_in = Prod(:,oa+1);
                prod_yrs = 1950:2100;
                inferRF =  1;

%                 if scen_ii == 1 && gas_ii == 3
%                     inferRF = 0; 
%                 else
%                     inferRF = 1; 
%                 end

                [bank_out, emiss_out, mf_out, bank_rf_out] = bank_emiss_mf(LT, Gg2ppt, mf_in, start_yr, end_yr, mf_obs_yr, bank_ref_yr, bank_ref, bank_rf, prod_in, prod_yrs, inferRF, Bank2015, oa);
                %bank_rf_out
                Bank_byscen(:, oa, scen_ii) = bank_out;
                Emiss_byscen(:, oa, scen_ii) = emiss_out;
                MF_byscen(:, oa, scen_ii) = mf_out;
                
            end
    
             MF_scen = OAmolefractions;
             MF_scen(:,2:5) = squeeze(MF_byscen(:, 1:4, scen_ii));
        elseif scen_ii == 2
            % Scenario 3: All lifetimes are the same in foreward simulation. 
            title_str= 'Scenario 2: Lifetime correction';
            Yend = [2006, 2010, 2014, 2018];
            clear MF_sim
            for oa = 1:4
    
                    mf_obs_yr = Yend(oa);  % Edit this for each scenario
                    LT = lt_2022;          % Edit this for each scenario
    
                    Gg2ppt = ppt_per_kg(oa)*(LT*(1-exp(-1/LT)));
                    mf_in = OAmolefractions(:, oa+1);
                    
                    start_yr =1950; 
                    end_yr = 2100;
                    
                    bank_ref_yr = Bank_refyear(oa);
                    bank_ref = Bank_ref_val(oa);
                    bank_rf = Bank_RF(oa);
                    
                    prod_in = Prod(:,oa+1);
                    prod_yrs = 1950:2100;
                    inferRF = 1; 

                    [bank_out, emiss_out, mf_out, bank_rf_out] = bank_emiss_mf(LT, Gg2ppt, mf_in, start_yr, end_yr, mf_obs_yr, bank_ref_yr, bank_ref, bank_rf, prod_in, prod_yrs, inferRF, Bank2015, oa);
    
                    Bank_byscen(:, oa, scen_ii) = bank_out;
                    Emiss_byscen(:, oa, scen_ii) = emiss_out;
                    MF_byscen(:, oa, scen_ii) = mf_out;
            end
            
    
             MF_scen = OAmolefractions;
             MF_scen(:,2:5) = squeeze(MF_byscen(:, 1:4, scen_ii));
    
        elseif scen_ii == 3
                title_str= 'Scenario 3: Mf to 2022';
                clear MF_sim
                for oa = 1:4%5
                    
                    mf_obs_yr = 2022-1;  % Edit this for each scenario
                    LT = lt_2022;          % Edit this for each scenario
    
                    Gg2ppt = ppt_per_kg(oa)*(LT*(1-exp(-1/LT)));
                    mf_in = OAmolefractions(:, 6);
    
                    
                    start_yr =1950; 
                    end_yr = 2100;
                    
                    bank_ref_yr = Bank_refyear(oa);
                    bank_ref = Bank_ref_val(oa);
                    bank_rf = Bank_RF(oa);

                    if oa == 1 
                        inferRF = 0; 
                    else
                        inferRF = 1; 
                    end
                    
                    prod_in = Prod(:,oa+1);
                    prod_yrs = 1950:2100;
                    

                    [bank_out, emiss_out, mf_out, bank_rf_out] = bank_emiss_mf(LT, Gg2ppt, mf_in, start_yr, end_yr, mf_obs_yr, bank_ref_yr, bank_ref, bank_rf, prod_in, prod_yrs, inferRF, Bank2015, oa);
    
                    Bank_byscen(:, oa, scen_ii) = bank_out;
                    Emiss_byscen(:, oa, scen_ii) = emiss_out;
                    MF_byscen(:, oa, scen_ii) = mf_out;
    
                    MF_scen = OAmolefractions;
                    MF_scen(:,2:5) = squeeze(MF_byscen(:, 1:4, scen_ii));
    
                    
                end
    
        elseif scen_ii == 4
        % Scenario 4: Emissions following publication follow 2022 emissions
        % assumptions
            title_str= 'Scenario 4: Emissions = 2022 OA';
            Yend = [2006, 2010, 2014, 2018];
            clear MF_sim
            for oa = 1:4%5
                    y1 = 1950;
                    y2 = 2022;
                    ind1 = find(OAmolefractions(:,1) == y1);
                    ind2 = find(OAmolefractions(:,1) == y2);
                    Ggtoppt = ppt_per_kg(oa)*(lt_2022*(1-exp(-1/lt_2022)));
                    MF_sim(ind1:ind2,oa) = OAmolefractions(ind1:ind2, 6);
                    
                for ii = 1:2100-y2
                    Inferred_emiss = (OAmolefractions(ind2+ii,6)-OAmolefractions(ind2+ii-1,6)*exp(-1/lt_2022))/Ggtoppt;
                    MF_sim(ind2 + ii, oa) = exp(-1/lt_2022)*MF_sim(ind2+ii-1,oa) + Ggtoppt*Inferred_emiss;
                end
            end
            
            MF_scen = OAmolefractions;
            MF_scen(:,2:5) = MF_sim;
            MF_byscen(:, 1:4, scen_ii) = MF_scen(:,2:5);
        end
    
    
    ind = find(MF_scen(:,1) == 1980);
    mf_1980 = MF_scen(ind,6);
        if halon1301
        else
            for oa = 1:5
                return_date(oa) = MF_scen(40,1) + find( MF_scen(41:end,oa+1) < mf_1980, 1 );
            end
            
            subplot(2,2,scen_ii);
            oa = 1;     p1 = plot(MF_scen(:,1),MF_scen(:,oa+1),'LineWidth', 2,'Color', LineCols(oa,:)); hold on;
            oa = 2;     p2 = plot(MF_scen(:,1),MF_scen(:,oa+1),'LineWidth', 2,'Color', LineCols(oa,:)); hold on;
            oa = 3;     p3 = plot(MF_scen(:,1),MF_scen(:,oa+1),'LineWidth', 2,'Color', LineCols(oa,:)); hold on;
            oa = 4;     p4 = plot(MF_scen(:,1),MF_scen(:,oa+1),'LineWidth', 2,'Color', LineCols(oa,:)); hold on;
            oa = 5;     p5 = plot(MF_scen(:,1),MF_scen(:,oa+1),'LineWidth', 2,'Color', LineCols(oa,:)); hold on;
            
            plot([1980,return_date(5)],[mf_1980, mf_1980],'--k','LineWidth',2);hold on;
            plot([1980,1980],[0, mf_1980],'--k','LineWidth',2); hold on; 
            for oa = 1:5
                hold on; 
                plot([return_date(oa), return_date(oa)],[0, mf_1980],'--','LineWidth',2, 'Color', LineCols(oa,:));
            end
            
            legend([p1, p2, p3, p4, p5], '2006 OA','2010 OA','2014 OA','2018 OA','2020 OA');
            xlim([1970,2080])
            ylabel('[CFC-11] pmol mol^{-1}')
            title(title_str)
            
            Return_Date2(scen_ii,:) = return_date;
        end
    end
     
    if halon1301
    else
        Return_Date2
    
        
        figure_width = 18; % in inches
        figure_height = 12; % in inches
        screen_ppi = 72; 
        
        screen_figure_width = round(figure_width*screen_ppi); % in pixels
        screen_figure_height = round(figure_height*screen_ppi); % in pixels
        set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [figure_width figure_height]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);
        
        FigName = scen_ii;
        Figstr = strcat(mole_str, 'IncrementalScenarios.pdf'); 
        print(gcf, '-dpdf', Figstr);
        
        fighandle = figure
        plot([Return_Date2(1,5), Return_Date2(1,5)], [0,5], '--k','LineWidth',2); 
        hold on; 
        p1 = plot(Return_Date2(1,1:4), [4,3,2,1], '.','MarkerSize',20, 'Color', LineCols(1,:)); hold on; 
        p2 = plot(Return_Date2(2,1:4), [4,3,2,1], '.','MarkerSize',20, 'Color', LineCols(2,:)); hold on; 
        p3 = plot(Return_Date2(3,1:4), [4,3,2,1], '.','MarkerSize',20, 'Color', LineCols(3,:)); hold on; 
        p4 = plot(Return_Date2(4,1:4), [4,3,2,1], '.','MarkerSize',20, 'Color', LineCols(4,:));
        legend([p1,p2,p3,p4],'initial OA','LT correction','MF correction','Banks correction')
        yticklabels('')
        xlim([min(min(Return_Date2))-2,max(max(Return_Date2))+2])
        xlabel('Return date'); box off
    
        figure_width = 18; % in inches
        figure_height = 6; % in inches
        screen_ppi = 72; 
        
        screen_figure_width = round(figure_width*screen_ppi); % in pixels
        screen_figure_height = round(figure_height*screen_ppi); % in pixels
        set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [figure_width figure_height]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);
        
        FigName = scen_ii;
        Figstr = strcat(mole_str,'IncrementalScenarios_ReturnDate.pdf'); 
        print(gcf, '-dpdf', Figstr);
    
    end
    
    if gas_ii == 1
        cfc11_banks = Bank_byscen;
        cfc11_mf = MF_byscen;
        cfc11_emiss = Emiss_byscen;
    elseif gas_ii == 2
        cfc12_banks = Bank_byscen; 
    elseif gas_ii == 3
        halon1301_banks = Bank_byscen;
    end
    
    fighandle = figure; 
    for oa = 1:4
        subplot(4,3,3*oa - 2)
        
        p1 = plot([1950:2100],0.001*squeeze(Bank_byscen(:,oa,1))','LineWidth',2,'Color', LineCols(1,:)); hold on;
        p2 = plot([1950:2100],0.001*squeeze(Bank_byscen(:,oa,2))','LineWidth',2,'Color', LineCols(2,:)); hold on;
        p3 = plot([1950:2100],0.001*squeeze(Bank_byscen(:,oa,3))','LineWidth',2,'Color', LineCols(3,:));
        title(strcat(num2str(2002+4*oa), ' OA Banks')); ylabel('Thousand Gg')
        xlim([1990,2050])
        
        subplot(4,3,3*oa-1)
        p1 = plot([1950:2100],squeeze(Emiss_byscen(:,oa,1))','LineWidth',2,'Color', LineCols(1,:)); hold on;
        p2 = plot([1950:2100],squeeze(Emiss_byscen(:,oa,2))','LineWidth',2,'Color', LineCols(2,:)); hold on;
        p3 = plot([1950:2100],squeeze(Emiss_byscen(:,oa,3))','LineWidth',2,'Color', LineCols(3,:));
        xlim([1990,2050]); ylabel('Thousand Gg')
        title(strcat(num2str(2002+4*oa), ' OA Emissions'))
        
        subplot(4,3,3*oa)
        p1 = plot([1950:2100],squeeze(MF_byscen(:,oa,1))','LineWidth',2,'Color', LineCols(1,:)); hold on;
        p2 = plot([1950:2100],squeeze(MF_byscen(:,oa,2))','LineWidth',2,'Color', LineCols(2,:)); hold on;
        p3 = plot([1950:2100],squeeze(MF_byscen(:,oa,3))','LineWidth',2,'Color', LineCols(3,:)); hold on; 
        p4 = plot([1950:2100], MF_scen(:,6), 'k');
        
        xlim([2040,2080]); ylabel('MF [pmol/mol]')
        title(strcat(num2str(2002+4*oa), ' OA MF'))
        legend([p1, p2, p3], 'Initial OA','LT correction','MF correction')
    end
    
    figure_width = 18; % in inches
    figure_height = 12; % in inches
    screen_ppi = 72; 
    
    screen_figure_width = round(figure_width*screen_ppi); % in pixels
    screen_figure_height = round(figure_height*screen_ppi); % in pixels
    set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [figure_width figure_height]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);
    
    FigName = scen_ii;
    Figstr = strcat(mole_str,'IncrementalScenarios_Bank_Emiss_MF.pdf'); 
    print(gcf, '-dpdf', Figstr);

    for oa = 1:4
        for scen_ii = 1:4
            mole_frac{4*(scen_ii-1)+oa}(:,gas_ii) = squeeze(MF_byscen(:, oa, scen_ii));
        end
    end
end




ozone_assessmentyear = {'OA2006','OA2010', 'OA2014','OA2018','OA2022'};
scenario_name = {'InitialOA','LTupdates','MFupdates','BKupdates'};
cd Scenarios_jan2024
for oa = 1:4
    for scen_ii = 1:4
        table_str = strcat(ozone_assessmentyear(oa),scenario_name(scen_ii),'.txt');
        str = string(table_str);
        cfc11 = mole_frac{4*(scen_ii-1)+oa}(:,1);
        cfc12 = mole_frac{4*(scen_ii-1)+oa}(:,2);
        halon1301 = mole_frac{4*(scen_ii-1)+oa}(:,3);
        year = [1950:2100]';
        T = table(year,cfc11,cfc12,halon1301);
        writetable(T, str);
    end
end
cd ..


%%  CCL4 analysis
load('ccl4_OA_data.mat')
lt = [26, 26, 26, 32, 30];
lt_2022 = 30; 
ppt_per_kg = 10^(-5)*[3.9505, 3.9505, 3.9505, 3.9505, 3.9505];
Yend = [2006, 2010, 2014, 2018];
gas_ii = 4; 

% 2006 WMO: Scenario 1; Initial Assessment
oa = 1; 
scen_ii = 1; 
lifetime = lt(oa);
Yend_scen = 2005;
[tmp] = oa_ccl4_methods(lifetime, ccl4_mf(:,oa+1), ccl4_mf(:,1), Yend_scen, oa);
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = tmp; 

% 2006 WMO: Scenario 2; Lifetime effect
scen_ii = 2;
lifetime = lt_2022;
[tmp] = oa_ccl4_methods(lifetime, ccl4_mf(:,oa+1), ccl4_mf(:,1), Yend_scen, oa);
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = tmp; 

% 2006 WMO: Scenario 3: Mole fraction effect
scen_ii = 3;
lifetime = lt_2022;
Yend_scen = 2022;
[tmp] = oa_ccl4_methods(lifetime, ccl4_mf(:,end), ccl4_mf(:,1), Yend_scen, oa);
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = tmp; 

% 2006 WMO: Scenario 4: Feedstock emission assumption
scen_ii = 4; 
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = ccl4_mf(:,end); 

% 2010 WMO: Scenario 1
oa = 2; 
scen_ii = 1; 
lifetime = lt(oa);
Yend_scen = 2009;
[tmp] = oa_ccl4_methods(lifetime, ccl4_mf(:,oa+1), ccl4_mf(:,1), Yend_scen, oa);
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = tmp; 

% 2010 WMO: Scenario 2
scen_ii = 2;
lifetime = lt_2022;
[tmp] = oa_ccl4_methods(lifetime, ccl4_mf(:,oa+1), ccl4_mf(:,1), Yend_scen, oa);
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = tmp; 

% 2010 WMO: Scenario 3
scen_ii = 3;
lifetime = lt_2022;
Yend_scen = 2022;
[tmp] = oa_ccl4_methods(lifetime, ccl4_mf(:,end), ccl4_mf(:,1), Yend_scen, oa);
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = tmp; 

% 2010 WMO: Scenario 4
scen_ii = 4; % Instead of Bank update this is the feedstock emission update
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = ccl4_mf(:,end); 

% 2014 WMO: Scenario 1
oa = 3; 
scen_ii = 1; 
lifetime = lt(oa);
Yend_scen = 2013;
[tmp] = oa_ccl4_methods(lifetime, ccl4_mf(:,oa+1), ccl4_mf(:,1), Yend_scen, oa);
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = tmp; 

% 2014 WMO: Scenario 2
scen_ii = 2;
lifetime = lt_2022;
[tmp] = oa_ccl4_methods(lifetime, ccl4_mf(:,oa+1), ccl4_mf(:,1), Yend_scen, oa);
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = tmp; 

% 2014 WMO: Scenario 3
scen_ii = 3;
lifetime = lt_2022;
Yend_scen = 2022;
[tmp] = oa_ccl4_methods(lifetime, ccl4_mf(:,end), ccl4_mf(:,1), Yend_scen, oa);
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = tmp; 

% 2014 WMO: Scenario 4
scen_ii = 4; % Instead of Bank update this is the feedstock emission update
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = ccl4_mf(:,end); 

% 2018 WMO: Scenario 1
oa = 4; % 2018 SAOD
scen_ii = 1; 
lifetime = lt(oa);
Yend_scen = 2017;
[tmp] = oa_ccl4_methods(lifetime, ccl4_mf(:,oa+1), ccl4_mf(:,1), Yend_scen, oa);
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = tmp; 

% 2018 WMO: Scenario 2
scen_ii = 2;
lifetime = lt_2022;
[tmp] = oa_ccl4_methods(lifetime, ccl4_mf(:,oa+1), ccl4_mf(:,1), Yend_scen, oa);
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = tmp; 

% 2018 WMO: Scenario 3
scen_ii = 3;
lifetime = lt_2022;
Yend_scen = 2022;
[tmp] = oa_ccl4_methods(lifetime, ccl4_mf(:,end), ccl4_mf(:,1), Yend_scen, oa);
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = tmp; 

% 2018 WMO: Scenario 4
scen_ii = 4; % Instead of Bank update this is the feedstock emission update
mole_frac{4*(scen_ii-1)+oa}(2:end,gas_ii) = ccl4_mf(:,end); 


ozone_assessmentyear = {'ccl4_2006','ccl4_2010', 'ccl4_2014','ccl4_2018'};
scenario_name = {'_originalSAOD_lifetime_molefraction_FullUpdate'};
cd Scenarios_jan2024
for oa = 1:4
    table_str = strcat(ozone_assessmentyear(oa),scenario_name,'.txt');
    str = string(table_str);
    scen_ii = 1;
    update1 = mole_frac{4*(scen_ii-1)+oa}(:, 4);
    scen_ii = 2;
    update2 = mole_frac{4*(scen_ii-1)+oa}(:, 4);
    scen_ii = 3;
    update3 = mole_frac{4*(scen_ii-1)+oa}(:, 4);
    scen_ii = 4;
    update4 = mole_frac{4*(scen_ii-1)+oa}(:, 4);
    year = [1950:2100]';
    T = table(year, update1, update2, update3, update4);
    writetable(T, str);

end
cd ..

end
%%
load('ccl4_OA_data.mat')
ppt_per_kg = 10^(3)*(1/25313.18968);
lt = 26; 
Gg2ppt = ppt_per_kg*(lt*(1-exp(-1/lt)));

Emiss_2006 = (1/Gg2ppt)*(ccl4_mf(2:end, 2) - exp(-1/lt)*ccl4_mf(1:end-1, 2));
Emiss_2010 = (1/Gg2ppt)*(ccl4_mf(2:end, 3) - exp(-1/lt)*ccl4_mf(1:end-1, 3));
Emiss_2014 = (1/Gg2ppt)*(ccl4_mf(2:end, 4) - exp(-1/lt)*ccl4_mf(1:end-1, 4));
lt = 32;
Emiss_2018 = (1/Gg2ppt)*(ccl4_mf(2:end, 5) - exp(-1/lt)*ccl4_mf(1:end-1, 5));
lt = 30;
Emiss_2022 = (1/Gg2ppt)*(ccl4_mf(2:end, 6) - exp(-1/lt)*ccl4_mf(1:end-1, 6));
     
fighandle = figure; 
ind1 = find(ccl4_mf(:,1) == 2005);
p1 = plot(ccl4_mf(2:ind1,1), Emiss_2006(1:ind1-1), 'Color','k',  'LineWidth',2)
hold on;
plot(ccl4_mf(ind1:end,1), Emiss_2006(ind1-1:end), '--', 'Color','k',  'LineWidth',2)
hold on; 
plot(ccl4_mf(ind1,1), Emiss_2006(ind1-1), '.', 'MarkerSize', 30, 'Color','k')
xlim([1995,2100])
legend([p1],'2006 OA');
ylabel('CCl4 Emissions [Gg yr^{-1}]')

figure_width = 10; % in inches
figure_height = 5; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);
set(gca, 'FontSize', 12)

%FigName = scen_ii;
Figstr = strcat('CCl4_OA_2006.pdf'); 
print(gcf, '-dpdf', Figstr);

ind1 = find(ccl4_mf(:,1) == 2009);
hold on;
p2 = plot(ccl4_mf(2:ind1,1), Emiss_2010(1:ind1-1), 'r', 'LineWidth',2)
hold on;
plot(ccl4_mf(ind1:end,1), Emiss_2010(ind1-1:end), 'r--', 'LineWidth',2)
xlim([1995,2100])
hold on; 
plot(ccl4_mf(ind1,1), Emiss_2010(ind1-1), '.', 'MarkerSize', 30, 'Color','r')
legend([p1, p2],'2006 OA','2010 OA')
Figstr = strcat('CCl4_OA_2010.pdf'); 
print(gcf, '-dpdf', Figstr);
       
ind1 = find(ccl4_mf(:,1) == 2013);
hold on; 
p3 = plot(ccl4_mf(2:ind1,1), Emiss_2014(1:ind1-1), 'b', 'LineWidth',2)
hold on;
plot(ccl4_mf(ind1:end,1), Emiss_2014(ind1-1:end), 'b--', 'LineWidth',2)
xlim([1995,2100])
hold on; 
plot(ccl4_mf(ind1,1), Emiss_2014(ind1-1), '.', 'MarkerSize', 30, 'Color','b')
legend([p1, p2, p3],'2006 OA','2010 OA', '2014 OA')
Figstr = strcat('CCl4_OA_2014.pdf'); 
print(gcf, '-dpdf', Figstr);

ind1 = find(ccl4_mf(:,1) == 2017);
hold on; 
p4 = plot(ccl4_mf(2:ind1,1), Emiss_2018(1:ind1-1), 'Color',"#FF8800", 'LineWidth',2)
hold on;
plot(ccl4_mf(ind1:end,1), Emiss_2018(ind1-1:end), '--', 'Color',"#FF8800",  'LineWidth',2)
hold on; 
plot(ccl4_mf(ind1,1), Emiss_2018(ind1-1), '.', 'MarkerSize', 30, 'Color',"#FF8800")
xlim([1995,2100])
legend([p1, p2, p3, p4],'2006 OA','2010 OA', '2014 OA', '2018 OA')
Figstr = strcat('CCl4_OA_2018.pdf'); 
print(gcf, '-dpdf', Figstr);

ind1 = find(ccl4_mf(:,1) == 2021);
hold on;
p5 = plot(ccl4_mf(2:ind1,1), Emiss_2022(1:ind1-1), 'Color',1/255*[148 0 211], 'LineWidth',2)
hold on;
plot(ccl4_mf(ind1:end,1), Emiss_2022(ind1-1:end), '--','Color',1/255*[148 0 211], 'LineWidth',2)
hold on; 
plot(ccl4_mf(ind1,1), Emiss_2022(ind1-1), '.', 'MarkerSize', 30, 'Color',1/255*[148 0 211])

xlim([1995,2100])

legend([p1, p2, p3, p4, p5],'WMO 2006','WMO 2010', 'WMO 2014', 'WMO 2018', 'WMO 2022')
Figstr = strcat('CCl4_OA_2022.pdf'); 
print(gcf, '-dpdf', Figstr);

%%
% 2022 Assessment
cd('/Users/meganlickley/Desktop/2022 Research/Kicking the Can/eesc_take2')
cd('2022')
file_name = strcat('eesc_3yr_old_air_engelFRFs.wmo2022_table7A-1.from_john_daniel.april_2023.dat');
            
%fidi = fopen(file_name,'rt');

a=readcell(file_name,"Delimiter"," ");
%fclose(fidi);
for kk = 8:150*12+7
    
    data_col = 1; 
    mm = 0; 
    while data_col <4
        if ismissing(a{kk, data_col + mm + 1})
            mm = mm+1;
        else
            data_col = data_col + 1;
        end
    end
    tmp(kk-7,1) = a{kk, 1};
    tmp(kk-7,2) = a{kk, data_col + mm};
end

indx1 = find(tmp(:,1)>1980, 1); 
indx2 = find(tmp(:,1)>1990,1); % Finding the first value after 1990 that dips below 1980 levels
indx3 = find(tmp(indx2+1:end,2)<tmp(indx1,2), 1)+indx2;
eesc_return_date(5,1) = tmp(indx3,1);
eesc_1980val(5,1) = tmp(indx1,2);
eesc_val(5,1)     = tmp(indx3,2);

if run_make_final_figures
    % in terminal window, download the files: wget -r -np -nH --cut-dirs=1 -R index.html https://www2.atmos.umd.edu/~rjs/eesc/    
    OA_dates = {'2006'; '2010'; '2014'; '2018'}
    
    scen = {'InitialOA','LTupdates_','MFupdates_','BKupdates_'}
    for oa = 1:4
        ll = 1;
        for scen_ii = 1:4
            cd('/Users/meganlickley/Desktop/2022 Research/Kicking the Can/eesc_jan2024')
            dir_name = strcat(OA_dates{oa},'/output_files');
            cd(dir_name)
            if scen_ii == 1 
                chem_num = 1; 
                Num_gases = 1;
            elseif scen_ii == 4
                chem_num = 2; 
                Num_gases = 7; 
            else
                chem_num = 2; 
                Num_gases = 4; 
            end
    
            for jj = chem_num:Num_gases
                if jj == 1 
                    gas_str = '';
                elseif jj == 2
                    gas_str = 'cfc11';
                elseif jj == 3
                    gas_str = 'cfc11_cfc12';
                elseif jj == 4
                    gas_str = 'cfc11_cfc12_halon1301';
                elseif jj == 5
                    gas_str = 'ccl4_update2';
                elseif jj == 6
                    gas_str = 'ccl4_update3';
                elseif jj == 7
                    gas_str = 'ccl4_update4';
                end
                file_name = strcat('eesc_3yr_old_air_engelFRFs.wmo',OA_dates{oa},'_megan_',scen{scen_ii},gas_str,'.dat');
                if jj > 4
                    file_name = strcat('eesc_3yr_old_air_engelFRFs.wmo',OA_dates{oa},'_megan_',gas_str,'.dat');
                end
                %fidi = fopen(file_name,'rt');
               
                a=readcell(file_name,"Delimiter"," ");
                %fclose(fidi);
                for kk = 8:150*12+7
                    
                    data_col = 1; 
                    mm = 0; 
                    while data_col <4
                        if ismissing(a{kk, data_col + mm + 1})
                            mm = mm+1;
                        else
                            data_col = data_col + 1;
                        end
                    end
                    tmp(kk-7,1) = a{kk, 1};
                    tmp(kk-7,2) = a{kk, data_col + mm};
                end

                indx1 = find(tmp(:,1)>1980, 1); 
                indx2 = find(tmp(:,1)>1990, 1);
                eesc_1980val(oa,ll) = 0.5*(tmp(indx1 - 1,2) + tmp(indx1,2));

                indx3 = find(tmp(indx2+1:end,2) < eesc_1980val(oa,ll),1)+indx2;
                eesc_return_date(oa,ll) = interp1([tmp(indx3-1,2), tmp(indx3,2)], [tmp(indx3-1,1), tmp(indx3,1)], eesc_1980val(oa,ll));
                eesc_val(oa,ll) = interp1([tmp(indx3-1,1), tmp(indx3,1)], [tmp(indx3-1,2), tmp(indx3,2)], eesc_return_date(oa,ll));
                ll = ll +1;
            end
        end
    end
end

eesc_return_date(5,:) = eesc_return_date(5,1);

%% Reading in crossing date directly
%cd 1980_crossing/
% Reading in 2022 return dates
dirname = '/Users/meganlickley/Desktop/2022 Research/Kicking the Can/eesc_take2/2022/1980_crossing_corrected';
cd(dirname)
file_name = strcat('eesc_3yr_old_air_engelFRFs.wmo2022_table7A-1.from_john_daniel.april_2023.1980_crossing.dat');
           

a=readcell(file_name,"Delimiter"," ");
eesc_return_date(5,1:13) = a{5,1};
cd('/Users/meganlickley/Desktop/2022 Research/Kicking the Can/eesc_jan2024') 
%cd eesc_jan2024
yr_opts = {'2006','2010','2014','2018'};
updates = { '_InitialOA',
            '_LTupdates_cfc11',
            '_LTupdates_cfc11_cfc12',
            '_LTupdates_cfc11_cfc12_halon1301',
            '_MFupdates_cfc11',
            '_MFupdates_cfc11_cfc12',
            '_MFupdates_cfc11_cfc12_halon1301',
            '_BKupdates_cfc11',
            '_BKupdates_cfc11_cfc12',
            '_BKupdates_cfc11_cfc12_halon1301',
            '_ccl4_update2',
            '_ccl4_update3',
            '_ccl4_update4',
            }
for oa = 1:4
    cd('/Users/meganlickley/Desktop/2022 Research/Kicking the Can') 
    cd eesc_jan2024
     dir_name = strcat(yr_opts{oa},'/output_files/1980_crossing');
     cd(dir_name)

     file_str1 = strcat('eesc_3yr_old_air_engelFRFs.wmo',yr_opts{oa},'_megan');
     for ii = 1:13
        file_str2 = strcat(updates{ii},'.1980_crossing.dat');
        file_name = strcat(file_str1, file_str2)
        a=readcell(file_name,"Delimiter"," ");
        eesc_return_date(oa,ii) = a{5,1};
     end
end

eesc_return_date(:,14) = eesc_return_date(5,1); % final update to current assessment
%%
fighandle = figure
x_tmp = [1,4];
curve1 = [2050, 2050];
curve2 = [2080, 2080];

x2 = [x_tmp, fliplr(x_tmp)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.95,0.95,0.95],'EdgeColor','none');

x_tmp = [4,7]
x2 = [x_tmp, fliplr(x_tmp)];
inBetween = [curve1, fliplr(curve2)];
hold on;
fill(x2, inBetween, [0.85,0.85,0.85],'EdgeColor','none');


x_tmp = [7,10]
x2 = [x_tmp, fliplr(x_tmp)];
inBetween = [curve1, fliplr(curve2)];
hold on;
fill(x2, inBetween, [0.8,0.8,0.8],'EdgeColor','none');

x_tmp = [10,11]
x2 = [x_tmp, fliplr(x_tmp)];
inBetween = [curve1, fliplr(curve2)];
hold on;
fill(x2, inBetween, [0.75, 0.75, 0.75],'EdgeColor','none');

x_tmp = [11,12]
x2 = [x_tmp, fliplr(x_tmp)];
inBetween = [curve1, fliplr(curve2)];
hold on;
fill(x2, inBetween, [0.7, 0.7, 0.7],'EdgeColor','none');

hold on;
p1 = plot(eesc_return_date(1,:)','-o','LineWidth', 2,'Color',[1, 0 0]);hold on; 
p2 = plot(eesc_return_date(2,:)','-o','LineWidth', 2,'Color',[0.7, 0 0.3]);hold on; 
p3 = plot(eesc_return_date(3,:)','-o','LineWidth', 2,'Color',[0.4, 0 0.6]);hold on; 
p4 = plot(eesc_return_date(4,:)','-o','LineWidth', 2,'Color',[0.3, 0 0.9]);hold on; 
p5 = plot(eesc_return_date(5,:)','--k','LineWidth', 2);
%legend([p1, p2, p3, p4, p5],'WMO 2006','WMO 2010','WMO 2014','WMO 2018','WMO 2022')
xlim([1,12])
ylim([2052,2070])
ylabel('EESC return date')
set(gca, 'FontSize',16)
figure_width = 8; % in inches
figure_height = 6; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = scen_ii;
Figstr = strcat('EESC_Line_Plot.pdf'); 
print(gcf, '-dpdf', Figstr);

%%  Summary EESC_return_date
% 2006 Color = 'b'
% 2010 Color = "#FF8800" 
% 2014 Color = [0 0.6 0.3]
% 2018 Color = 'r'
% 2022 Color = 1/255*[148 0 211]

LineCols = [0,      0.4470, 0.7410;
            0.8500, 0.3250, 0.0980;
            0.9290, 0.6940, 0.1250;
            0.4940, 0.1840, 0.5560;
            0.4660, 0.6740, 0.1880];

fighandle = figure
y0 = [2048.4, 2046.1, 2047.7, 2049.5];

for oa = 1:4
    
    inc_diff = eesc_return_date(oa, 2:end) - eesc_return_date(oa, 1:end-1);
    y1 = eesc_return_date(oa, 1);
    yfinal = eesc_return_date(5, 1);

     curve1 = [5 - oa - 0.2,  5 - oa - 0.2, 5 - oa]; 
     curve2 = [5 - oa + 0.2,  5 - oa + 0.2, 5 - oa]; 
     x_tmp = [y0(oa),  y1-0.25, y1];

     hold on; 
     x2 = [x_tmp, fliplr(x_tmp)];
     inBetween = [curve1, fliplr(curve2)];
     p6 = fill(x2, inBetween, [1,1,1],'EdgeColor','k', 'LineWidth',2);


    hold on;
    if oa == 1
        plot( [y1, y1], [0, 5 - oa], '--b')
    elseif oa == 2
        plot( [y1, y1], [0, 5 - oa], '--','Color',"#FF8800" );
    elseif oa == 3
        plot( [y1, y1], [0, 5 - oa], '--','Color',[0 0.6 0.3]);
    elseif oa == 4
        plot( [y1, y1], [0, 5 - oa], '--','Color','r');
    end
    
    for jj = 1:3
        for ii = 1:3
            last_pos_val = 3 - find(fliplr(inc_diff(3*(jj-1)+1:3*(jj-1)+3))>0, 1) + 1;
            first_neg_val = find(inc_diff(3*(jj-1)+1:3*(jj-1)+3)<0, 1);
            if inc_diff(3*(jj-1)+ii)>0
                
                if ii == last_pos_val 
                    curve1 = [5 - oa - 0.2,  5 - oa - 0.2, 5 - oa]; 
                    curve2 = [5 - oa + 0.2,  5 - oa + 0.2, 5 - oa]; 
                    x_tmp = [y1, max(y1, y1 + inc_diff(ii+3*(jj-1)) - 0.25), y1+inc_diff(3*(jj-1)+ii)];
                    

                else
                    curve1 = [5 - oa - 0.2,  5 - oa - 0.2]; 
                    curve2 = [5 - oa + 0.2,  5 - oa + 0.2]; 
                    x_tmp = [y1, y1+inc_diff(3*(jj-1)+ii)];

                end                
                y1 = y1+inc_diff(3*(jj-1)+ii);

                hold on; 
                x2 = [x_tmp, fliplr(x_tmp)];
                inBetween = [curve1, fliplr(curve2)];
                fill(x2, inBetween, LineCols(ii,:),'EdgeColor','none');

                if oa == 1 && jj == 1
                    if ii == 1
                        hold on;
                        p1 = fill(x2, inBetween, LineCols(ii,:),'EdgeColor','none');
                    elseif ii == 2
                        hold on;
                        p2 = fill(x2, inBetween, LineCols(ii,:),'EdgeColor','none');
                    elseif ii == 3
                        hold on;
                        p3 = fill(x2, inBetween, LineCols(ii,:),'EdgeColor','none');
                    end
                end

            elseif inc_diff(3*(jj-1)+ii)<0
                if ii == first_neg_val 
                    curve1 = [5 - oa, 5 - oa - 0.2,  5 - oa - 0.2]-0.4; 
                    curve2 = [5 - oa, 5 - oa + 0.2,  5 - oa + 0.2]-0.4; 
                    x_tmp = [yfinal, min(yfinal+0.25, yfinal+abs(inc_diff(3*(jj-1)+ii))), yfinal+abs(inc_diff(3*(jj-1)+ii))];
 

                else
                    curve1 = [5 - oa - 0.2,  5 - oa - 0.2]-0.4; 
                    curve2 = [5 - oa + 0.2,  5 - oa + 0.2]-0.4; 
                    x_tmp = [yfinal, yfinal+abs(inc_diff(3*(jj-1)+ii))];

                end
                yfinal = yfinal+abs(inc_diff(3*(jj-1)+ii));

                 hold on; 
                x2 = [x_tmp, fliplr(x_tmp)];
                inBetween = [curve1, fliplr(curve2)];
                fill(x2, inBetween, LineCols(ii,:),'EdgeColor','none');
                    
            end

        end
    end
    
    for ii = 10:13
        if inc_diff(ii)>0
            curve1 = [5 - oa - 0.2,  5 - oa - 0.2, 5 - oa]; 
            curve2 = [5 - oa + 0.2,  5 - oa + 0.2, 5 - oa]; 
            x_tmp = [y1, max(y1, y1 + inc_diff(ii) - 0.25), y1+inc_diff(ii)];
            y1 = y1+inc_diff(ii);
            
            hold on; 
            x2 = [x_tmp, fliplr(x_tmp)];
            inBetween = [curve1, fliplr(curve2)];
            if ii == 10
                fill(x2, inBetween, LineCols(ii-6,:),'EdgeColor','none','FaceAlpha',0.9);
            elseif ii == 11
                if oa == 1
                p4 = fill(x2, inBetween, LineCols(ii-7,:),'EdgeColor','none','FaceAlpha',0.6);
                else
                fill(x2, inBetween, LineCols(ii-7,:),'EdgeColor','none','FaceAlpha',0.6);
                end
            elseif ii == 12
                fill(x2, inBetween, LineCols(ii-8,:),'EdgeColor','none','FaceAlpha',0.3);
            elseif ii == 13
                if oa == 1
                p5 = fill(x2, inBetween, LineCols(ii-8,:),'EdgeColor','none');
                else
                    fill(x2, inBetween, LineCols(ii-8,:),'EdgeColor','none');
                end
            end
                         
%              if oa == 1
%                  if ii == 12
%                  hold on;
%                  p4 = fill(x2, inBetween, LineCols(ii-8,:),'EdgeColor','none','FaceAlpha',0.3);
%                  elseif ii == 13
%                  hold on;
%                  p5 = fill(x2, inBetween, LineCols(ii-8,:),'EdgeColor','none');  
%                  end
%              end
        elseif inc_diff(ii) < 0 
    
             curve1 = [5 - oa, 5 - oa - 0.2,  5 - oa - 0.2] - 0.4; 
             curve2 = [5 - oa, 5 - oa + 0.2,  5 - oa + 0.2] - 0.4; 
             x_tmp = [yfinal, min(yfinal+0.25, yfinal+abs(inc_diff(ii))), yfinal+abs(inc_diff(ii))];
             yfinal = yfinal+abs(inc_diff(ii));
             
             hold on; 
             x2 = [x_tmp, fliplr(x_tmp)];
             inBetween = [curve1, fliplr(curve2)];
             if ii == 10
                fill(x2, inBetween, LineCols(ii-6,:),'EdgeColor','none','FaceAlpha',0.9);
            elseif ii == 11
                fill(x2, inBetween, LineCols(ii-7,:),'EdgeColor','none','FaceAlpha',0.6);
            elseif ii == 12
                fill(x2, inBetween, LineCols(ii-8,:),'EdgeColor','none','FaceAlpha',0.3);
             elseif ii == 13
                fill(x2, inBetween, LineCols(ii-8,:),'EdgeColor','none');  
            end
        end



    end
end



    hold on; plot( [eesc_return_date(5, 1), eesc_return_date(5, 1)], [0, 5], '--k')
    xlabel('EESC return date')
    box on; 

 
legend([p6, p1, p2, p3, p4, p5],'EESC Formulation', 'CFC-11', 'CFC-12', 'Halon-1301', 'CCl4', 'Other12')
    xlim([2040,2072])
    figure_width = 16; % in inches

    figure_height = 6; % in inches

    screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = scen_ii;
Figstr = strcat('EESC_Summary_Arrows_ccl4_contributions.pdf'); 
print(gcf, '-dpdf', Figstr);
                
%%






fighandle = figure; 
x_tmp = 1:4;
curve1 = 2050*ones(size(x_tmp)); 
curve2 = 2070*ones(size(x_tmp)); 
x2 = [x_tmp, fliplr(x_tmp)];
hold on; 
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.9,0.9,0.9],'EdgeColor','none');

hold on; 

x_tmp = 7:10;
x2 = [x_tmp, fliplr(x_tmp)];
hold on; 
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.9,0.9,0.9],'EdgeColor','none');


hold on; 
p1 = plot(eesc_return_date(1, :), 'LineWidth', 2);
p2 = plot(eesc_return_date(2, :), 'LineWidth', 2);
p3 = plot(eesc_return_date(3, :), 'LineWidth', 2);
p4 = plot(eesc_return_date(4, :), 'LineWidth', 2);
p5 = plot([1,10], [eesc_return_date(5, 1),eesc_return_date(5, 1)], '--k', 'LineWidth', 2)
xticks(1:10)
xticklabels({'Initial Assessment','CFC-11','CFC-12','Other','CFC-11','CFC-12','Other','CFC-11','CFC-12','Other'})
set(gca, 'FontSize', 12); box on;
legend([p1, p2, p3, p4, p5],'2006 Assessment','2010  Assessment','2014 Assessment','2018 Assessment','2022 Assessment')
ylim([2050,2068])

figure_width = 8; % in inches
figure_height = 8; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = scen_ii;
Figstr = strcat('EESC_Summary_Fig1.pdf'); 
print(gcf, '-dpdf', Figstr);


fighandle = figure; 
plot([eesc_return_date(5, 1),eesc_return_date(5, 1)], [0,5], '--k','LineWidth',2); 
hold on; 
p1 = plot(eesc_return_date(1:4,1)', [4,3,2,1], '.','MarkerSize',20, 'Color', LineCols(1,:)); hold on; 
p2 = plot(eesc_return_date(1:4,4)', [4,3,2,1], '.','MarkerSize',20, 'Color', LineCols(2,:)); hold on; 
p3 = plot(eesc_return_date(1:4,7)', [4,3,2,1], '.','MarkerSize',20, 'Color', LineCols(3,:)); hold on; 
p4 = plot(eesc_return_date(1:4,10)', [4,3,2,1], '.','MarkerSize',20, 'Color', LineCols(4,:));
legend([p1,p2,p3,p4],'Initial OA','LT correction','MF correction','Banks correction')
set(gca, 'FontSize', 12); box on;
yticklabels('')
xlim([2052,2068])
xlabel('EESC return date'); box off
    
figure_width = 8; % in inches
figure_height = 6; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = scen_ii;
Figstr = strcat('EESC_Summary_Fig2.pdf'); 
print(gcf, '-dpdf', Figstr);


fighandle = figure; 
clear corr_by_chem
corr_by_chem(1:4,1) = eesc_return_date(1:4,1);
diff_matrix = eesc_return_date(1:4,2:end) - eesc_return_date(1:4,1:end - 1);
cfc11_diff = sum(diff_matrix(:,1:3:end)')';
cfc12_diff = sum(diff_matrix(:,2:3:end)')';
allother_diff = sum(diff_matrix(:,3:3:end)')';

corr_by_chem(1:4,2) = eesc_return_date(1:4,1) + cfc11_diff ;
corr_by_chem(1:4,3) = eesc_return_date(1:4,1) + cfc11_diff +cfc12_diff;
corr_by_chem(1:4,4) = eesc_return_date(1:4,1) + cfc11_diff +cfc12_diff + allother_diff;

plot([eesc_return_date(5, 1),eesc_return_date(5, 1)], [0,5], '--k','LineWidth',2); 
hold on; 
p1 = plot(corr_by_chem(:,1)', [4,3,2,1], '.','MarkerSize',20, 'Color', LineCols(1,:)); hold on; 
p2 = plot(corr_by_chem(:,2)', [4,3,2,1], '.','MarkerSize',20, 'Color', LineCols(2,:)); hold on; 
p3 = plot(corr_by_chem(:,3)', [4,3,2,1], '.','MarkerSize',20, 'Color', LineCols(3,:)); hold on; 
p4 = plot(corr_by_chem(:,4), [4,3,2,1], '.','MarkerSize',20, 'Color', LineCols(4,:)); hold on;
%p5 = plot(eesc_return_date(1:4,10), [4,3,2,1], '*k')
legend([p1,p2,p3,p4],'Initial OA','cfc-11 update','cfc-12 update','all others')
set(gca, 'FontSize', 12); box on;
yticklabels('')
xlim([2052,2068])
xlabel('EESC return date'); box off
    
figure_width = 8; % in inches
figure_height = 6; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = scen_ii;
Figstr = strcat('EESC_Summary_Fig3.pdf'); 
print(gcf, '-dpdf', Figstr);


%% Explanatory Figure 1: 


%% Explanatory Figure 2: 
 
        %cfc11_banks = Bank_byscen;
        %cfc11_mf = MF_byscen; % years x ozone assessment x scenario
        %cfc11_emiss = Emiss_byscen;

 
 load('/Users/meganlickley/Dropbox (MIT)/Research/Ozone Assessment/Code/July2021/CFC11_fugitive.mat'); 
 
 ShortBank.cfc11 = median(CFC11post.shortBank); 
 MediumBank.cfc11 = median(CFC11post.medBank); 
 LongBank.cfc11 = median(CFC11post.longBank); 
 
 MediumRF.cfc11 = CFC11post.rf_m; 
 LongRF.cfc11 = CFC11post.rf_l;
 
 ShortBank.cfc11(65) = 0;
 MediumBank.cfc11(65) = median(CFC11post.medBank(:,64).*(1 - CFC11post.rf_m'));
 LongBank.cfc11(65) = median(CFC11post.longBank(:,64).*(1 - CFC11post.rf_l'));
 
 TotalBanks.cfc11 = CFC11post.shortBank + CFC11post.medBank + CFC11post.longBank;

 TotalEmiss.cfc11 = CFC11post.emiss_s + CFC11post.emiss_l + CFC11post.emiss_m
 MF.cfc11 = CFC11post.MF;
 
 F_factor = 1.07;
 ppt_millionmols = (1.06*10^32)/(6.02*10^23*10^6*F_factor);
 Molecular_weight = 137.37;
 ppt_to_tonnes =  ppt_millionmols*Molecular_weight;
 LT = 52; 

 for ii = 65:145
    CFC11post.medBank(:,ii) = CFC11post.medBank(:,ii-1).*(1 - CFC11post.rf_m');
    CFC11post.longBank(:,ii) = CFC11post.longBank(:,ii-1).*(1 - CFC11post.rf_l');
    TotalBanks.cfc11(:,ii) = CFC11post.medBank(:,ii) + CFC11post.longBank(:,ii);

    TotalEmiss.cfc11(:,ii) = CFC11post.medBank(:,ii-1).*CFC11post.rf_m' + CFC11post.longBank(:,ii-1).*CFC11post.rf_l';
    MF.cfc11(:,ii) = exp(-1/LT)*MF.cfc11(:,ii-1)+(1/ppt_to_tonnes)*TotalEmiss.cfc11(:,ii-1);
 end
load('/Users/meganlickley/Dropbox (MIT)/Research/Ozone Assessment/Code/July2021/CFC12.mat'); 
 
ShortBank.cfc12a = median(CFC12post.shortBank1); 
ShortBank.cfc12b = median(CFC12post.shortBank2); 
MediumBank.cfc12 = median(CFC12post.medBank); 
LongBank.cfc12 = median(CFC12post.longBank); 

MediumRF.cfc12 = CFC12post.rf_m; 
LongRF.cfc12 = CFC12post.rf_l;

ShortBank.cfc12a(65) = 0; 
ShortBank.cfc12b(65) = 0; 
MediumBank.cfc12(65) = median(CFC12post.medBank(:,64).*(1 - CFC12post.rf_m'));
LongBank.cfc12(65) = median(CFC12post.longBank(:,64).*(1 - CFC12post.rf_l'));

TotalBanks.cfc12 = CFC12post.shortBank1 + CFC12post.shortBank2 + CFC12post.medBank + CFC12post.longBank;
TotalBanks.cfc12(:,65) = CFC12post.medBank(:,64).*(1 - CFC12post.rf_m') + CFC12post.longBank(:,64).*(1 - CFC12post.rf_l');

load('/Users/meganlickley/Dropbox (MIT)/Research/Ozone Assessment/Code/July2021/Halon1301_originalAGAGE_published.mat'); 

HalonRF.H1301 = median(halon1301post.rf);
HalonEmiss.H1301 = median(halon1301post.Emiss);
HalonBank.H1301 = median(halon1301post.Bank);

TotalBanks.H1301 = halon1301post.Bank;

fighandle = figure
subplot(1,3,1)
 MED = 0.001*prctile(TotalBanks.cfc11, 50);
 LB = 0.001*prctile(TotalBanks.cfc11, 5);
 UB = 0.001*prctile(TotalBanks.cfc11, 95);
 boundedline([1955:2019], MED(1:65)',[MED(1:65)'- LB(1:65)',UB(1:65)' - MED(1:65)'],'alpha','cmap',1/255*[148 0 211]);
 hold on; 
 p1 = plot([1955:2019], MED(1:65)', 'Color',1/255*[148 0 211],'LineWidth', 2);
 hold on; 
 p2 = plot([1950:2019], 0.001*squeeze(cfc11_banks(1:70,1,1)),'k','LineWidth',2);
 hold on; 
 plot([2002,2015], 0.001*[squeeze(cfc11_banks(53,1,1)),squeeze(cfc11_banks(66,1,1))], '.','MarkerSize',20, 'Color', 'k'); 
        
 hold on; 
 p3 = plot([1950:2019], 0.001*squeeze(cfc11_banks(1:70,2,1)),'Color', 'r' ,'LineWidth',2);
 %hold on;
 %plot(2008, 0.001*squeeze(cfc11_banks(59,2,1)), '.','MarkerSize',20, 'Color', "#FF8800" )
 
 hold on; 
 p4 = plot([1950:2019], 0.001*squeeze(cfc11_banks(1:70,3,1)),'Color','b','LineWidth',2);
 %hold on;
 %plot(2008, 0.001*squeeze(cfc11_banks(59,3,1)), '.','MarkerSize',20, 'Color', [0 0.6 0.3] )
 
 
 hold on; 
 p5 = plot([1950:2019], 0.001*squeeze(cfc11_banks(1:70,4,1)),'Color',"#FF8800", 'LineWidth',2);
 hold on;
 plot(2008, 0.001*squeeze(cfc11_banks(59,4,1)), '.','MarkerSize',20, 'Color', "#FF8800")

 xlim([1990,2019]); ylim([0,3000]); box on; 
 set(gca, 'FontSize', 12); 
 legend([p2, p3, p4, p5, p1], 'WMO 2006','WMO 2010','WMO 2014', 'WMO 2018','WMO 2022')
ylabel('CFC-11 Banks[Gg]')
%title('CFC-11 Bank Estimates')

subplot(1,3,2)
 MED = 0.001*prctile(TotalBanks.cfc12, 50);
 LB = 0.001*prctile(TotalBanks.cfc12, 5);
 UB = 0.001*prctile(TotalBanks.cfc12, 95);
 boundedline([1955:2019], MED(1:65)',[MED(1:65)'- LB(1:65)',UB(1:65)' - MED(1:65)'],'alpha','cmap',1/255*[148 0 211]);
 hold on; 
 p1 = plot([1955:2019], MED(1:65)', 'Color',1/255*[148 0 211],'LineWidth', 2);
 hold on; 
 p2 = plot([1950:2019], 0.001*squeeze(cfc12_banks(1:70,1,1)),'k','LineWidth',2);
 hold on; 
 plot([2002,2015], 0.001*[squeeze(cfc12_banks(53,1,1)),squeeze(cfc12_banks(66,1,1))], '.','MarkerSize',20, 'Color', 'b'); 
        
 hold on; 
 p3 = plot([1950:2019], 0.001*squeeze(cfc12_banks(1:70,2,1)),'Color', 'r' ,'LineWidth',2);

 hold on; 
 p4 = plot([1950:2019], 0.001*squeeze(cfc12_banks(1:70,3,1)),'Color','b','LineWidth',2);
 
 hold on; 
 p5 = plot([1950:2019], 0.001*squeeze(cfc12_banks(1:70,4,1)),'Color',"#FF8800", 'LineWidth',2);
 hold on;
 plot(2008, 0.001*squeeze(cfc12_banks(59,4,1)), '.','MarkerSize',20, 'Color', "#FF8800")

 xlim([1990,2019])
 set(gca, 'FontSize', 12); ylim([0,3000]); box on; 
 legend([p2, p3, p4, p5, p1], 'WMO 2006','WMO 2010','WMO 2014', 'WMO 2018','WMO 2022')
 ylabel('CFC-12 Banks[Gg]')
%title('CFC-12 Bank Estimates')

 subplot(1,3,3)
 MED = 0.001*prctile(TotalBanks.H1301, 50);
 LB = 0.001*prctile(TotalBanks.H1301, 5);
 UB = 0.001*prctile(TotalBanks.H1301, 95);
 boundedline(halon1301post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',1/255*[148 0 211]);
 hold on; 
 p1 = plot(halon1301post.Years, MED', 'Color',1/255*[148 0 211],'LineWidth', 2);
 hold on; 
 p2 = plot([1950:2019], 0.001*squeeze(halon1301_banks(1:70,1,1)),'k','LineWidth',2);
 hold on; 
 plot([2002,2015], 0.001*[squeeze(halon1301_banks(53,1,1)),squeeze(halon1301_banks(66,1,1))], '.','MarkerSize',20, 'Color', 'b'); 
        
 hold on; 
 p3 = plot([1950:2019], 0.001*squeeze(halon1301_banks(1:70,2,1)),'Color', 'r' ,'LineWidth',2);

 hold on; 
 p4 = plot([1950:2019], 0.001*squeeze(halon1301_banks(1:70,3,1)),'Color','b','LineWidth',2);
 
 hold on; 
 p5 = plot([1950:2019], 0.001*squeeze(halon1301_banks(1:70,4,1)),'Color',"#FF8800", 'LineWidth',2);
 hold on;
 plot(2008, 0.001*squeeze(halon1301_banks(59,4,1)), '.','MarkerSize',20, 'Color', "#FF8800")

 xlim([1990,2019]); ylim([0,150]); box on; 
 set(gca, 'FontSize', 12); 
 legend([p2, p3, p4, p5, p1], 'WMO 2006','WMO 2010','WMO 2014', 'WMO 2018','WMO 2022')
ylabel('Halon-1301 Banks[Gg]');
%title('Halon-1301 Bank Estimates')

figure_width = 18; % in inches
figure_height = 5; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = scen_ii;
Figstr = strcat('Bank_comparison.pdf'); 
print(gcf, '-dpdf', Figstr);




 fighandle = figure
 subplot(3,1,1)
 MED = 0.001*prctile(TotalBanks.cfc11, 50);
 LB = 0.001*prctile(TotalBanks.cfc11, 5);
 UB = 0.001*prctile(TotalBanks.cfc11, 95);
 boundedline([1955:2099], MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.3,0.3,0.3]);
 hold on; 
 p1 = plot([1955:2099], MED', 'Color',[0.3,0.3,0.3],'LineWidth', 2);
 hold on; 
 p2 = plot([1950:2100], 0.001*squeeze(cfc11_banks(:,1,1)),'LineWidth',2);
 hold on; 
 p3 = plot([1950:2100], 0.001*squeeze(cfc11_banks(:,1,2)),'LineWidth',2);
 hold on; 
 p4 = plot([1950:2100], 0.001*squeeze(cfc11_banks(:,1,3)),'LineWidth',2);
 set(gca, 'FontSize', 12); 
 legend([p2, p3, p4, p1], '2006 Assessment','LT update','MF update', 'Bank update')
 xlim([2000,2070]); ylabel('CFC-11 Banks [Gg]')

 subplot(3,1,2)  
 MED = prctile(0.001*TotalEmiss.cfc11, 50);
 LB = prctile(0.001*TotalEmiss.cfc11, 5);
 UB = prctile(0.001*TotalEmiss.cfc11, 95);
 boundedline([1955:2099], MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.3,0.3,0.3]);
 hold on; 
 p1 = plot([1955:2099], MED', 'Color',[0.3,0.3,0.3],'LineWidth', 2);
hold on; 
p2 =  plot([1950:2100], 0.001*squeeze(cfc11_emiss(:,1,1)),'LineWidth',2);
hold on; 
p3 =  plot([1950:2100], 0.001*squeeze(cfc11_emiss(:,1,2)),'LineWidth',2);
hold on; 
p4 =  plot([1950:2100], 0.001*squeeze(cfc11_emiss(:,1,3)),'LineWidth',2);
set(gca, 'FontSize', 12); 
legend([p2, p3, p4, p1], '2006 Assessment','LT update','MF update', 'Bank update')
 xlim([2000,2070]); ylabel('CFC-11 Emissions [Gg yr^{-1}]')

 subplot(3,1,3)  
 MED = prctile(MF.cfc11, 50);
 LB = prctile(MF.cfc11, 5);
 UB = prctile(MF.cfc11, 95);
 boundedline([1955:2099], MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.3,0.3,0.3]);
 hold on; 
 p1 = plot([1955:2099], MED', 'Color',[0.3,0.3,0.3],'LineWidth', 2);
 hold on; 
 p2 =  plot([1950:2100], squeeze(cfc11_mf(:,1,1)),'LineWidth',2);
 hold on; 
 p3 = plot([1950:2100], squeeze(cfc11_mf(:,1,2)),'LineWidth',2);
 hold on; 
 p4 = plot([1950:2100], squeeze(cfc11_mf(:,1,3)),'LineWidth',2);
 set(gca, 'FontSize', 12); 
 legend([p2, p3, p4, p1], '2006 Assessment','LT update','MF update', 'Bank update')
 xlim([2000,2070]); ylabel('Mole Fraction [pmol mol^{-1}')
% This Represents the posterior samples for each of the 

figure_width = 8; % in inches
figure_height = 14; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = scen_ii;
Figstr = strcat('CFC11_explanatory_Fig1.pdf'); 
print(gcf, '-dpdf', Figstr);
