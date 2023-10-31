clear all; close all;
addpath('../functions/');

% Empirical example
% Based on Jarocinski and Karadi `Deconstructing Monetary Policy Surprises - The Role of Information Shocks', American Economic Journal: Macroeconomics
% Written by Mátyás Farkas

  % Reporting
    resp_var = 6;                   % Index of response variable
    innov = 1;                      % Index of innovation
    maxhorz = 36;                   % Maximum horizon to plot
    alpha = 0.1;                    % Significance level
    plot_ylim = [-10 5];             % y-axis limits
    plot_xtick = [1 6:6:maxhorz];   % x-axis ticks
        horzs = 1:maxhorz;

    
%% Run the replication file of JK
 
% JK_rep
load JK_res
irs_jk_orig = squeeze(quantile(irfs_draws(6,1,horzs,:),[0.5 alpha 1-alpha],4));


%%
rhod= linspace(0.8,1,21);
for k = 1:21
    create_HF_dataset(rhod(k));
    
    % Data files for VAR
    file_var = '../data/jk/data.csv';      % VAR series
    
    % Original code from GK example
    % Specification
    
    vars = {'ff4_hf','sp500_hf','gs1','logsp500','us_rgdp','us_gdpdef','ebpnew'};  % Variables in VAR
    p = 12;                                                     % Lag length
    
  
    
    % Bootstrap
    boot_num = 500;                     % # of repetitions
    poolobj = gcp;
    boot_workers = poolobj.NumWorkers;  % # of parallel workers
    rng(20230427, 'twister');           % Set random number seed
    
    
    %% Load data
    
    %dat = innerjoin(readtable(file_var), readtable(file_iv));
    dat = readtable(file_var);
    Y = dat{:,vars}; % Select variables
    Y = Y(~any(isnan(Y),2),:); % Remove missing
    Y = detrend(Y,0); % Remove mean
    
    JK_loop;
    irs_jk_loop = squeeze(quantile(irfs_draws(6,1,horzs,:),[0.5 alpha 1-alpha],4));

    
    %% Lag-augmented local projection
    
    disp('Lag-augmented LP: bootstrapping...');
    [irs_lp_la, ~, cis_dm_lp_la, cis_boot_lp_la] = ...
        ir_estim(Y, p, horzs, ...
        'estimator', 'lp', 'lag_aug', true, ...
        'resp_var', resp_var, 'innov', innov, 'alpha', alpha, ...
        'bootstrap', 'var', 'boot_num', boot_num, ...
        'boot_workers', boot_workers, 'verbose', true);
    
    figure('visible','off');
    subplot(4,1,1);
    plot_band(horzs, irs_lp_la, ...
        cis_dm_lp_la(1,:), cis_dm_lp_la(2,:), ...
        ['Lag-augmented LP, delta method interval with daily AR(1) =' num2str(rhod(k)) ], ...
        plot_ylim, plot_xtick);
    subplot(4,1,2);
    plot_band(horzs, irs_lp_la, ...
        cis_boot_lp_la(1,:,3), cis_boot_lp_la(2,:,3), ...
        ['Lag-augmented LP, percentile-t interval with daily AR(1) =' num2str(rhod(k)) ], ...
        plot_ylim, plot_xtick);
    
          subplot(4,1,3);
    plot_band(horzs, irs_jk_orig(horzs,1)'*100, ...
        irs_jk_orig(horzs,2)'*100, irs_jk_orig(horzs,3)'*100, ...
        ['JK baseline two sign restrictions - Published IRF'], ...
        plot_ylim, plot_xtick);
    
    
     subplot(4,1,4);
    plot_band(horzs, irs_jk_loop(horzs,1)'*100, ...
        irs_jk_loop(horzs,2)'*100, irs_jk_loop(horzs,3)'*100, ...
        ['JK baseline identification with daily AR(1) =' num2str(rhod(k)) ' shocks'], ...
         plot_ylim, plot_xtick);
    
%     drawnow;
      saveas(gcf,['/figures/Rho_', num2str(rhod(k)),'_Lag_augmented_LP.png'])

    
%     %% Non-augmented local projection
%     
%     disp('Non-augmented LP: bootstrapping...');
%     [irs_lp_na, ~, cis_dm_lp_na, cis_boot_lp_na] = ...
%         ir_estim(Y, p, horzs, ...
%         'estimator', 'lp', 'lag_aug', false, ...
%         'resp_var', resp_var, 'innov', innov, 'alpha', alpha, ...
%         'har', @(Y,bw) ewc(Y,bw), ...              % HAR estimator
%         'har_bw', @(T) round(0.4*T.^(2/3)), ...    % HAR bandwidth
%         'har_cv', @(bw) tinv(1-alpha/2,bw), ...    % HAR CV
%         'bootstrap', 'var', 'boot_num', boot_num, ...
%         'boot_workers', boot_workers, 'verbose', true);
%     
%     figure('visible','off');
%     subplot(2,1,1);
%     plot_band(horzs, irs_lp_na, ...
%         cis_dm_lp_na(1,:), cis_dm_lp_na(2,:), ...
%         ['Non-augmented LP, delta method interval with daily AR(1) =' num2str(rhod(k)) ], ...
%         plot_ylim, plot_xtick);
%     subplot(2,1,2);
%     plot_band(horzs, irs_lp_na, ...
%         cis_boot_lp_na(1,:,3), cis_boot_lp_na(2,:,3), ...
%         ['Non-augmented LP, percentile-t interval with daily AR(1) =' num2str(rhod(k)) ], ...
%         plot_ylim, plot_xtick);
% %     drawnow;
%     saveas(gcf,['Rho_', num2str(rhod(k)), '_Non_augmented_LP.png'])
%     

    %% Lag-augmented VAR
    
    disp('Lag-augmented VAR: bootstrapping...');
    [irs_var_na, ~, cis_dm_var_na, cis_boot_var_na] = ...
        ir_estim(Y, p, horzs, ...
        'estimator', 'var', 'lag_aug', true, ...
        'resp_var', resp_var, 'innov', innov, 'alpha', alpha, ...
        'bootstrap', 'var', 'boot_num', boot_num, ...
        'boot_workers', boot_workers, 'verbose', true);
    
    figure('visible','off');
    subplot(4,1,1);
    plot_band(horzs, irs_var_na, ...
        cis_dm_var_na(1,:), cis_dm_var_na(2,:), ...
        ['Lag-augmented VAR, delta method interval with daily AR(1) =' num2str(rhod(k)) ], ...
        plot_ylim, plot_xtick);
    subplot(4,1,2);
    plot_band(horzs, irs_var_na, ...
        cis_boot_var_na(1,:,3), cis_boot_var_na(2,:,3), ...
        ['Lag-augmented VAR, percentile-t interval with daily AR(1) =' num2str(rhod(k)) ], ...
        plot_ylim, plot_xtick);
    
      
          subplot(4,1,3);
    plot_band(horzs, irs_jk_orig(horzs,1)'*100, ...
        irs_jk_orig(horzs,2)'*100, irs_jk_orig(horzs,3)'*100, ...
        ['JK baseline two sign restrictions - Published IRF'], ...
        plot_ylim, plot_xtick);
    
    
     subplot(4,1,4);
    plot_band(horzs, irs_jk_loop(horzs,1)'*100, ...
        irs_jk_loop(horzs,2)'*100, irs_jk_loop(horzs,3)'*100, ...
        ['JK baseline identification with daily AR(1) =' num2str(rhod(k)) ], ...
        plot_ylim, plot_xtick);
    
    
%     drawnow;
    saveas(gcf,['/figures/Rho_', num2str(rhod(k)) ,'_Lag_augmented_VAR.png'])
    
%     %% Non-augmented VAR
%     
%     disp('Non-augmented VAR: bootstrapping...');
%     [irs_var_na, ~, cis_dm_var_na, cis_boot_var_na] = ...
%         ir_estim(Y, p, horzs, ...
%         'estimator', 'var', 'lag_aug', false, ...
%         'resp_var', resp_var, 'innov', innov, 'alpha', alpha, ...
%         'bootstrap', 'var', 'boot_num', boot_num, ...
%         'boot_workers', boot_workers, 'verbose', true);
%     
%     figure('visible','off');
%     subplot(2,1,1);
%     plot_band(horzs, irs_var_na, ...
%         cis_dm_var_na(1,:), cis_dm_var_na(2,:), ...
%         ['Non-augmented VAR, delta method interval with daily AR(1) =' num2str(rhod(k)) ], ...
%         plot_ylim, plot_xtick);
%     subplot(2,1,2);
%     plot_band(horzs, irs_var_na, ...
%         cis_boot_var_na(1,:,3), cis_boot_var_na(2,:,3), ...
%         ['Non-augmented VAR, percentile-t interval with daily AR(1) =' num2str(rhod(k)) ], ...
%         plot_ylim, plot_xtick);
% %     drawnow;
% 
%     
% 
%     saveas(gcf,['Rho_', num2str(rhod(k)) ,'_Non_augmented_VAR.png'])

end
% delete(poolobj);


%% Auxiliary plot function

function plot_band(x, y, lower, upper, plot_title, plot_ylim, plot_xtick)

% Plot IRF and error bands

plot(x, y, '-k', 'LineWidth', 2);
hold on;
plot(x, [lower; upper], 'LineStyle', '--', 'Color', 'k');
line([min(x) max(x)], [0 0], 'Color', 'k');
hold off;

xlim([min(x) max(x)]);
ylim(plot_ylim);
set(gca, 'XTick', plot_xtick);
set(gca, 'FontSize', 12);
title(plot_title, 'FontSize', 12);

end


%% Settings from JK - unchanged
% 
% % SAMPLE AND DATA: spl, modname, addvar
% spl = [1984 2; 2016 12];
% %spl = [1984 2; 2008 12];% Dec2008 ZLB reached
% %spl = [1990 2; 2016 12];% Feb1999 surprises start
% %spl = [1979 7; 2016 12];% GertlerKaradi2015 sample
% 
% modname = 'us1'; %'us1','ea1','us2','ea2'
% addvar = ''; %'exp_gdp_12m','exp_cpi_12m','bkeven05','gs10','sven5f5'
% 
% % IDENTIFICATION: idscheme, mnames
% idscheme = 'sgnm2'; %'chol','sgnm2','sgnm2strg','supdem'
% 
% mnames = {'ff4_hf','sp500_hf'}; % US baseline
% %mnames = {'ff4_hf'};
% %mnames = {'pmnegm_ff4sp500','pmposm_ff4sp500'}; % poor man's sign restrictions
% %mnames = {'ff4_hf','sp500_hf','dbkeven02_d'}; % for supdem identification
% %mnames = {'pc1ff1_hf','usstocks1_hf'}; % VAR with factors (Online Appendix C.4)
% %mnames = {'pmnegm_pc1ff1usstocks1','pmposm_pc1ff1usstocks1'}; % VAR with factors, poor man's shocks (Online Appendix C.4)
% 
% %mnames = {'eureon3m_hf','stoxx50_hf'}; % euro area baseline
% %mnames = {'eureon3m_hf'};
% %mnames = {'pmnegm_eureon3mstoxx50','pmposm_eureon3mstoxx50'}; % poor man's sign restrictions
% %mnames = {'eureon3m_hf','stoxx50_hf','deurinflswap2y_d'}; % for supdem identification
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % PRIOR
% prior.lags = 12;
% prior.minnesota.tightness = .2;
% prior.minnesota.decay = 1;
% prior.Nm = length(mnames);
% 
% % create the output folder based on the calling filename
% st = dbstack; pathout = st(end).name; clear st
% mkdir(pathout)
% pathout = [pathout '/'];
% 
% % functions for time operations
% ym2t = @(x) x(1)+x(2)/12 - 1/24; % convert [year month] into time
% t2datestr = @(t) [num2str(floor(t)) char(repmat(109,length(t),1)) num2str(12*(t-floor(t)+1/24),'%02.0f')];
% t2ym = @(t) [floor(t) round(12*(t-floor(t)+1/24))]; % convert time into [year month]
% ymdif = @(x1,x2) (x2(1)-x1(1))*12+x2(2)-x1(2);
% findym = @(x,t) find(abs(t-ym2t(x))<1e-6); % find [year month] in time vector t
% 
% % Gibbs sampler settings
% gssettings.ndraws = 4000;
% gssettings.burnin = 4000;
% gssettings.saveevery = 4;
% gssettings.computemarglik = 0;
% 
% % detect poorman shocks and if yes override the identification to choleski
% if (~isempty(strfind(mnames{1},'neg')) && ~isempty(strfind(mnames{2},'pos')))
%     idscheme = 'chol';
% end
% 
% % nice names
% mdict = {'eureon3m_', '  surprise in\newline3m Eonia swaps';
%     'stoxx50_', 'surprise in\newline Euro Stoxx 50';
%     'ff4_', '  surprise in\newline3m ff futures';
%     'sp500_', 'surprise in\newlineS&P500';
%     'pc1ff1_', 'surprise in\newlinepolicy ind.';
%     'usstocks1_', 'surprise in\newline1pc of stocks'};
% mnames_nice = applydictregexp(mnames, mdict);
% mylimits = nan(length(mnames),2);
% 
% % define y
% ny1 = 5;
% switch modname
%     case 'us1'
%         ynames = {'gs1','logsp500','us_rgdp','us_gdpdef','ebpnew'};
%     case 'us2'
%         ynames = {'gs1','logsp500','us_ip','us_cpi','ebpnew'};
%     case 'ea1'
%         ynames = {'de1y_haver','stoxx50','ea_rgdp','ea_gdpdef','ea_bbb_oas_all_fred'};
%     case 'ea2'
%         ynames = {'de1y_haver','stoxx50','ea_ip_excl_constr','hicp','ea_bbb_oas_all_fred'};
%     otherwise
%         disp(modname), error('unknown modname');
% end
% 
% if ~isempty(addvar)
%     try
%         findstrings(addvar,ynames);
%     catch
%         ynames = [ynames {addvar}]; modname = [modname '_' addvar];
%     end
% end
% 
% % nice_names, yylimits, nonst
% dictfname = '../data/jk/ydict.csv';
% fileID = fopen(dictfname);
% ydict = textscan(fileID,'%s %q %d %f %f','Delimiter',',','HeaderLines',1);
% fclose(fileID);
% ynames_nice = applydict(ynames, [ydict{1} ydict{2}]);%ynames_nice = ynames;
% yylimits = [ydict{4} ydict{5}];
% yylimits = yylimits(findstrings(ynames,ydict{1}),:);
% nonst = ydict{3};
% nonst = nonst(findstrings(ynames,ydict{1}),:);
% 
% % load data
% datafname = '../data/jk/data.csv';
% data.Nm = length(mnames);
% data.names = [mnames ynames];
% d = importdata(datafname); dat = d.data; txt = d.colheaders;
% tbeg = find(dat(:,1)==spl(1,1) & dat(:,2)==spl(1,2)); if isempty(tbeg), tbeg=1; end
% tend = find(dat(:,1)==spl(2,1) & dat(:,2)==spl(2,2)); if isempty(tend), tend=size(dat,1); end
% ysel = findstrings(data.names, txt(1,:));
% data.y = dat(tbeg:tend, ysel);
% data.w = ones(size(data.y,1),1);
% data.time = linspace(ym2t(dat(tbeg,1:2)), ym2t(dat(tend,1:2)), size(data.y,1))';
% clear d dat txt tbeg tend ysel
% 
% % check ydata for missing values, determine the sample
% datatemp = checkdata(data, t2datestr, 1:data.Nm);
% idspl = [t2datestr(datatemp.time(1)) '-' t2datestr(datatemp.time(end))];
% clear datatemp
% 
% % output file names
% fname = [pathout modname '_' strjoin(mnames,'_') '_' idspl '_' idscheme];
% diary([fname '.txt'])
% data = checkdata(data, t2datestr, 1:data.Nm); % again, for the diary
% plot_y;
% 
% % print the correlation matrix of m
% table = corr(data.y(:,1:data.Nm),'rows','pairwise');
% in.cnames = strvcat(mnames);
% in.rnames = strvcat(['correl.:' mnames]);
% in.width = 200;
% mprint(table,in)
% 
% % complete the minnesota prior
% prior.minnesota.mvector = [zeros(data.Nm,1); nonst];
% 
% % replace NaNs with zeros in the initial condition
% temp = data.y(1:prior.lags,:); temp(isnan(temp)) = 0; data.y(1:prior.lags,:) = temp;
% 
% % drop the shocks before February 1994
% %id = data.time<ym2t([1994 2])-1e-6; data.y(id,1:data.Nm) = NaN;
% 
% % estimate the VAR
% %data.y(isnan(data.y)) = 0; res = VAR_dummyobsprior(data.y,data.w,gssettings.ndraws,prior);
% res = VAR_withiid1kf(data, prior, gssettings);
% 
% savedata([fname '_data.csv'], data, t2ym)
% %% identification
% MAlags = 36;
% N = length(data.names);
% 
% switch idscheme
%     case 'chol'
%         shocknames = data.names;
%         irfs_draws = NaN(N,N,MAlags,gssettings.ndraws);
%         for i = 1:gssettings.ndraws
%             betadraw = res.beta_draws(1:end-size(data.w,2),:,i);
%             sigmadraw = res.sigma_draws(:,:,i);
%             response = impulsdtrf(reshape(betadraw',N,N,prior.lags), chol(sigmadraw), MAlags);
%             irfs_draws(:,:,:,i) = response;
%         end
%         ss = 1;
%         if length(mnames)>1 && ((~isempty(strfind(mnames{1},'neg')) && ~isempty(strfind(mnames{2},'pos'))) || ~isempty(strfind(mnames{2},'_signrestr'))), ss = 1:2; end
%     case 'sgnm2' % baseline two sign restrictions
%         shocknames = [{'mon.pol.', 'CBinfo'} mnames(2+1:end) ynames];
%         dims = {[1 2]};
%         imonpol = 1; inews = 2;
%         test_restr = @(irfs)... %% restrictions by shock (i.e. by column):
%             irfs(1,imonpol,1) > 0 && irfs(2,imonpol,1) < 0 &&... % mp
%             irfs(1,inews,1) > 0 && irfs(2,inews,1) > 0; % cbi
%         b_normalize = ones(1,N);
%         max_try = 1000;
%         disp(test_restr)
%         irfs_draws = resirfssign(res, MAlags, dims, test_restr, b_normalize, max_try);
%         %[irfs_draws, irfs_l_draws, irfs_u_draws] = resirfssign_robust(res, MAlags, dims, test_restr, b_normalize, max_try);
%         ss = 1:2;
%     case 'sgnm2strg' % strong instrument restriction
%         % imposes that the THIRD variable goes up after mp shock
%         shocknames = [{'mon.pol.', 'CBinfo'} ynames];
%         dims = {[1 2]};
%         iyld = 3;
%         imonpol = 1; inews = 2;
%         test_restr = @(irfs)... %% restrictions by shock (i.e. by column):
%             irfs(1,imonpol,1) > 0 && irfs(2,imonpol,1) < 0 && irfs(iyld,imonpol,1)>0.01 &&... % mp
%             irfs(1,inews,1) > 0 && irfs(2,inews,1) > 0; % cbi
%         b_normalize = ones(1,N);
%         max_try = 1000;
%         disp(test_restr)
%         irfs_draws = resirfssign(res, MAlags, dims, test_restr, b_normalize, max_try);
%         ss = 1:2;
%     case 'supdem' % disentangle CB info about supply and demand
%         % requires: 1. interest rate; 2. stock price; 3. break-even inflation
%         % CBinfosup shock moves stock price down but break-even inflation up
%         shocknames = [{'mon.pol.', 'CBinfodem', 'CBinfosup'} mnames(3+1:end) ynames];
%         dims = {[1 2 3]};
%         imonpol = 1; inews = 2; isup = 3;
%         test_restr = @(irfs)... %% restrictions by shock (i.e. by column):
%             irfs(1,imonpol,1) > 0 && irfs(2,imonpol,1) < 0 && irfs(3,imonpol,1) < 0 &&... % mp
%             irfs(1,inews,1) > 0 && irfs(2,inews,1) > 0 && irfs(3,inews,1) > 0 &&... % info demand
%             irfs(1,isup,1) > -100 && irfs(2,isup,1) > 0 && irfs(3,isup,1) < 0; % info supply
%         b_normalize = ones(N,1); b_normalize(3) = -1;
%         max_try = 5000;
%         disp(test_restr)
%         irfs_draws = resirfssign(res, MAlags, dims, test_restr, b_normalize, max_try);
%         ss = 1:3;
% end
% 
% %% reporting
% 
% % Reporting
% resp_var = 4;                   % Index of response variable
% innov = 5;                      % Index of innovation
% maxhorz = 48;                   % Maximum horizon to plot
% alpha = 0.1;                    % Significance level
% plot_ylim = [-4 4];             % y-axis limits
% plot_xtick = [1 6:6:maxhorz];   % x-axis ticks
% 
% % Bootstrap
% boot_num = 500;                     % # of repetitions
% poolobj = gcp;
% boot_workers = poolobj.NumWorkers;  % # of parallel workers
% rng(20200711, 'twister');           % Set random number seed
% 
% 
% %% Load data
% 
% dat = innerjoin(readtable(file_var), readtable(file_iv));
% Y = dat{:,vars}; % Select variables
% Y = Y(~any(isnan(Y),2),:); % Remove missing
% Y = detrend(Y,0); % Remove mean
% horzs = 1:maxhorz;
% 
% 
% %% Lag-augmented local projection
% 
% disp('Lag-augmented LP: bootstrapping...');
% [irs_lp_la, ~, cis_dm_lp_la, cis_boot_lp_la] = ...
%     ir_estim(Y, p, horzs, ...
%     'estimator', 'lp', 'lag_aug', true, ...
%     'resp_var', resp_var, 'innov', innov, 'alpha', alpha, ...
%     'bootstrap', 'var', 'boot_num', boot_num, ...
%     'boot_workers', boot_workers, 'verbose', true);
% 
% figure;
% subplot(2,1,1);
% plot_band(horzs, irs_lp_la, ...
%     cis_dm_lp_la(1,:), cis_dm_lp_la(2,:), ...
%     'Lag-augmented LP, delta method interval', ...
%     plot_ylim, plot_xtick);
% subplot(2,1,2);
% plot_band(horzs, irs_lp_la, ...
%     cis_boot_lp_la(1,:,3), cis_boot_lp_la(2,:,3), ...
%     'Lag-augmented LP, percentile-t interval', ...
%     plot_ylim, plot_xtick);
% drawnow;
% 
% 
% %% Non-augmented local projection
% 
% disp('Non-augmented LP: bootstrapping...');
% [irs_lp_na, ~, cis_dm_lp_na, cis_boot_lp_na] = ...
%     ir_estim(Y, p, horzs, ...
%     'estimator', 'lp', 'lag_aug', false, ...
%     'resp_var', resp_var, 'innov', innov, 'alpha', alpha, ...
%     'har', @(Y,bw) ewc(Y,bw), ...              % HAR estimator
%     'har_bw', @(T) round(0.4*T.^(2/3)), ...    % HAR bandwidth
%     'har_cv', @(bw) tinv(1-alpha/2,bw), ...    % HAR CV
%     'bootstrap', 'var', 'boot_num', boot_num, ...
%     'boot_workers', boot_workers, 'verbose', true);
% 
% figure;
% subplot(2,1,1);
% plot_band(horzs, irs_lp_na, ...
%     cis_dm_lp_na(1,:), cis_dm_lp_na(2,:), ...
%     'Non-augmented LP, delta method interval', ...
%     plot_ylim, plot_xtick);
% subplot(2,1,2);
% plot_band(horzs, irs_lp_na, ...
%     cis_boot_lp_na(1,:,3), cis_boot_lp_na(2,:,3), ...
%     'Non-augmented LP, percentile-t interval', ...
%     plot_ylim, plot_xtick);
% drawnow;
% 
% 
% %% Non-augmented VAR
% 
% disp('Non-augmented VAR: bootstrapping...');
% [irs_var_na, ~, cis_dm_var_na, cis_boot_var_na] = ...
%     ir_estim(Y, p, horzs, ...
%     'estimator', 'var', 'lag_aug', false, ...
%     'resp_var', resp_var, 'innov', innov, 'alpha', alpha, ...
%     'bootstrap', 'var', 'boot_num', boot_num, ...
%     'boot_workers', boot_workers, 'verbose', true);
% 
% figure;
% subplot(2,1,1);
% plot_band(horzs, irs_var_na, ...
%     cis_dm_var_na(1,:), cis_dm_var_na(2,:), ...
%     'Non-augmented VAR, delta method interval', ...
%     plot_ylim, plot_xtick);
% subplot(2,1,2);
% plot_band(horzs, irs_var_na, ...
%     cis_boot_var_na(1,:,3), cis_boot_var_na(2,:,3), ...
%     'Non-augmented VAR, percentile-t interval', ...
%     plot_ylim, plot_xtick);
% drawnow;
% 
% delete(poolobj);
% 
% 
% %% Auxiliary plot function
% 
% function plot_band(x, y, lower, upper, plot_title, plot_ylim, plot_xtick)
% 
% % Plot IRF and error bands
% 
% plot(x, y, '-k', 'LineWidth', 2);
% hold on;
% plot(x, [lower; upper], 'LineStyle', '--', 'Color', 'k');
% line([min(x) max(x)], [0 0], 'Color', 'k');
% hold off;
% 
% xlim([min(x) max(x)]);
% ylim(plot_ylim);
% set(gca, 'XTick', plot_xtick);
% set(gca, 'FontSize', 12);
% title(plot_title, 'FontSize', 14);
% 
% end