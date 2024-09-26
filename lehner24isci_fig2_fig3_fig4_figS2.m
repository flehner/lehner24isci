close all
clear all

% ---------------------------
% lehner24isci_fig2_fig3_fig4_figS2.m
% ---------------------------
% Code to produce Figs. 2-4 and Fig. S2 in Lehner (2024)
%
% Citation:
% Lehner, F. (2024):
% Climate model large ensembles as test beds for applied compound event research
% iScience, DOI: TBD
%
% Notes on code:
% - Only does calculations and plotting, not pre-processing
% - Requires pre-processed data, which is provided in the 'data' subfolder:
%   - Upper Colorado River Basin precipitation and temperature (PRISM)
%   - Upper Colorado River naturalized flow (Bureau of Reclamation)
%   - Global average surface temperature (BEST)
%   - Combined Lake Mead and Lake Powell reservoir storage (Bureau of Reclamation), transferred from elevation to volume
%   - Temperature from SMILEs (CMIP5+6), averaged globally
%   - Temperature and precipitation from SMILEs (CMIP5+6), averaged over the Upper Colorado River Basin
% ---------------------------

% -- select which figures to plot:
fig2d = 1; % 0=no, 1=yes
fig2c = 0; % 0=no, 1=yes
fig3  = 0; % 0=no, 1=yes
fig4  = 0; % 0=no, 1=yes
figS2 = 0; % 0=no, 1=yes

if fig4 == 1
  recalc1 = 1; % recalculate all risk surfaces
else
  recalc1 = 0; % don't recalculate all risk surfaces
end

pathin = '~/Dropbox/publication/lehner22_smile_perspective/lehner24isci/data/';

refstart  = 1950;
refende   = 2099;
time      = refstart+1:refende; % +1 because of water year notation


% == OBSERVATIONS =================================================================
% -- Upper Colorado precip
obs_data_name = 'PRISM';
pr_obs 	= squeeze(ncread([pathin '/prism_precip_189501-202309_UCRB_ts.nc'],['precip']));
time_Pobs	= 1896:2023;
for i = 2:ceil(length(pr_obs)/12)
  pr_obs_wy(i-1) = sum(pr_obs((i-2)*12+10:(i-1)*12+9));
end
Pobs            = pr_obs_wy;
Pobs_anom       = Pobs-nanmean(Pobs(1950-time_Pobs(1)+1:end));
Pobs_anom_rel   = (Pobs/nanmean(Pobs(1950-time_Pobs(1)+1:end))) * 100; % precip as percent of normal

% -- Upper Colorado temperature
% -- PRISM
tas_obs   = squeeze(ncread([pathin '/prism_tmean_189501-202309_UCRB_ts.nc'],['tmean']));
time_Tobs	= 1896:2023;
% --
for i = 2:ceil(length(tas_obs)/12)
  tas_obs_wy(i-1) = nanmean(tas_obs((i-2)*12+10:(i-1)*12+9));
end
Tobs_anom      = tas_obs_wy-nanmean(tas_obs_wy(refstart-time_Tobs(1)+1:end)); % anomaly

% -- Upper Colorado (Lees Ferry) naturalized flow (https://www.usbr.gov/lc/region/g4000/NaturalFlow/current.html):
data              = load([pathin '/NaturalFlows1906-2019_20210420_edit.csv']);
timem_Qobs        = 1905+10/12:1/12:2022+9/12;
time_Qobs         = 1906:2019;
ro_obs            = data*1e-6; % convert to million acre feet (MAF)
% -- make water year totals:
for i = 1:length(time_Qobs)
  ro_obs_wy(i)  = nansum(ro_obs((i-1)*12+1:i*12));
end
Qobs      = ro_obs_wy;
% -- provisional additional years from https://www.usbr.gov/lc/region/g4000/NaturalFlow/provisional.html
time_Qobs = 1906:2023;
Qobs      = [Qobs 9.558765 7.818000 10.110000 18.603000];

% -- Global temperature
% -- BEST (https://berkeley-earth-temperature.s3.us-west-1.amazonaws.com/Global/Land_and_Ocean_complete.txt)
data = readtable([pathin '/Land_and_Ocean_complete_edit.txt']);
tas_obs_global    = data.Var3;
time_Tobs_global  = unique(data.Var1);
time_Tobs_global  = time_Tobs_global(1:end-1); % 2024 not complete/needed
for i = 2:floor(length(tas_obs_global)/12)
  Tobs_anom_global(i-1) = nanmean(tas_obs_global((i-2)*12+10:(i-1)*12+9)); % water year mean
end
Tobs_anom_global  = Tobs_anom_global-nanmean(Tobs_anom_global(refstart-time_Tobs_global(1)+1:end));
Tobs_anom_global  = Tobs_anom_global(refstart-time_Tobs_global(1):end);


% -- find common time period
[time_common,idxQ,idxP] = intersect(time_Qobs,time_Pobs);
Pobs = Pobs(idxP);
Pobs_anom_rel = Pobs_anom_rel(idxP);
Pobs_anom = Pobs_anom(idxP);
Qobs = Qobs(idxQ);
[tmp98,tmp99,idxT] = intersect(time_common,time_Tobs);
Tobs_anom = Tobs_anom(idxT);
% -- limit to 1950-present to be consistent with models
ind = find(time_common==refstart);
Tobs_anom = Tobs_anom(ind:end);
Pobs_anom_rel = Pobs_anom_rel(ind:end);
Pobs_anom = Pobs_anom(ind:end);
Pobs = Pobs(ind:end);
Qobs = Qobs(ind:end);
time_common = refstart:time_common(end);

% -- detrend Tobs and Pobs with global mean temperature:
b                 = regress(Tobs_anom',[ones(size(Tobs_anom')) Tobs_anom_global']);
Tobs_anom_y       = b(1) + b(2).*Tobs_anom_global;
Tobs_anom_dtr     = Tobs_anom-Tobs_anom_y;
b                 = regress(Pobs_anom_rel',[ones(size(Pobs_anom_rel')) Tobs_anom_global']);
Pobs_anom_rel_y   = b(1) + b(2).*Tobs_anom_global;
Pobs_anom_rel_dtr = Pobs_anom_rel-Pobs_anom_rel_y;
b                 = regress(Pobs_anom',[ones(size(Pobs_anom')) Tobs_anom_global']);
Pobs_anom_y       = b(1) + b(2).*Tobs_anom_global;
Pobs_anom_dtr     = Pobs_anom-Pobs_anom_y;


% -- actual Lake Mead+Powell storage:
mead_powell_combo = load([pathin '/mead_and_powell_combined_updated.csv']);

% -- run reservoir model
[Sobs,Sfracobs,Robs,Eobs,Cobs,Dobs,Defobs] = resmod(Qobs,Pobs,Tobs_anom,Qobs); % 'resmod' takes absolute Q, absolute P, and anomalous T



% == SMILEs data =================================================================
model_scen_pairs = ...
{'CESM1-CAM5','rcp85','rcp85';...
'CanESM2','rcp85','rcp85';...
'CSIRO-Mk3-6-0','rcp85','rcp85';...
'EC-EARTH','rcp85','rcp85';...
'GFDL-CM3','rcp85','rcp85';...
'GFDL-ESM2M','rcp85','rcp85';...
'MPI-ESM','rcp85','rcp85';...
'CESM2','ssp370','ssp370';...
'ACCESS-ESM1-5','ssp245','ssp370';...
'CanESM5','ssp370','ssp370';...
'EC-Earth3','ssp245','ssp370';...
'IPSL-CM6A-LR','ssp245','ssp370';...
'MIROC6','ssp126','ssp370';...
'MIROC-ES2L','ssp245','ssp370';...
'GFDL-SPEAR-MED','ssp585','ssp370';...
'E3SMv2','ssp370','ssp370';...
'MPI-ESM1-2-LR','ssp126','ssp370'};
models        = cellstr(model_scen_pairs(:,1));
scen_smiles   = cellstr(model_scen_pairs(:,2));
scen2_smiles  = cellstr(model_scen_pairs(:,3));
models_select = [1:3 5:11 13:14 16:17]; %1:length(models)
models      = models(models_select);
scen_smiles = scen_smiles(models_select);
ensmem  = [40,50,30,16,20,30,100,100,40,25,21,11,50,30,30,20,30];
ensmem  = zeros(1,length(models)) + min(ensmem(models_select)); % limit all ensemble sizes to e.g. 16, 20, 100
ensmem2 = [40,50,30,20,30,100,100,40,25,8,50,10,20,30]; % actual ensemble size of rcp85 and ssp370 simulations (only used for model eval)
ensmem2(ensmem2>min(ensmem)) = min(ensmem); % cap ensemble sizes of rcp85 and ssp370 to same as other SMILEs
start   = [1920,1950,1850,1860,1920,1950,1850,1850,1850,1850,1850,1850,1850,1850,1920,1850,1850];
ende    = [2100,2100,2100,2100,2100,2100,2099,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100];
start   = start(models_select);
ende    = ende(models_select);
Tanom   = NaN(length(models),max(ensmem),refende-refstart);
Traw    = Tanom;
Praw    = Tanom;
Praw_anom = Tanom;
for m = 1:length(models)
  % ---
  vari  = 'tas';
  pathin = '/Users/fl439/Dropbox/work/reservoir_model/';
  if strcmp(scen_smiles{m},'rcp85')==1 || strcmp(models{m},'CESM2')==1 || strcmp(models{m},'E3SMv2')==1
    in    = dir([pathin vari '_Amon_' models{m} '_*ens*_g025_UCRB_ts_anom.nc']);
    in2   = in;
  else
    in    = dir([pathin '/cmip6-ng/' vari '_Amon_' models{m} '_' scen_smiles{m} '_*ens*_g025_UCRB_ts_anom.nc']); % scenario selection to recover CMIP6 spread
    in2   = dir([pathin '/cmip6-ng/' vari '_Amon_' models{m} '_' scen2_smiles{m} '_*ens*_g025_UCRB_ts_anom.nc']); % consistent scenario (rcp85 and ssp370)
  end
  % -- scenario selection to recover CMIP6 spread
  tmp0	= squeeze(ncread([in.folder '/' in.name],vari));
  for e = 1:ensmem(m)
    tmp = wy_mean(tmp0(e,:));
    tmp = tmp(refstart-start(m)+1:refende-start(m)); % cut to common length
    Traw(m,e,:) = tmp;
    Tanom(m,e,:)= tmp-nanmean(tmp(1:time_Tobs(end)-refstart)); % reference to 1950-1979
  end
  % -- scenario selection to ensure consistency (rcp85 and ssp370)
  tmp0	= squeeze(ncread([in2.folder '/' in2.name],vari));
  ensmem_tmp = length(tmp0(:,1));
  % ['m=' num2str(m) ' / ' num2str(ensmem_tmp)]
  for e = 1:ensmem2(m)
    % -- fill up with NaN if ensemble is less than 20
    if e > ensmem_tmp
      tmp = wy_mean(tmp0(ensmem_tmp,:))*NaN;
    else
      tmp = wy_mean(tmp0(e,:));
    end
    tmp = tmp(refstart-start(m)+1:refende-start(m)); % cut to common length
    Traw2(m,e,:) = tmp;
    Tanom2(m,e,:)= tmp-nanmean(tmp(1:time_Tobs(end)-refstart)); % reference to 1950-1979
  end
  % ---
  vari  = 'pr';
  tf    = 86400*30.4; % --> mm/month (like obs)
  if strcmp(models{m},'CESM2')==1 || strcmp(models{m},'E3SMv2')==1
    tf  = 86400*30.4 * 1e3; % --> mm/month (like obs)
  end
  if strcmp(scen_smiles{m},'rcp85')==1 || strcmp(models{m},'CESM2')==1 || strcmp(models{m},'E3SMv2')==1
    in    = dir([pathin vari '_Amon_' models{m} '_*ens*_g025_UCRB_ts_anom.nc']);
    in2   = in;
  else
    in    = dir([pathin '/cmip6-ng/' vari '_Amon_' models{m} '_' scen_smiles{m} '_*ens*_g025_UCRB_ts_anom.nc']); % scenario selection to recover CMIP6 spread
    in2   = dir([pathin '/cmip6-ng/' vari '_Amon_' models{m} '_' scen2_smiles{m} '_*ens*_g025_UCRB_ts_anom.nc']); % consistent scenario (rcp85 and ssp370)
  end
  % -- scenario selection to recover CMIP6 spread
  tmp0	= squeeze(ncread([in.folder '/' in.name],vari)) * tf;
  for e = 1:ensmem(m)
    tmp                   = wy(tmp0(e,:));
    tmp                   = tmp(refstart-start(m)+1:refende-start(m)); % cut to common length
    Praw(m,e,:)           = tmp;
    Praw_anom(m,e,:)      = tmp-nanmean(tmp(1:time_Pobs(end)-refstart)); % absolute anomalies to 1950-2023 mean
    Padj(m,e,:)           = Praw_anom(m,e,:) * (std(Pobs)/std(tmp(1:time_Pobs(end)-refstart))) + mean(Pobs); % adjust variance and mean to match observations
    Padj_anom(m,e,:)      = Padj(m,e,:)-nanmean(Padj(m,e,1:time_Pobs(end)-refstart)); % absolute anomalies to 1950-2023 mean
    Padj_anom_rel(m,e,:)  = (Padj(m,e,:)/nanmean(Padj(m,e,1:time_Pobs(end)-refstart)))*100; % relative anomalies to 1950-2023 mean
  end
  % -- scenario selection to ensure consistency (rcp85 and ssp370)
  tmp0	= squeeze(ncread([in2.folder '/' in2.name],vari)) * tf;
  ensmem_tmp = length(tmp0(:,1));
  for e = 1:ensmem2(m)
    % -- fill up with NaN if ensemble is less than 20
    if e > ensmem_tmp
      tmp = wy(tmp0(ensmem_tmp,:))*NaN;
    else
      tmp = wy(tmp0(e,:));
    end
    tmp                   = tmp(refstart-start(m)+1:refende-start(m)); % cut to common length
    Praw2(m,e,:)          = tmp;
    Praw2_anom(m,e,:)     = tmp-nanmean(tmp(1:time_Pobs(end)-refstart)); % absolute anomalies to 1950-2023 mean
    Padj2(m,e,:)          = Praw2_anom(m,e,:) * (std(Pobs)/std(tmp(1:time_Pobs(end)-refstart))) + mean(Pobs); % adjust variance and mean to match observations
    Padj2_anom(m,e,:)     = Padj2(m,e,:)-nanmean(Padj2(m,e,1:time_Pobs(end)-refstart)); % absolute anomalies to 1950-2023 mean
    Padj2_anom_rel(m,e,:) = (Padj2(m,e,:)/nanmean(Padj2(m,e,1:time_Pobs(end)-refstart)))*100; % relative anomalies to 1950-2023 mean
  end
  % ---
  % -- global temperature
  vari  = 'tas';
  if strcmp(scen_smiles{m},'rcp85')==1 || strcmp(models{m},'CESM2')==1 || strcmp(models{m},'E3SMv2')==1
    in    = dir([pathin vari '_Amon_' models{m} '_*ens*_g025_ts_anom.nc']);
    in2   = in;
  else
    in    = dir([pathin '/cmip6-ng/' vari '_Amon_' models{m} '_' scen_smiles{m} '_*ens*_g025_ts_anom.nc']); % scenario selection to recover CMIP6 spread
    in2   = dir([pathin '/cmip6-ng/' vari '_Amon_' models{m} '_' scen2_smiles{m} '_*ens*_g025_ts_anom.nc']); % consistent scenario (rcp85 and ssp370)
  end
  % -- scenario selection to recover CMIP6 spread
  tmp0	= squeeze(ncread([in.folder '/' in.name],vari));
  for e = 1:ensmem(m)
    tmp = wy_mean(tmp0(e,:));
    tmp = tmp(refstart-start(m)+1:refende-start(m)); % cut to common length
    Traw_global(m,e,:) = tmp;
    Tanom_global(m,e,:)= tmp-nanmean(tmp(1:time_Tobs(end)-refstart)); % reference to 1950-1979
  end
  % -- scenario selection to ensure consistency (rcp85 and ssp370)
  tmp0	= squeeze(ncread([in2.folder '/' in2.name],vari));
  ensmem_tmp = length(tmp0(:,1));
  for e = 1:ensmem2(m)
    % -- fill up with NaN if ensemble is less than 20
    if e > ensmem_tmp
      tmp = wy_mean(tmp0(ensmem_tmp,:))*NaN;
    else
      tmp = wy_mean(tmp0(e,:));
    end
    tmp = tmp(refstart-start(m)+1:refende-start(m)); % cut to common length
    Traw2_global(m,e,:) = tmp;
    Tanom2_global(m,e,:)= tmp-nanmean(tmp(1:time_Tobs(end)-refstart)); % reference to 1950-1979
  end
  % ---
  vari = 'mrro';
  tmp0 = squeeze(ncread([pathin vari '_Lmon_CESM1-CAM5_ens_g025_UCRB_ts_anom.nc'],vari)) * tf;
  for e = 1:40
    tmp = wy(tmp0(e,:));
    tmp = tmp*(std(Qobs)/std(tmp(1:time_Tobs(end)-refstart)))+mean(Qobs); % adjust variance and mean to match observations
    Q0(e,:) = tmp(refstart-1920+1:refende-1920);
  end
end

return

% -- DETRENDING --
% -- detrend P and T with smoothed ensemble mean:
for m = 1:length(models)
  % -- unadjusted model P, detrended and mean of observations added back in
  em = squeeze(nanmean(Praw(m,:,:)));
  Praw_dtr(m,:,:) = squeeze(Praw(m,:,:)) - polyval(polyfit(1:length(em),em,4),1:length(em)); % + mean(Pobs);
  % -- adjusted model P (variance and mean from observations), detrended
  em = squeeze(nanmean(Padj_anom(m,:,:)));
  Padj_anom_dtr(m,:,:) = squeeze(Padj_anom(m,:,:)) - polyval(polyfit(1:length(em),em,4),1:length(em));
  % -- relative adjusted model P (variance and mean from observations), detrended
  em = squeeze(nanmean(Padj_anom_rel(m,:,:)));
  Padj_anom_rel_dtr(m,:,:) = squeeze(Padj_anom_rel(m,:,:)) - polyval(polyfit(1:length(em),em,4),1:length(em)) + 100;
  % -- T is already in anomaly space, so no need to add anything back in after detrending
  em = squeeze(nanmean(Tanom(m,:,:)));
  Tanom_dtr(m,:,:) = squeeze(Tanom(m,:,:)) - polyval(polyfit(1:length(em),em,4),1:length(em));
end

% -- regression model to calculate Q from P and T and previous-year Q
% -- determine CESM1-CAM5 runoff coefficients to be used for other models
m   = 1;
e   = 2;
y2  = sqrt(Q0(e,1:end)');
y   = y2(2:end);
x1  = squeeze(Praw(m,e,2:end));
x2  = squeeze(Tanom(m,e,2:end));
Ptmp = (x1/nanmean(x1(1:70)))*100;
Ttmp = x2-nanmean(x2(1:70));
% b      = regress(y,[ones(size(y)) x1 x2]);
b      = regress(y,[ones(size(y)) Ptmp Ttmp y2(1:end-1)])
b_abs  = regress(y,[ones(size(y)) x1 Ttmp y2(1:end-1)]) % regression with absolute P values

% -- using observation-based coefficients instead
y  = sqrt(Qobs(2:end)'); % normalize runoff with squareroot
y2 = sqrt(Qobs(1:end-1)'); % normalize runoff with squareroot
x1 = Pobs_anom_rel(2:end)';
x2 = Tobs_anom(2:end)';
b_obs     = regress(y,[ones(size(y)) x1 x2 y2]);
b_obs(3)  = -0.091; % what if we instead substitute the Milly & Dunne (2020) T sensitivity of -9.1%/°C -- NOTE: only works if you do prediction in anomaly space
b_obs

% -- create a P/T based Qobs (predict Q based on P and T):
Qobs2 = ( b_obs(1) + b_obs(2)*x1 + b_obs(3)*x2 + b_obs(4)*y2 ).^2;
Qobs3 = qmod(Pobs_anom_rel,Tobs_anom,b_obs);


% -- estimate Q from climate models' P and T:
Q = NaN(length(models),max(ensmem),length(Padj(1,1,:)));
S = Q;
Sfrac = Q;
R = Q;
E = Q;
D = Q;
Def = Q;
for m = 1:length(models)
  for e = 1:ensmem(m)
    % -- using raw time series --
    Ptmp  = squeeze(Padj_anom_rel(m,e,:)); % P = variance- and mean-adjusted model precip
    Ttmp  = squeeze(Tanom(m,e,:));
    % -- calculate Q and reservoir terms
    Q(m,e,:) = qmod(Ptmp,Ttmp,b_obs);
    [S(m,e,:),Sfrac(m,e,:),R(m,e,:),E(m,e,:),C,D(m,e,:),Def(m,e,:)] = resmod(Q(m,e,:),Ptmp,Ttmp,Qobs);
    % -- using detrended T and P time series --
    Ptmp  = squeeze(Padj_anom_rel_dtr(m,e,:));
    Ttmp  = squeeze(Tanom_dtr(m,e,:));
    % -- calculate Q and reservoir terms
    Q2(m,e,:) = qmod(Ptmp,Ttmp,b_obs);
    [S2(m,e,:),Sfrac2(m,e,:),R2(m,e,:),E2(m,e,:),C2,D2(m,e,:),Def2(m,e,:)] = resmod(Q2(m,e,:),Ptmp,Ttmp,Qobs);
    % -- using raw T and detrended P time series --
    Ptmp  = squeeze(Padj_anom_rel_dtr(m,e,:));
    Ttmp  = squeeze(Tanom(m,e,:));
    % -- calculate Q and reservoir terms
    Q3(m,e,:) = qmod(Ptmp,Ttmp,b_obs);
    [S3(m,e,:),Sfrac3(m,e,:),R3(m,e,:),E3(m,e,:),C3,D3(m,e,:),Def3(m,e,:)] = resmod(Q3(m,e,:),Ptmp,Ttmp,Qobs);
    % -- using detrended T and raw P time series --
    Ptmp  = squeeze(Padj_anom_rel(m,e,:));
    Ttmp  = squeeze(Tanom_dtr(m,e,:));
    % -- calculate Q and reservoir terms
    Q4(m,e,:) = qmod(Ptmp,Ttmp,b_obs);
    [S4(m,e,:),Sfrac4(m,e,:),R4(m,e,:),E4(m,e,:),C4,D4(m,e,:),Def4(m,e,:)] = resmod(Q4(m,e,:),Ptmp,Ttmp,Qobs);
  end
end



if recalc1 == 1
  % -- calculation of risk surfaces (plotting follows later) --
  xlim = [-4 10];
  ylim = [-100 210];
  incr_x = .3; % box size over which to collect and average data (P, T, Sfrac)
  incr_y = 4;

  % -- calculate risk surfaces for various combinations:
  %    - raw T & P = risk1
  %    - detrended T & P = risk2
  %    - raw T & detrended P = risk3
  %    - expanded T & P = risk4
  xi = xlim(1):incr_x:xlim(2)+2*incr_x;
  yi = ylim(1):incr_y:ylim(2)+2*incr_y;
  clear('stor1','stor2','stor3','stor4','stor5')
  clear('risk1','risk2','risk3','risk4','risk5')
  risk_th = 0.2; % assess risk of falling below this storage fraction
  for i = 1:length(xi)-1
    ['' num2str(i) '/' num2str(length(xi)-1) '']
    for j = 1:length(yi)-1
      xv = [xi(i) xi(i) xi(i+1) xi(i+1) xi(i)];
      yv = [yi(j) yi(j+1) yi(j+1) yi(j) yi(j)];
      [in,on]     = inpolygon(Tanom(:),Padj_anom_rel(:),xv,yv);
      count1(i,j) = nansum(in);
      stor1(i,j)  = nanmedian(Sfrac(in));
      risk1(i,j)  = nansum(Sfrac(in)<risk_th)/count1(i,j);
      [in,on]     = inpolygon(Tanom_dtr(:),Padj_anom_rel_dtr(:),xv,yv);
      count2(i,j) = nansum(in);
      stor2(i,j)  = nanmedian(Sfrac2(in));
      risk2(i,j)  = nansum(Sfrac2(in)<risk_th)/count2(i,j);
      [in,on]     = inpolygon(Tanom(:),Padj_anom_rel_dtr(:),xv,yv);
      count3(i,j) = nansum(in);
      stor3(i,j)  = nanmedian(Sfrac3(in));
      risk3(i,j)  = nansum(Sfrac3(in)<risk_th)/count3(i,j);
      [in,on]     = inpolygon(Tanom_dtr(:),Padj_anom_rel(:),xv,yv);
      count4(i,j) = nansum(in);
      stor4(i,j)  = nanmedian(Sfrac4(in));
      risk4(i,j)  = nansum(Sfrac4(in)<risk_th)/count4(i,j);
    end
  end

  % -- smooth risk surfaces:
  swl = 2; % smoothing parameter
  risk1s = smooth2a(risk1,swl,swl);
  risk2s = smooth2a(risk2,swl,swl);
  risk3s = smooth2a(risk3,swl,swl);
  risk4s = smooth2a(risk4,swl,swl);
  stor1s = smooth2a(stor1,swl,swl);
  stor2s = smooth2a(stor2,swl,swl);
  stor3s = smooth2a(stor3,swl,swl);
  stor4s = smooth2a(stor4,swl,swl);
  % -- end of calculation of risk surfaces --
end


% -- calculate density of Storage traces over time:
incr = 0.1;
Sfrac_dens = NaN(length(Padj_anom(1,1,:)),length(0:incr:1-incr));
j = 1;
for i = 0:incr:1-(2*incr)
  tmp = squeeze(nansum(nansum(Sfrac>=i & Sfrac<(i+incr))));
  Sfrac_dens(:,j) = tmp./sum(ensmem);
  j = j+1;
end
% -- special case of Sfrac<=1
tmp = squeeze(nansum(nansum(Sfrac>=(i+incr) & Sfrac<=1)));
Sfrac_dens(:,j) = tmp./sum(ensmem);
% -- add extra row and column for pcolor plotting
Sfrac_dens(:,end+1) = NaN;
Sfrac_dens(end+1,:) = NaN;

Sfrac_dens(Sfrac_dens==0) = NaN;
irange = 0:incr:1; %-incr;


% -- Compound hot-dry events --
% -- using Bevacqua et al. 2022 definition:
%   - 10th percentile historical precip
%   - 90th percentile historical temperature
nl = 140;
for m = 1:length(models)
  tmp     = squeeze(Padj_anom_rel(m,:,1:50));
  P10(m)  = prctile(tmp(:),[10]);
  tmp     = squeeze(Tanom(m,:,1:50));
  T90(m)  = prctile(tmp(:),[90]);
end
% - set up empty matrices to store events (H=hot, D=dry, HD=hot-dry)
n1 = 3; % number of years to consider in epoch analysis BEFORE event
n2 = 5; % number of years to consider in epoch analysis AFTER event
ll = 2e3;
Sfrac_casesH  = NaN(length(models),ll,n1+n2+1);
Sfrac_casesD  = Sfrac_casesH;
Sfrac_casesHD = Sfrac_casesH;
for m = 1:length(models)
  countH  = 1;
  countD  = 1;
  countHD = 1;
  for e = 1:ensmem(m)
    casesH  = find(Tanom(m,e,1:nl-3)>T90(m) & Padj_anom_rel(m,e,1:nl-3)>=P10(m)); % hot cases that are not also dry
    casesD  = find(Padj_anom_rel(m,e,1:nl-3)<P10(m) & Tanom(m,e,1:nl-3)<=T90(m)); % dry cases that are not also hot
    casesHD = find(Tanom(m,e,1:nl-3)>T90(m) & Padj_anom_rel(m,e,1:nl-3)<P10(m)); % hot-dry cases
    % -- associated flow and storage changes
    for c = 1:length(casesH)
      if c > 3 && all(squeeze(Sfrac(m,e,casesH(c)-n1:casesH(c)+n2))>0)==1 % only cases that do not result in zero Storage
        Sfrac_casesH(m,countH,:)  = squeeze(Sfrac(m,e,casesH(c)-n1:casesH(c)+n2)); % select years +/- 3 years around year with extreme event (hot, dry, or hot-dry)
        Q_casesH(m,countH,:)      = squeeze(Q(m,e,casesH(c)-n1:casesH(c)+n2));
        countH = countH+1;
      end
    end
    for c = 1:length(casesD)
      if c > 3 && all(squeeze(Sfrac(m,e,casesD(c)-n1:casesD(c)+n2))>0)==1 % only cases that do not result in zero Storage
        Sfrac_casesD(m,countD,:)  = squeeze(Sfrac(m,e,casesD(c)-n1:casesD(c)+n2)); % select years +/- 3 years around year with extreme event (hot, dry, or hot-dry)
        Q_casesD(m,countD,:)      = squeeze(Q(m,e,casesD(c)-n1:casesD(c)+n2));
        countD = countD+1;
      end
    end
    for c = 1:length(casesHD)
      if c > 3 && all(squeeze(Sfrac(m,e,casesHD(c)-n1:casesHD(c)+n2))>0)==1 % only cases that do not result in zero Storage
        Sfrac_casesHD(m,countHD,:)  = squeeze(Sfrac(m,e,casesHD(c)-n1:casesHD(c)+n2)); % select years +/- 3 years around year with extreme event (hot, dry, or hot-dry)
        Q_casesHD(m,countHD,:)      = squeeze(Q(m,e,casesHD(c)-n1:casesHD(c)+n2));
        countHD = countHD+1;
      end
    end
  end
end







%% ================================================================================
% -- plotting
close all

cols = viridis(length(models));



if fig2c == 1

  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [40 10 20 10]) % on home screen
  hold on
  [h] = pcolor([time time(end)+1],irange,Sfrac_dens'*100);
  set(h, 'EdgeColor', 'none');
  % -- custom colormaps
  colorMap = [linspace(247,39,10)', linspace(238,117,10)', linspace(237,242,10)']/256; % light blue scale
  colorMap = [linspace(247,0,10)', linspace(238,0,10)', linspace(237,256,10)']/256; % blue scale
  % -- gray colormap
  colorMap = gray; colorMap = colorMap(end:-1:1,:);
  colormap(figure1,colorMap);
  a = colorbar;
  ylabel(a,'Density of all traces (%)')%,'FontSize',12)
  caxis([0 50])
  for m = 1:length(models)
    h2 = plot(time,squeeze(nanmean(Sfrac(m,:,:))),'Color',[.3 .3 .3],'LineWidth',1)
  end
  s = size(Sfrac);
  tmp = reshape(Sfrac,[s(1)*s(2),s(3)]);
  plot(time,prctile(tmp,[5 95],1),'k:','LineWidth',1)
  h0 = plot(1935:2023,mead_powell_combo(8:12:end,2)/max(mead_powell_combo(:,2)),'b','LineWidth',2) % actual Mead+Powell data (Capacity assumed to be the 1983 max level)
  h1 = plot(time_common,Sfracobs,'k','LineWidth',2)
  legend([h0 h1 h2],'Observed reservoir storage','Simulation with observations',['Simulations with SMILEs' char(10) '(only ensemble means shown)'],'location','SouthWest')
  legend boxoff
  ylabel('Storage fraction')
  xlabel('Time (Year)')
  set(gca,'Layer','top','YLim',[0 1],'XLim',[1950 2100])
  box on

  % tightfig

  return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = ['~/Dropbox/publication/lehner22_smile_perspective/fig/reservoir_model_smiles_traces_storage'];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  % return

end




if fig2d == 1

  close all

  dim = size(Sfrac_casesHD);

  % -- for each event, subtract the mean of the first pre-event years (e.g., first 3) from the whole time series (e.g., 7 years)
  %    to isolate the drop in reservoir storage as a function of hot, dry, or hot-dry events
  rangeH    = Sfrac_casesH - squeeze(nanmean(Sfrac_casesH(:,:,1:n1),3));
  rangeD    = Sfrac_casesD - squeeze(nanmean(Sfrac_casesD(:,:,1:n1),3));
  rangeHD   = Sfrac_casesHD - squeeze(nanmean(Sfrac_casesHD(:,:,1:n1),3));
  % -- create range for H+D from randomly sampling and summing individual H and D events
  clear('rangeHpD')
  for m = 1:length(models)
    tmp1  = squeeze(Sfrac_casesH(m,~isnan(Sfrac_casesH(m,:,1)),:)); % select non-NaN values
    tmp2  = squeeze(Sfrac_casesD(m,~isnan(Sfrac_casesD(m,:,1)),:)); % select non-NaN values
    tmp1  = tmp1 - squeeze(nanmean(tmp1(:,1:n1),2));
    tmp2  = tmp2 - squeeze(nanmean(tmp2(:,1:n1),2));
    nn1   = size(tmp1);
    nn1   = nn1(1);
    nn2   = size(tmp2);
    nn2   = nn2(1);
    nn    = max(nn1,nn2);
    iter  = 500;
    tmp4 = [];
    for i = 1:iter
      if nn1>nn2
        tmp3 = tmp1(randi(nn,nn2,1),:)+tmp2;
      end
      if nn1<nn2
        tmp3 = tmp1+tmp2(randi(nn,nn1,1),:);
      end
      tmp4 = [tmp4; tmp3];
    end
    rangeHpD(m,:,:) = tmp4(1:ll,:);
  end
  % -- test if difference between samples is significant
  [h,p]   = ttest2(squeeze(rangeHD(1,1:ll,:)),squeeze(rangeHpD(1,1:ll,:)),'alpha',.1)
  h(1:n1) = 0;
  h(h==0) = NaN;

  xj      = -n1:n2;

  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [40 10 10 11]) % on home screen
  hold on
  title(['Reservoir storage changes' char(10) 'during extreme events'],'FontSize',12)

  % -- shading uncertainty range
  p1 = 75;
  p2 = 25;
  tmp1 = reshape(rangeH,[length(models)*ll,n1+n2+1]);
  jbfill(xj,prctile(tmp1,p1,1),prctile(tmp1,p2,1),[1 .7 .7],'none',1,.5)
  tmp2 = reshape(rangeD,[length(models)*ll,n1+n2+1]);
  jbfill(xj,prctile(tmp2,p1,1),prctile(tmp2,p2,1),[.7 .7 1],'none',1,.5)
  tmp3 = reshape(rangeHD,[length(models)*ll,n1+n2+1]);
  jbfill(xj,prctile(tmp3,p1,1),prctile(tmp3,p2,1),[.7 1 .7],'none',1,.5)
  tmp4 = reshape(rangeHpD,[length(models)*ll,n1+n2+1]);
  jbfill(xj,prctile(tmp4,p1,1),prctile(tmp4,p2,1),[.7 1 .7],'none',1,.5)
  % --
  hold on
  mean_responseH    = squeeze(nanmean(nanmean(Sfrac_casesH,1),2))-squeeze(nanmean(nanmean(nanmean(Sfrac_casesH(:,:,1:n1)))));
  mean_responseD    = squeeze(nanmean(nanmean(Sfrac_casesD,1),2))-squeeze(nanmean(nanmean(nanmean(Sfrac_casesD(:,:,1:n1)))));
  mean_responseHD   = squeeze(nanmean(nanmean(Sfrac_casesHD,1),2))-squeeze(nanmean(nanmean(nanmean(Sfrac_casesHD(:,:,1:n1)))));
  mean_responseHpD  = squeeze(nanmean(nanmean(rangeHpD,1),2));
  h1 = plot(xj,mean_responseH,'r','LineWidth',3)
  h2 = plot(xj,mean_responseD,'b','LineWidth',3)
  h3 = plot(xj,mean_responseHD,'Color',[0 .5 0],'LineWidth',3)
  h4 = plot(xj,mean_responseHpD,':','Color',[0 .5 0],'LineWidth',3)
  plot(xj,mean_responseHD.*h','o','Color',[0 .5 0],'LineWidth',2,'MarkerSize',10)
  hline(0,'k')
  ax = gca(figure1);
  ax.XAxis.FontSize = 12;
  ax.YAxis.FontSize = 12;
  set(gca,'Layer','top','XLim',[-2 n2],'XTick',[-2:1:n2],'YLim',[-.4 .1])
  box on
  xlabel('Year relative to event year')%,'FontSize',12)
  ylabel('\DeltaStorage fraction')%,'FontSize',12)
  legend([h1 h2 h4 h3],'Hot events','Dry events','Hot events + dry events','Hot-dry events','location','SouthWest','FontSize',12)
  legend boxoff

  % tightfig

  return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = ['~/Dropbox/publication/lehner22_smile_perspective/fig/reservoir_model_smiles_hot_dry_event_composite'];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  return

end




if fig3 == 1
% -- model evaluation
  close all
  clear('f0','x0','x00')

  % -- (a) and (b): Formal tests of whether distributions are different --
  for m = 1:length(models)
    for e = 1:ensmem(m)
      kstest_T(m,e) = kstest2(Tobs_anom_dtr,squeeze(Tanom_dtr(m,e,1:length(Tobs_anom))),'Alpha',.05);
      kstest_P(m,e) = kstest2(Pobs_anom_dtr,squeeze(Praw_dtr(m,e,1:length(Pobs_anom_rel))),'Alpha',.05);
    end
  end
  % -- Fraction of ensemble members with significantly different Tanom or Panom distributions from obs
  for m = 1:length(models)
    Tanom_test(m) = nansum(kstest_T(m,:),2)/ensmem(m);
    Panom_test(m) = nansum(kstest_P(m,:),2)/ensmem(m);
  end

  % -- (c): Correlation between P and T over time (correlation within each model)
  wl = length(Tobs_anom); % 72 30 50
  for m = 1:length(models)
    for e = 1:ensmem(m)
      for i = 1:length(Padj_anom_rel)-wl
        % rPT(m,e,i) = corr(squeeze(P(m,e,i:i+wl-1)),squeeze(T(m,e,i:i+wl-1)));
        rPT(m,e,i) = corr(squeeze(Padj_anom_rel_dtr(m,e,i:i+wl-1)),squeeze(Tanom_dtr(m,e,i:i+wl-1)));
      end
    end
  end
  % -- correlation between deltaP and deltaT (across models! not something where we can compare with obs...)
  for e = 1:20
    dP = squeeze(nanmean(Padj_anom_rel(:,e,end-wl+1:end))-nanmean(Padj_anom_rel(:,e,1:wl)));
    dT = squeeze(nanmean(Tanom(:,e,end-wl+1:end))-nanmean(Tanom(:,e,1:wl)));
    corr_dP_dT(e) = corr(dP,dT);
  end

  % -- (d) and (e): Rank of observed trend within model ensemble of trends
  trend_length = length(refstart:time_Tobs(end));
  [tmp,ci] = regress(Tobs_anom',[ones(size(Tobs_anom')) (refstart:time_Tobs(end))'], (1-0.834))
  Tobs_anom_trend    = tmp(2) * trend_length;
  Tobs_anom_trend_ci = ci(2,:) .* trend_length;
  [tmp,ci] = regress(Pobs',[ones(size(Pobs')) (refstart:time_Pobs(end))'], (1-0.834))
  Pobs_trend    = tmp(2) * trend_length;
  Pobs_trend_ci = ci(2,:) .* trend_length;
  clear('Ttrend','Ptrend')
  for m = 1:length(models)
    for e = 1:ensmem(m)
      tmp = polyfit(refstart:time_Tobs(end),Traw(m,e,1:length(Tobs_anom)),1);
      Ttrend(m,e) = tmp(1) * trend_length;
      tmp = polyfit(refstart:time_Pobs(end),Praw(m,e,1:length(Pobs_anom_rel)),1);
      Ptrend(m,e) = tmp(1) * trend_length;
    end
    % -- check whether any model trend is within 83.4% confidence interval of obs trend
    if sum(floor(((Ttrend(m,:) >= min(Tobs_anom_trend_ci)) + (Ttrend(m,:) <= max(Tobs_anom_trend_ci))) / 2)) > 0
      Ttrend_rank(m) = 1;
    else
      Ttrend_rank(m) = 0;
    end
    if sum(floor(((Ptrend(m,:) >= min(Pobs_trend_ci)) + (Ptrend(m,:) <= max(Pobs_trend_ci))) / 2)) > 0
      Ptrend_rank(m) = 1;
    else
      Ptrend_rank(m) = 0;
    end
  end
  % -- Does observed T-P correlation fall outside 95% of ensemble members?
  rPT_obs = corr(Pobs_anom_rel_dtr',Tobs_anom_dtr');
  for m = 1:length(models)
    (sum(squeeze(rPT(m,:,1))>rPT_obs)/ensmem(m))
    if (sum(squeeze(rPT(m,:,1))>rPT_obs)/ensmem(m)) > 0.95 || (sum(squeeze(rPT(m,:,1))<rPT_obs)/ensmem(m)) > 0.95
      corr_test(m) = 1; % bad
    else
      corr_test(m) = 0; % good
    end
  end

  % -- print for sanity check
  % (Tanom_test>0.05)'
  % (Panom_test>0.05)'
  % corr_test'
  % (Tanom_test>0.05)' + (Panom_test>0.05)' + corr_test'


  % --------------------------------------------------------
  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [40 10 30 15])

  subplot(2,3,1)
  hold on
  title('(a) Temperature distribution')
  for m = 1:length(models)
    for e = 1:ensmem(m)
      tmp = detrend(squeeze(Traw(m,e,1:length(Tobs_anom))));
      [f,x,u] = ksdensity(tmp(:));
      f00(e,:) = f;
      x00(e,:) = x;
      line(x,f,'Color',[1 1 1]*(1-(m/40)),'LineWidth',.5)
    end
    f0(m,:) = nanmean(f00);
    x0(m,:) = nanmean(x00);
  end
  for m = 1:length(models)
    if Tanom_test(m) >= 0.05
      h2 = line(x0(m,:),f0(m,:),'Color',[.4 .4 .4], 'LineWidth',2,'LineStyle',':');
    else
      h2 = line(x0(m,:),f0(m,:),'Color',[.4 .4 .4], 'LineWidth',2,'LineStyle',':');
      h1 = line(x0(m,:),f0(m,:),'Color',[.4 .4 .4], 'LineWidth',2);
    end
  end
  [f0,x0,u] = ksdensity(Tobs_anom_dtr);
  h0 = line(x0,f0,'Color',[0 0 1], 'LineWidth',2);
  set(gca,'ytick',[])
  xlabel('Temperature anomaly (\circC)')
  box on
  legend([h0 h1 h2],'Observations','Model (good)','Model (bad)','location','northwest','FontSize',9)


  subplot(2,3,2)
  hold on
  title('(b) Precipitation distribution')
  for m = 1:length(models)
    for e = 1:ensmem(m)
      tmp = detrend(squeeze(Praw(m,e,1:length(Pobs_anom_rel))));
      [f,x,u] = ksdensity(tmp(:));
      f00(e,:) = f;
      x00(e,:) = x;
      line(x,f,'Color',[1 1 1]*(1-(m/40)),'LineWidth',.5)
    end
    f0(m,:) = nanmean(f00);
    x0(m,:) = nanmean(x00);
  end
  for m = 1:length(models)
    if Panom_test(m) >= 0.05
      h2 = line(x0(m,:),f0(m,:),'Color',[.4 .4 .4], 'LineWidth',2,'LineStyle',':');
    else
      h2 = line(x0(m,:),f0(m,:),'Color',[.4 .4 .4], 'LineWidth',2,'LineStyle',':');
      h1 = line(x0(m,:),f0(m,:),'Color',[.4 .4 .4], 'LineWidth',2);
    end
  end
  [f0,x0,u] = ksdensity(Pobs_anom_dtr);
  h1 = line(x0,f0,'Color',[0 0 1], 'LineWidth',2);
  set(gca,'ytick',[])
  xlabel('Precipitation anomaly (mm/month)')
  box on


  subplot(2,3,4)
  hold on
  title('(d) Temperature trends')
  clear('Tchange_global','Tchange','Pchange')
  yl = [3 8];
  jbfill([Tobs_anom_trend_ci(1) Tobs_anom_trend_ci(2)],[yl(2) yl(2)],[yl(1) yl(1)],[.7 .7 1],'none',1,.5)
  hold on
  h0 = plot([Tobs_anom_trend Tobs_anom_trend],[yl(1) yl(2)],'b','LineWidth',2)
  for m = 1:length(models)
    Tchange_global(m,:) = nanmean(Tanom2_global(m,:,100:149),3)-nanmean(Tanom2_global(m,:,1:50),3);
    Tchange(m,:) = (nanmean(Tanom2(m,:,100:149),3)-nanmean(Tanom2(m,:,1:50),3)); %./Tchange2_global(m,:);
    plot(Ttrend(m,1:ensmem2(m)),Tchange(m,1:ensmem2(m)),'o','Color',[.7 .7 .7])
    if Ttrend_rank(m) == 0
      h2 = plot_ellipse(Ttrend(m,1:ensmem2(m))',Tchange(m,1:ensmem2(m))')
      h2.Color = [.5 .5 .5]
      h2.LineStyle = ':'
      h2.LineWidth = 1.5
    else
      h1 = plot_ellipse(Ttrend(m,1:ensmem2(m))',Tchange(m,1:ensmem2(m))')
      h1.Color = [.5 .5 .5]
      h1.LineStyle = '-'
      h1.LineWidth = 1.5
    end
  end
  xl = xlim;
  text(xl(1)+(xl(2)-xl(1))*0.7, yl(1)+(yl(2)-yl(1))*0.1, ['R^2 = ' num2str(round(corr(nanmean(Ttrend,2),nanmean(Tchange,2))^2,2))],'FontSize',10)
  set(gca,'Layer','top','YLim',yl,'XLim',xl)
  xlabel(['Historical T trend, ' num2str(refstart) '-' num2str(time_Tobs(end)) ' (\circC/' num2str(trend_length) 'yr)'])
  ylabel(['Future T change (\circC)' char(10) '(2050-2099)-(1950-1999)'])
  box on
  legend([h0 h1 h2],'Observations','Model (good)','Model (bad)','location','northwest')%,'FontSize',10)

  subplot(2,3,5)
  hold on
  title('(e) Precipitation trends')
  yl = [-25 40];
  xl = [-130 110];
  set(gca,'Layer','top','YLim',yl,'XLim',xl)
  hline(0,'k')
  vline(0,'k')
  jbfill([Pobs_trend_ci(1) Pobs_trend_ci(2)],[yl(2) yl(2)],[yl(1) yl(1)],[.7 .7 1],'none',1,.5)
  hold on
  h0 = plot([Pobs_trend Pobs_trend],[yl(1) yl(2)],'b','LineWidth',2)
  for m = 1:length(models)
    Pchange(m,:) = (nanmean(Padj2_anom(m,:,100:149),3)-nanmean(Padj2_anom(m,:,1:50),3))./Tchange(m,:);
    plot(Ptrend(m,1:ensmem2(m)),Pchange(m,1:ensmem2(m)),'o','Color',[.7 .7 .7])
    if Ptrend_rank(m) == 0
      h2 = plot_ellipse(Ptrend(m,1:ensmem2(m))',Pchange(m,1:ensmem2(m))')
      h2.Color = [.5 .5 .5]
      h2.LineStyle = ':'
      h2.LineWidth = 1.5
    else
      h1 = plot_ellipse(Ptrend(m,1:ensmem2(m))',Pchange(m,1:ensmem2(m))')
      h1.Color = [.5 .5 .5]
      h1.LineStyle = '-'
      h1.LineWidth = 1.5
    end
  end
  xl = xlim;
  text(xl(1)+(xl(2)-xl(1))*0.7, yl(1)+(yl(2)-yl(1))*0.1, ['R^2 = ' num2str(round(corr(nanmean(Ptrend,2),nanmean(Pchange,2))^2,2))],'FontSize',10)
  set(gca,'Layer','top','YLim',yl,'XLim',xl)
  xlabel(['Historical P trend, ' num2str(refstart) '-' num2str(time_Tobs(end)) ' (mm/month/' num2str(trend_length) 'yr)'])
  ylabel(['Future P change (mm/month/\circC)' char(10) '(2050-2099)-(1950-1999)'])
  box on

  subplot(2,3,3)
  title(['(c) Correlation of detrended' char(10) 'temperature and precipitation'])
  for m = 1:length(models)
    jbfill(refstart+ceil(wl/2):2099-floor(wl/2)-1,prctile(squeeze(rPT(m,:,:)),95),prctile(squeeze(rPT(m,:,:)),5),[1 1 1]*(1-(m/30)),'none',1,.5)
    hold on
  end
  for m = 1:length(models)
    if corr_test(m)==1
      h2 = plot(refstart+ceil(wl/2):2099-floor(wl/2)-1, squeeze(nanmean(rPT(m,:,:),2))','Color',[.4 .4 .4],'LineWidth',2,'LineStyle',':')
    else
      h1 = plot(refstart+ceil(wl/2):2099-floor(wl/2)-1, squeeze(nanmean(rPT(m,:,:),2))','Color',[.4 .4 .4],'LineWidth',2)
    end
  end
  h0 = plot([refstart+ceil(wl/2) refstart+ceil(wl/2)],[rPT_obs rPT_obs],'b.','MarkerSize',40)
  set(gca,'XLim',[refstart+ceil(wl/2)-5 2099-floor(wl/2)+5],'YLim',[-.8 .3])%, 'YAxisLocation','right')%,'Color','none','XColor','b','YColor','b');
  box on
  xlabel('Time (Year)')
  ylabel('Correlation coefficient')
  legend([h0 h1 h2],'Observations','Model (good)','Model (bad)','location','northeast','FontSize',9)
  legend boxoff


  subplot(2,3,6)
  title(['(f) Summary'])
  hold on
  all = [Tanom_test<0.05; Panom_test<0.05; corr_test==0; Ttrend_rank; Ptrend_rank];
  all = all(end:-1:1,:);
  all(:,end+1) = NaN;
  all(end+1,:) = NaN;
  pcolor(all)
  colorMap = [232, 60, 66; 21, 140, 39]/256;
  colormap(figure1,colorMap);
  set(gca,'XLim',[1 length(models)+1])
  set(gca,'yticklabel',[],'xtick',[1:1:length(models)]+.5,'xticklabel',[1:1:length(models)])
  xtickangle(45)
  xlabel('Model')
  xx = -.4;
  text(xx, 1.5, '(e)')
  text(xx, 2.5, '(d)')
  text(xx, 3.5, '(c)')
  text(xx, 4.5, '(b)')
  text(xx, 5.5, '(a)')
  text(xx-.2, 6.3, 'Test:')


  tightfig

  return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = ['~/Dropbox/publication/lehner22_smile_perspective/fig/model_validation_v2'];
  print('-painters','-dpdf', ['' fileo '.pdf'])
  % return

end




if fig4 == 1
  close all

  wl = 50; % length of climatology in years

  xlim  = [-4 9.5];
  ylim  = [130 900];
  ylim2 = [-75 110];

  figure2 = figure;
  set(figure2, 'units', 'centimeters', 'pos', [40 10 45 8]) % on home screen

  subplot(1,4,1)
  title(['(a) Actual T, actual P'])
  hold on
  h = pcolor(xi(1:end-1),yi(1:end-1)-100, risk1s');
  set(h, 'EdgeColor', 'none');
  colorMap = [linspace(247,39,10)', linspace(238,117,10)', linspace(237,242,10)']/256; % 247, 238, 237 / 39, 117, 242
  colorMap = othercolor('BuOr_10',10);
  colormap(colorMap);
  caxis([0 1])
  colorbar
  % -- multivariate ksdensity
  x1 = xlim(1):.5:xlim(2);
  y1 = ylim(1):incr_y*2:ylim(2);
  y1 = ((y1/nanmean(Pobs(1:wl)))-1)*100;
  [gridx1,gridy1] = ndgrid(x1,y1);
  gridx1 = gridx1';
  gridy1 = gridy1';
  xy = [gridx1(:) gridy1(:)];
  bw = .7; % bandwidth parameter for ksdensity
  cq = .95; % contour quantile (e.g., .9=90% of points within contour)
  x         = Tanom(:,:,1:wl);
  y         = Padj_anom_rel(:,:,1:wl)-100;
  h = plot_ellipse(x(:),y(:))
  h.Color = 'k'
  h.LineWidth = 3
  h.LineStyle = '-'
  x         = Tanom(:,:,end-wl+1:end);
  y         = Padj_anom_rel(:,:,end-wl+1:end)-100;
  h = plot_ellipse(x(:),y(:))
  h.Color = 'k'
  h.LineWidth = 3
  h.LineStyle = ':'
  set(gca,'Layer','top','YLim',ylim2,'XLim',xlim,'XTick',[-4:2:10])
  hline(0,'k')
  vline(0,'k')
  % -- plot model means into climatology space
  wl = 50; % running mean window length
  for m = 1:length(models)
    tmp1 = rm(squeeze(nanmean(Tanom(m,:,:))),wl);
    tmp2 = rm(squeeze(nanmean(Padj_anom_rel(m,:,:)-100)),wl);
    h3 = plot([0,tmp1(end-wl/2)],[0,tmp2(end-wl/2)],'k.-','LineWidth',1,'MarkerSize',16)
  end
  % -- plot obs years
  split_obs = 0; % 0=no, 1=yes
  if split_obs == 1
    h1 = plot(Tobs_anom(1:wl)-nanmean(Tobs_anom(1:wl)),((Pobs_anom_rel(1:wl)/nanmean(Pobs_anom_rel(1:wl)))-1)*100,'o','Color',[.7 .7 .7],'MarkerSize',3,'LineWidth',1)
    h2 = plot(Tobs_anom(wl+1:end)-nanmean(Tobs_anom(1:wl)),((Pobs_anom_rel(wl+1:end)/nanmean(Pobs_anom_rel(1:wl)))-1)*100,'o','Color',[1 .3 .3],'MarkerSize',3,'LineWidth',1)
    legend([h1(1) h2(1) h3(1)],'Observations 1950-1999',['Observations 2000-' num2str(time_Tobs(end))],'SMILE ensemble means','location','best')
  else
    h1 = plot(Tobs_anom-nanmean(Tobs_anom),((Pobs_anom_rel/nanmean(Pobs_anom_rel))-1)*100,'o','Color',[.7 .7 .7],'MarkerSize',3,'LineWidth',1)
    legend([h1(1) h3(1)],['Observations 1950-' num2str(time_Tobs(end))],'SMILE ensemble means','location','best')
  end
  legend boxoff
  box on
  ylabel(['\DeltaP (%)'])
  xlabel(['\DeltaT (\circC)'])


  subplot(1,4,2)
  title(['(b) Actual T, detrended P'])
  hold on
  h = pcolor(xi(1:end-1),yi(1:end-1)-100, risk3s');
  set(h, 'EdgeColor', 'none');
  colorMap = [linspace(247,39,10)', linspace(238,117,10)', linspace(237,242,10)']/256; % 247, 238, 237 / 39, 117, 242
  colorMap = othercolor('BuOr_10',10);
  colormap(colorMap);
  caxis([0 1])
  colorbar% -- multivariate ksdensity
  [gridx1,gridy1] = ndgrid(x1,y1);
  gridx1 = gridx1';
  gridy1 = gridy1';
  xy = [gridx1(:) gridy1(:)];
  x         = Tanom(:,:,1:wl);
  y         = Padj_anom_rel_dtr(:,:,1:wl)-100;
  h = plot_ellipse(x(:),y(:))
  h.Color = 'k'
  h.LineWidth = 3
  h.LineStyle = '-'
  x         = Tanom(:,:,end-wl+1:end);
  y         = Padj_anom_rel_dtr(:,:,end-wl+1:end)-100;
  h = plot_ellipse(x(:),y(:))
  h.Color = 'k'
  h.LineWidth = 3
  h.LineStyle = ':'
  set(gca,'YLim',ylim2,'XLim',xlim,'XTick',[-4:2:10])
  hline(0,'k')
  vline(0,'k')
  box on
  set(gca,'Layer','top')
  ylabel(['\DeltaP (%)'])
  xlabel(['\DeltaT (\circC)'])

  subplot(1,4,3)
  title(['(c) Detrended T, actual P'])
  hold on
  h = pcolor(xi(1:end-1),yi(1:end-1)-100, risk4s');
  set(h, 'EdgeColor', 'none');
  colorMap = [linspace(247,39,10)', linspace(238,117,10)', linspace(237,242,10)']/256; % 247, 238, 237 / 39, 117, 242
  colorMap = othercolor('BuOr_10',10);
  colormap(colorMap);
  caxis([0 1])
  colorbar% -- multivariate ksdensity
  [gridx1,gridy1] = ndgrid(x1,y1);
  gridx1 = gridx1';
  gridy1 = gridy1';
  xy = [gridx1(:) gridy1(:)];
  x         = Tanom_dtr(:,:,1:wl);
  y         = Padj_anom_rel(:,:,1:wl)-100;
  h = plot_ellipse(x(:),y(:))
  h.Color = 'k'
  h.LineWidth = 3
  h.LineStyle = '-'
  x         = Tanom_dtr(:,:,end-wl+1:end);
  y         = Padj_anom_rel(:,:,end-wl+1:end)-100;
  h = plot_ellipse(x(:),y(:))
  h.Color = 'k'
  h.LineWidth = 3
  h.LineStyle = ':'
  set(gca,'YLim',ylim2,'XLim',xlim,'XTick',[-4:2:10])
  hline(0,'k')
  vline(0,'k')
  box on
  set(gca,'Layer','top')
  ylabel(['\DeltaP (%)'])
  xlabel(['\DeltaT (\circC)'])

  subplot(1,4,4)
  title(['(d) Detrended T, detrended P'])
  hold on
  h = pcolor(xi(1:end-1),yi(1:end-1)-100, risk2s');
  set(h, 'EdgeColor', 'none');
  caxis([0 1])
  colorbar
  % -- multivariate ksdensity
  x1 = xlim(1):.5:xlim(2);
  y1 = ylim(1):incr_y*2:ylim(2);
  y1 = ((y1/nanmean(Pobs(1:wl)))-1)*100;
  [gridx1,gridy1] = ndgrid(x1,y1);
  gridx1 = gridx1';
  gridy1 = gridy1';
  xy = [gridx1(:) gridy1(:)];
  x         = Tanom_dtr(:,:,1:wl);
  y         = Padj_anom_rel_dtr(:,:,1:wl)-100;
  h = plot_ellipse(x(:),y(:))
  h.Color = 'k'
  h.LineWidth = 3
  h.LineStyle = '-'
  x         = Tanom_dtr(:,:,end-wl+1:end);
  y         = Padj_anom_rel_dtr(:,:,end-wl+1:end)-100;
  h = plot_ellipse(x(:),y(:))
  h.Color = 'k'
  h.LineWidth = 3
  h.LineStyle = ':'
  set(gca,'YLim',ylim2,'XLim',xlim,'XTick',[-4:2:10])
  hline(0,'k')
  vline(0,'k')
  box on
  ylabel(['\DeltaP (%)'])
  xlabel(['\DeltaT (\circC)'])

  % tightfig

  return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = ['~/Dropbox/publication/lehner22_smile_perspective/fig/reservoir_model_smiles_risk_surfaces_ellipses'];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  return

end





if figS2 == 1

  close all

  m = 6; % MPI-ESM
  % ---
  vari  = 'tas';
  in    = dir([pathin vari '_Amon_' models{m} '_*ens_g025_UCRB_ts_anom.nc']);
  tmp0	= squeeze(ncread([in.folder '/' in.name],vari));
  for e = 1:100
    tmp = wy_mean(tmp0(e,:));
    tmp = tmp(refstart-start(m)+1:refende-start(m)); % cut to common length
    Traw_mpi(e,:) = tmp;
    Tanom_mpi(e,:)= tmp-nanmean(tmp(1:time_Tobs(end)-refstart)); % reference to 1950-2023
  end
  % ---
  vari  = 'pr';
  in    = dir([pathin vari '_Amon_' models{m} '_*ens_g025_UCRB_ts_anom.nc']);
  tmp0	= squeeze(ncread([in.folder '/' in.name],vari)) * tf;
  for e = 1:100
    tmp                   = wy(tmp0(e,:));
    tmp                   = tmp(refstart-start(m)+1:refende-start(m)); % cut to common length
    Praw_mpi(e,:)           = tmp;
    Praw_mpi_anom(e,:)      = tmp-nanmean(tmp(1:time_Pobs(end)-refstart)); % absolute anomalies to 1950-2023 mean
    Padj_mpi(e,:)           = Praw_mpi_anom(e,:) * (std(Pobs)/std(tmp(1:time_Pobs(end)-refstart))) + mean(Pobs); % adjust variance and mean to match observations
    Padj_mpi_anom(e,:)      = Padj_mpi(e,:)-nanmean(Padj_mpi(e,1:time_Pobs(end)-refstart)); % absolute anomalies to 1950-2023 mean
    Padj_mpi_anom_rel(e,:)  = (Padj_mpi(e,:)/nanmean(Padj_mpi(e,1:time_Pobs(end)-refstart)))*100-100; % relative anomalies to 1950-2023 mean
  end


  close all
  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [40 10 20 15])

  nn = 20; % which ensemble size to highlight with vertical dotted line?
  iter = 500; % how many times to resample the full ensemble to estimate error?

  data = Tanom_mpi;
  truth1 = nanmean(data,1); % true ensemble mean
  truth2 = nanstd(data,1); % true stddev
  truth3 = nanmean(truth1(100:149))-nanmean(truth1(1:50)); % true mean change signal
  tmp    = polyfit(time,truth2,1); % true stddev change signal
  truth4 = tmp(1)*length(time); % true stddev change signal
  count1 = 1;
  clear('tmp0','tmp1','tmp2')
  clear('pseudo_truth1','pseudo_truth2','pseudo_truth3','pseudo_truth4')
  for nm = 2:100
    for i = 1:iter
      tmp0 = data(randperm(100,nm),:);
      pseudo_truth1(count1,i,:) = nanmean(tmp0,1);
      tmp1(i) = rmse(truth1,squeeze(pseudo_truth1(count1,i,:))');
      pseudo_truth2(count1,i,:) = nanstd(tmp0,1);
      tmp2(i) = rmse(truth2,squeeze(pseudo_truth2(count1,i,:))');
      pseudo_truth3(count1,i) = nanmean(squeeze(pseudo_truth1(count1,i,100:149)))-nanmean(squeeze(pseudo_truth1(count1,i,1:50)));
      tmp    = polyfit(time,pseudo_truth2(count1,i,:),1);
      pseudo_truth4(count1,i) = tmp(1)*length(time);
    end
    data_rmse(count1) = nanmean(tmp1);
    std_rmse(count1)  = nanmean(tmp2);
    count1 = count1+1;
  end

  ensmem_ext  = [40,50,30,16,20,30,100,100,40,25,21,11,50,30,30,20,30]; % extended ensemble sizes
  nums(1) = 29; % first entry: number of CMIP6 models
  for e = 10:10:100
    nums(e/10+1) = sum(ensmem_ext>=e); % number of models with a certain ensemble size
  end


  subplot(2,2,1)
  title(['(a) Error in estimating' char(10) 'T \mu change'])
  hold on
  % -- y-axis 1
  plot(2:100,prctile(((pseudo_truth3-truth3)/truth3)*100,[95],2),'k')
  ylabel('Potential error (%)')%,'FontSize',12)
  % -- y-axis 2
  yyaxis right
  plot([1 10:10:100],nums,'b.-')
  ylabel('# of models available')%,'FontSize',12)
  xlabel('Ensemble size')
  set(gca,'XLim',[1 100],'YColor','b')
  vline(nn,'k:')
  box on

  subplot(2,2,2)
  title(['(b) Error in estimating' char(10) 'T \sigma change'])
  hold on
  % -- y-axis 1
  plot(2:100,prctile(((pseudo_truth4-truth4)/truth4)*100,[95],2),'k')
  ylabel('Potential error (%)')%,'FontSize',12)
  % -- y-axis 2
  yyaxis right
  plot([1 10:10:100],nums,'b.-')
  ylabel('# of models available')%,'FontSize',12)
  xlabel('Ensemble size')
  set(gca,'XLim',[1 100],'YColor','b')
  vline(nn,'k:')
  box on


  data = Padj_mpi_anom_rel;
  truth1 = nanmean(data,1); % true ensemble mean
  truth2 = nanstd(data,1); % true stddev
  truth3 = nanmean(truth1(100:149))-nanmean(truth1(1:50)); % true mean change signal
  tmp    = polyfit(time,truth2,1); % true stddev change signal
  truth4 = tmp(1)*length(time); % true stddev change signal
  count1 = 1;
  clear('tmp0','tmp1','tmp2')
  clear('pseudo_truth1','pseudo_truth2','pseudo_truth3','pseudo_truth4')
  for nm = 2:100
    for i = 1:iter
      tmp0 = data(randperm(100,nm),:);
      pseudo_truth1(count1,i,:) = nanmean(tmp0,1);
      tmp1(i) = rmse(truth1,squeeze(pseudo_truth1(count1,i,:))');
      pseudo_truth2(count1,i,:) = nanstd(tmp0,1);
      tmp2(i) = rmse(truth2,squeeze(pseudo_truth2(count1,i,:))');
      pseudo_truth3(count1,i) = nanmean(squeeze(pseudo_truth1(count1,i,100:149)))-nanmean(squeeze(pseudo_truth1(count1,i,1:50)));
      tmp    = polyfit(time,pseudo_truth2(count1,i,:),1);
      pseudo_truth4(count1,i) = tmp(1)*length(time);
    end
    data_rmse(count1) = nanmean(tmp1);
    std_rmse(count1)  = nanmean(tmp2);
    count1 = count1+1;
  end

  subplot(2,2,3)
  title(['(c) Error in estimating' char(10) 'P \mu change'])
  hold on
  % -- y-axis 1
  plot(2:100,prctile(((pseudo_truth3-truth3)/truth3)*100,[95],2),'k')
  ylabel('Potential error (%)')%,'FontSize',12)
  % -- y-axis 2
  yyaxis right
  plot([1 10:10:100],nums,'b.-')
  ylabel('# of models available')%,'FontSize',12)
  xlabel('Ensemble size')
  set(gca,'XLim',[1 100],'YColor','b')
  vline(nn,'k:')
  box on

  subplot(2,2,4)
  title(['(d) Error in estimating' char(10) 'P \sigma change'])
  hold on
  % -- y-axis 1
  plot(2:100,prctile(((pseudo_truth4-truth4)/truth4)*100,[95],2),'k')
  ylabel('Potential error (%)')%,'FontSize',12)
  % -- y-axis 2
  yyaxis right
  plot([1 10:10:100],nums,'b.-')
  ylabel('# of models available')%,'FontSize',12)
  xlabel('Ensemble size')
  set(gca,'XLim',[1 100],'YColor','b')
  vline(nn,'k:')
  box on

  tightfig

  return
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'renderer','Painters')
  fileo = ['~/Dropbox/publication/lehner22_smile_perspective/fig/ensemble_size_test4'];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  save2pdf(['' fileo '.pdf'])
  % return

end



% ===============================================================================================



% -- functions -------------
function [S,Sfrac,R,E,C,D,Def] = resmod(Q,P,T,Qobs)
  % -- reservoir model
  % -- takes absolute Q, relative P, anomalous T, and absolute Qobs
  % --------------------------------------------
  C   = 4*mean(Qobs); % reservoir capacity, as a multiple of OBSERVED annual average flow
  n   = length(Q); % length of data in years
  % -- Demand:
  Dmax = 1;
  D   = [linspace(.85,Dmax,50) linspace(Dmax,.95,23) ones(1,length(P)-73)*.95] * mean(Q(1:length(Qobs))); % mimic the real-world ramp-up in demand on the Colorado River - need to select
  % --
  S   = .35*C;  % initial storage proportion (start with 35% full)
  Sfrac = S/C;
  E   = (0.025*T(1)+0.08)*(184/28.6); % https://edepot.wur.nl/183107 % initial evaporation
  for i = 2:n
    % -- evaporation from the reservoir (set to fraction of storage initially but also temperature-dependent)
    E(i)  = (0.025*T(i-1)+0.08)*(184/28.6); % https://edepot.wur.nl/183107
    % -- demand D dependent on predifined ramp-up of D and variations in T and P
    D(i)  = D(i) + (D(i)*((T(i)-nanmean(T(1:length(Qobs))))/50)*0) - (D(i)*((P(i)/nanmean(P(1:50)))-1)*0.1); % makes demand dependent on T and P; scaled to only modify ramped-up D defined above
    % --
    R(i)  = max( min(D(i),S(i-1)+Q(i)-E(i)) , (S(i-1)+Q(i)-E(i))-C); % Release; can either be exactly the Demand or Storage+Inflow Q (if there's not enough water in S+Q)
    S(i)  = S(i-1) + Q(i) - R(i) - E(i); % new Storage = current Storage+Inflow(Q)-Release-Evaporation
    Sfrac(i) = S(i)/C; % Storage as fraction of Capacity
    Def(i)  = R(i)-D(i); % subtract Demand from Release to get Deficit
  end
end


% -- runoff model
function Qhat = qmod(P,T,b) % takes relative P (in % of normal) and anomalous T
  x1 = P;
  x2 = T;
  Qhat(1) = b(1) + b(2)*x1(1) + b(3)*x2(1) + b(4)*sqrt(16); % hard-coded first year of prediction, using climatological Q (16M acre-feet) for preceding year
  for i = 2:length(x1)
    Qhat(i) = b(1) + b(2)*x1(i) + b(3)*x2(i) + b(4)*Qhat(i-1); % use previous year Q in this year's Q calculation
  end
  Qhat = Qhat.^2; % transfer back to real space (Q was scaled via square root during the regression model development)
end
