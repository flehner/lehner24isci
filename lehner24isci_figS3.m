close all
clear all

% ---------------------------
% lehner24isci_figS3.m
% ---------------------------
% Code to produce Fig. S3 in Lehner (2024)
%
% Citation:
% Lehner, F. (2024):
% Climate model large ensembles as test beds for applied compound event research
% iScience, DOI: TBD
%
% Notes on code:
% - Only does calculations and plotting, not pre-processing
% - Requires pre-processed data, which is provided in the 'data' subfolder:
%   - Temperature and precipitation from SMILEs an CMIP6 (CMIP5+6), averaged over the Upper Colorado River Basin
% ---------------------------


pathin = '/Users/fl439/Dropbox/work/reservoir_model/';

refstart  = 1950;
refende   = 2099;
time      = refstart+1:refende; % +1 because of water year notation

time_Tobs	= 1896:2023;
time_Pobs	= 1896:2023;

% == SMILEs data =================================================================
model_scen_pairs = ...
{'CESM1-CAM5','rcp85';...
'CanESM2','rcp85';...
'CSIRO-Mk3-6-0','rcp85';...
'EC-EARTH','rcp85';...
'GFDL-CM3','rcp85';...
'GFDL-ESM2M','rcp85';...
'MPI-ESM','rcp85';...
'CESM2','ssp370';...
'ACCESS-ESM1-5','ssp245';...
'CanESM5','ssp370';...
'EC-Earth3','ssp245';...
'IPSL-CM6A-LR','ssp245';...
'MIROC6','ssp126';...
'MIROC-ES2L','ssp245';...
'GFDL-SPEAR-MED','ssp585';...
'E3SMv2','ssp370';...
'MPI-ESM1-2-LR','ssp126'};
models      = cellstr(model_scen_pairs(:,1));
scen_smiles = cellstr(model_scen_pairs(:,2));
models_select = [1:3 5:11 13:14 16:17]; %1:length(models)
models      = models(models_select);
scen_smiles = scen_smiles(models_select);
ensmem  = [40,50,30,16,20,30,100,100,40,25,21,11,50,30,30,20,30];
ensmem  = zeros(1,length(models)) + min(ensmem(models_select)); % limit all ensemble sizes to e.g. 16, 20, 100
start   = [1920,1950,1850,1860,1920,1950,1850,1850,1850,1850,1850,1850,1850,1850,1920,1850,1850];
ende    = [2100,2100,2100,2100,2100,2100,2099,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100];
start   = start(models_select);
ende    = ende(models_select);
Tanom_smiles   = NaN(length(models),max(ensmem),refende-refstart);
Traw_smiles    = Tanom_smiles;
Praw_smiles    = Tanom_smiles;
Praw_anom_smiles = Tanom_smiles;
Praw_anom_rel_smiles = Tanom_smiles;
for m = 1:length(models)
  % ---
  vari  = 'tas';
  if strcmp(scen_smiles{m},'rcp85')==1 || strcmp(models{m},'CESM2')==1 || strcmp(models{m},'E3SMv2')==1
    in    = dir([pathin vari '_Amon_' models{m} '_*ens*_g025_UCRB_ts_anom.nc']);
  else
    in    = dir([pathin '/cmip6-ng/' vari '_Amon_' models{m} '_' scen_smiles{m} '_*ens*_g025_UCRB_ts_anom.nc']);
  end
  tmp0	= squeeze(ncread([in.folder '/' in.name],vari));
  for e = 1:ensmem(m)
    tmp = wy_mean(tmp0(e,:));
    tmp = tmp(refstart-start(m)+1:refende-start(m)); % cut to common length
    Traw_smiles(m,e,:) = tmp;
    Tanom_smiles(m,e,:)= tmp-nanmean(tmp(1:time_Tobs(end)-refstart)); % reference to 1950-1979
  end
  % ---
  vari  = 'pr';
  tf    = 86400*30.4; % --> mm/month (like obs)
  if strcmp(models{m},'CESM2')==1 || strcmp(models{m},'E3SMv2')==1
    tf  = 86400*30.4 * 1e3; % --> mm/month (like obs)
  end
  if strcmp(scen_smiles{m},'rcp85')==1 || strcmp(models{m},'CESM2')==1 || strcmp(models{m},'E3SMv2')==1
    in    = dir([pathin vari '_Amon_' models{m} '_*ens*_g025_UCRB_ts.nc']);
  else
    in    = dir([pathin '/cmip6-ng/' vari '_Amon_' models{m} '_' scen_smiles{m} '_*ens*_g025_UCRB_ts.nc']);
  end
  tmp0	= squeeze(ncread([in.folder '/' in.name],vari)) * tf;
  for e = 1:ensmem(m)
    tmp                   = wy(tmp0(e,:));
    tmp                   = tmp(refstart-start(m)+1:refende-start(m)); % cut to common length
    Praw_smiles(m,e,:)           = tmp;
    Praw_anom_smiles(m,e,:)      = tmp-nanmean(tmp(1:time_Pobs(end)-refstart)); % absolute anomalies to 1950-2023 mean
    Praw_anom_rel_smiles(m,e,:)  = (tmp/nanmean(tmp(1:time_Pobs(end)-refstart)))*100-100; % relative anomalies to 1950-2023 mean
  end
end



% == CMIP6 data =================================================================
% -- 35 models
% models_cmip6 = {'ACCESS-CM2','ACCESS-ESM1-5','AWI-CM-1-1-MR','BCC-CSM2-MR','CAMS-CSM1-0','CESM2','CESM2-WACCM','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1','CanESM5','CanESM5-CanOE','EC-Earth3','EC-Earth3-Veg','FGOALS-f3-L','FGOALS-g3','FIO-ESM-2-0','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G','HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MCM-UA-1-0','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NESM3','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};
% -- 29 models
models_cmip6 = {'ACCESS-CM2','ACCESS-ESM1-5','AWI-CM-1-1-MR','BCC-CSM2-MR','CAMS-CSM1-0',...
                'CESM2','CESM2-WACCM','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1',...
                'CanESM5','CanESM5-CanOE','EC-Earth3','EC-Earth3-Veg','FGOALS-f3-L',...
                'FGOALS-g3','FIO-ESM-2-0','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G',...
                'HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G',...
                'MCM-UA-1-0','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR',...
                'MRI-ESM2-0','NESM3','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};
models_cmip6(ismember(models_cmip6,'GFDL-CM4')) = []; % what's wrong with this one?
models_cmip6(ismember(models_cmip6,'NorESM2-LM')) = []; % pr: ssp585 r1i1p1f1 exists, but the historical r1i1p1f1 only starts in 1950 (?)
models_cmip6(ismember(models_cmip6,'FIO-ESM-2-0')) = []; % no ssp370
models_cmip6(ismember(models_cmip6,'HadGEM3-GC31-LL')) = []; % no ssp370
models_cmip6(ismember(models_cmip6,'NESM3')) = []; % no ssp370
models_cmip6(ismember(models_cmip6,'KACE-1-0-G')) = []; % piControl only 150 years

% scen    = {'ssp119','ssp126','ssp245','ssp370','ssp585'};
scen    = {'ssp126','ssp245','ssp370','ssp585'};
if length(scen)==5
  % -- the following models have ssp119:
  models_cmip6 = {'CAMS-CSM1-0','CNRM-ESM2-1','CanESM5','EC-Earth3-Veg-LR','EC-Earth3-Veg','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','GISS-E2-1-H','IPSL-CM6A-LR','MIROC-ES2L','MIROC6','MPI-ESM1-2-LR','MRI-ESM2-0','UKESM1-0-LL'};
end
models_select = 1:length(models_cmip6)

ensmem  = zeros(1,length(models_cmip6)) + 50; % limit all ensemble sizes to 50 (CanESM5)
start   = 1850;
ende    = 2100;

Tanom         = NaN(length(models_cmip6),length(scen),max(ensmem),refende-refstart);
Traw          = Tanom;
Praw          = Tanom;
Praw_anom     = Tanom;
Praw_anom_rel = Tanom;

var1  = 'tas';
var2  = 'pr';

ensmem0 = NaN([length(models_cmip6),length(scen)]);
tic
for m = 1:length(models_cmip6)
  [num2str(m) '/' num2str(length(models_cmip6))]
  for s = 1:length(scen)
    int       = dir([pathin '/cmip6-ng/' var1 '_Amon_' models_cmip6{m} '_' scen{s} '_r*_g025_UCRB_ts_anom.nc']);
    inp       = dir([pathin '/cmip6-ng/' var2 '_Amon_' models_cmip6{m} '_' scen{s} '_r*_g025_UCRB_ts.nc']);
    int_e = [];
    for e = 1:length(int)
      int_e     = [int_e {int(e).name(end-28:end-21)}];
    end
    inp_e = [];
    for e = 1:length(inp)
      inp_e     = [inp_e {inp(e).name(end-23:end-16)}];
    end
    ensmem_tmp    = length(int);
    ensmem0(m,s)  = ensmem_tmp;
    ee = 0;
    for e = 1:ensmem_tmp
      idx = find(matches(inp_e, int_e{e}));
      if ~isempty(idx) == 1
        ee = ee+1;
        tmp0	= squeeze(ncread([int(e).folder '/' int(e).name], var1));
        tmp = wy_mean(tmp0);
        tmp = tmp(refstart-start+1:refende-start); % cut to common length
        Traw(m,s,ee,:) = tmp;
        Tanom(m,s,ee,:)= tmp-nanmean(tmp(1:time_Tobs(end)-refstart)); % reference to 1950-1979
        % ---
        vari  = 'pr';
        tf    = 86400*30.4; % --> mm/month (like obs)
        if strcmp(models_cmip6{m},'CESM2')==1 || strcmp(models_cmip6{m},'E3SMv2')==1
          tf  = 86400*30.4 * 1e3; % --> mm/month (like obs)
        end
        tmp0	= squeeze(ncread([inp(idx).folder '/' inp(idx).name], var2));
        tmp                      = wy(tmp0);
        tmp                      = tmp(refstart-start+1:refende-start); % cut to common length
        Praw(m,s,ee,:)           = tmp;
        Praw_anom(m,s,ee,:)      = tmp-nanmean(tmp(1:time_Pobs(end)-refstart)); % absolute anomalies to 1950-2023 mean
        Praw_anom_rel(m,s,ee,:)  = (tmp/nanmean(tmp(1:time_Pobs(end)-refstart)))*100-100; % relative anomalies to 1950-2023 mean
      end
    end
  end
end
toc
return




%% ================================================================================
% -- plotting
close all

figure1 = figure;
set(figure1, 'units', 'centimeters', 'pos', [50 10 20 12]) % on home screen

subplot(3,2,[1 3])
hold on
% -- CMIP6 (one member per model)
tmp1 = reshape(Tanom(:,:,1,:),[length(models_cmip6)*length(scen),149]);
plot(time,polyval(polyfit(time,prctile(tmp1,[5]),4),time),'b')
plot(time,polyval(polyfit(time,prctile(tmp1,[95]),4),time),'b')
h1 = plot(time,polyval(polyfit(time,nanmean(tmp1),4),time),'b','LineWidth',2)
% -- SMILEs (all members per model)
tmp2 = reshape(Tanom_smiles,[length(models)*20,149])
plot(time,polyval(polyfit(time,prctile(tmp2,[5]),4),time),'r')
plot(time,polyval(polyfit(time,prctile(tmp2,[95]),4),time),'r')
h2 = plot(time,polyval(polyfit(time,nanmean(tmp2),4),time),'r','LineWidth',2)
% -- SMILEs (one member per model)
tmp3 = reshape(Tanom_smiles(:,1,:),[length(models),149])
plot(time,polyval(polyfit(time,prctile(tmp3,[5]),4),time),'r--')
plot(time,polyval(polyfit(time,prctile(tmp3,[95]),4),time),'r--')
h3 = plot(time,polyval(polyfit(time,nanmean(tmp3),4),time),'r--','LineWidth',2)

box on
legend([h1 h2 h3],['CMIP6 (' num2str(length(models_cmip6)) ' x 1 x 4 = ' num2str(length(models_cmip6)*4) ')'],['SMILEs (' num2str(length(models)) ' x 20 x 1 = ' num2str(length(models)*20) ')'],['SMILEs (' num2str(length(models)) ' x 1 x 1 = ' num2str(length(models)) ')'],'location','northwest')
legend boxoff
title('(a) Temperature projections')
ylabel('\Delta T (\circC)')

subplot(3,2,5)
hold on
tmp1p1 = polyval(polyfit(time,prctile(tmp1,[5]),4),time);
tmp1p2 = polyval(polyfit(time,prctile(tmp1,[95]),4),time);
tmp2p1 = polyval(polyfit(time,prctile(tmp2,[5]),4),time);
tmp2p2 = polyval(polyfit(time,prctile(tmp2,[95]),4),time);
tmp3p1 = polyval(polyfit(time,prctile(tmp3,[5]),4),time);
tmp3p2 = polyval(polyfit(time,prctile(tmp3,[95]),4),time);
h1 = plot(time,(tmp2p2-tmp2p1)./(tmp1p2-tmp1p1)*100,'r')
h2 = plot(time,(tmp3p2-tmp3p1)./(tmp1p2-tmp1p1)*100,'r--')
hline(100,'k:')
box on
set(gca,'YLim',[75 115])
legend([h1 h2],['SMILEs (' num2str(length(models)) ' x 20 x 1 = ' num2str(length(models)*20) ')'],['SMILEs (' num2str(length(models)) ' x 1 x 1 = ' num2str(length(models)) ')'],'location','southwest')
legend boxoff
title('(c) Range of CMIP6 covered by SMILEs')
ylabel(['Range (%)'])
xlabel('Time (Year)')


subplot(3,2,[2 4])
hold on
tmp1 = reshape(Praw_anom_rel(:,:,1,:),[length(models_cmip6)*length(scen),149]);
plot(time,polyval(polyfit(time,prctile(tmp1,[5]),4),time),'b')
plot(time,polyval(polyfit(time,prctile(tmp1,[95]),4),time),'b')
h1 = plot(time,polyval(polyfit(time,nanmean(tmp1),4),time),'b','LineWidth',2)
tmp2 = reshape(Praw_anom_rel_smiles,[length(models)*20,149])
plot(time,polyval(polyfit(time,prctile(tmp2,[5]),4),time),'r')
plot(time,polyval(polyfit(time,prctile(tmp2,[95]),4),time),'r')
h2 = plot(time,polyval(polyfit(time,nanmean(tmp2),4),time),'r','LineWidth',2)
% -- SMILEs (one member per model)
tmp3 = reshape(Praw_anom_rel_smiles(:,1,:),[length(models),149])
plot(time,polyval(polyfit(time,prctile(tmp3,[5]),4),time),'r--')
plot(time,polyval(polyfit(time,prctile(tmp3,[95]),4),time),'r--')
h3 = plot(time,polyval(polyfit(time,nanmean(tmp3),4),time),'r--','LineWidth',2)
box on
title('(b) Precipitation projections')
ylabel('\Delta P (%)')

subplot(3,2,6)
hold on
tmp1p1 = polyval(polyfit(time,prctile(tmp1,[5]),4),time);
tmp1p2 = polyval(polyfit(time,prctile(tmp1,[95]),4),time);
tmp2p1 = polyval(polyfit(time,prctile(tmp2,[5]),4),time);
tmp2p2 = polyval(polyfit(time,prctile(tmp2,[95]),4),time);
tmp3p1 = polyval(polyfit(time,prctile(tmp3,[5]),4),time);
tmp3p2 = polyval(polyfit(time,prctile(tmp3,[95]),4),time);
plot(time,(tmp2p2-tmp2p1)./(tmp1p2-tmp1p1)*100,'r')
plot(time,(tmp3p2-tmp3p1)./(tmp1p2-tmp1p1)*100,'r--')
hline(100,'k:')
box on
set(gca,'YLim',[75 115])
title('(d) Range of CMIP6 covered by SMILEs')
ylabel(['Range (%)'])
xlabel('Time (Year)')

% tightfig

return
set(gcf,'PaperPositionMode','auto');
set(gcf,'renderer','Painters')
fileo = ['~/Dropbox/publication/lehner22_smile_perspective/fig/reservoir_model_cmip6_vs_smiles'];
print('-r300','-loose', '-depsc', ['' fileo '.eps'])
save2pdf(['' fileo '.pdf'])
return
