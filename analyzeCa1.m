% for analysis of gcamp6s and 7s calcium datasets, malgamal@mgh.harvard.edu
% with neuropil correction for excitatory cells

function [n_spikes] = analyzeCa(Project_path, name)
load(sprintf('%s/data/%s.mat', Project_path, name));
oasis_setup;
smin = evalin('base','smin');
lambda = evalin('base','lambda');
thr = evalin('base','thr_inactive');
if exist('dt') ==1
    dt=dt;
else
    dt = 1/ops.fs;
end
%% preprocessing and removing dim cells for traces with neuropil ROIs and S2P files
if dt/0.2 <=0.9
   dt =dt;
elseif dt/0.2>1.1;
   dt=dt;
else dt =0.2;
end

dt =dt;
if exist('Fneu') ==1
    iscell2=iscell(:,1);
    Fb= transpose(mean(F')./mean(Fneu'));% brightness relative to neruopil
    Remove2=find(Fb<0.8);% only cells brighter than neuropil are included
    %Remove2=find(Fb>=0.95);% only cells brighter than neuropil are included
    Remove1=find(iscell2==0);% Remove non-cell ROIs
    Remove = sortrows([Remove2;Remove1]);
    Remove = unique(Remove);
    F(Remove,:) = [];
    Fneu(Remove,:) = [];
    %redcell(Remove,:) = [];
    for u= [Remove]; % removing iscell 0 from stat
    stat(u) = [];
    end
    data_original = F;
else
    data = data'; % if not using S2P
    data_original = data'; % if not using S2P
end
%% ROI centers for distance estimation, only if usinf

fs = 1/ dt;% original sample rate

%% neuropil correction and eliminating negative values
if exist('Fneu') ==1
   data = F -(0.7*Fneu);
   %data = F -(0.0*Fneu);% for SST and PV
    for n = 1:size (data,1);
    ndata3 = data (n,:);
    mindata2 = min (ndata3);
        if mindata2 <0 % adding an offset for negative data
           ndata2 = max(ndata3,0);% to convert negatives to 0
           warn = {'negative F-Fneu', name, n};
           disp (warn)
        else
           ndata2 = ndata3;
        end
    ndata (n,:) = ndata2 (:)';
    dim= stat{1,n};
    dim2 = dim.med;
    dim3= [dim2, n];
    center(n, :) = dim3(:)';% ROI centers and their S2p #
    end
else
    ndata= data;
    center = NaN;
end

  

%% Moving averge baselining, adpted fron feng paper , (chen 2020) without high passfilter

% now take a moving average of a few frames for a first-pass low-pass
% filter
%framesToMovAvg = 3;
%ndata = nPointMean(ndata,framesToMovAvg);

% Now let's estimate the baselines. The quantile cut-off is calculated
% based on the statistics of the F distributions.
baselineWindow = 250;
blCutOffs = computeQunatileCutoffs(ndata);
somaticF_BLs=slidingBaseline2(ndata,baselineWindow,blCutOffs);

% now the DF is the (F-BL)/BL
bdata = (ndata - somaticF_BLs)./abs(somaticF_BLs);
bdata = bdata'*100;% %df/f

for n = 1:size (bdata,2); % one more baseline correction using Oasis
    y2 = bdata(:,n);
    miny =abs (min(y2));
    %miny=1000;
    y = y2 + miny;
    bmin= prctile (y2, 5);
    [b1, sn1] = estimate_baseline_noise2(y, bmin);
    if b1 > prctile (y, 65)
       b2 = 0;
       warn = {'Possible Baseline Error', name, n};
       disp (warn)
    elseif b1 < prctile (y, 5)
       b2 = 0;
       warn = {'Possible Baseline Error', name, n};
       disp (warn)
    else
       b2 = b1-miny;
    end
    miny2(:,n)=miny(:)';
    btemp(:,n)= b1(:)';
    b(:,n)= b2(:)';
    sn (:,n)= sn1(:)';
    bmin2 (:,n)= bmin(:)';
    ytemp(:,n)= y(:)';
    if max (y2)> 800;
    warn = {'high df', name, n};
    disp (warn)
    elseif min (y2)< -100;
    warn = {'high df', name, n};
    disp (warn)
    else
    end
end
%% deconvolution using Oasis
bdata = bdata-b;

n_traces = size (bdata,2);
t_time = (length (bdata))/fs;
bdata = double(bdata);
g = exp(-(1/(1.25*fs))); % AR coefficient 
for j =  1:n_traces;
    y = bdata(:, j);
    [c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', g, 'foopsi', 'lambda', lambda, 'smin', smin);
    
    %%
    spikes(j, :) = s_oasis(:)';
    options2(j, :) = options(:)';
    traces(j, :) = c_oasis(:)';
    dffilt (j, :) = y (:)';   
end
%% df averaging
dffilt = transpose(dffilt);
spikes = transpose(spikes);
traces = transpose(traces);
spikes2 = spikes; % to be used to remove zero elements
%% number of spikes and spike rate
n_traces2 = size(spikes, 2);
for j =  1:n_traces2
    spikes3 = spikes(:, j);
    Nspikes = nnz(spikes3);
    Nspikes= Nspikes/ t_time;
    n_spikes(j, :) = Nspikes(:)';
end
n_spikes = double(n_spikes'); % spikes per min
toriginal= 0:dt:t_time-dt;
%% ploting
%
figure;plot(toriginal,data,'DisplayName','original data');
t= title('original');
t.FontSize = 48;
xlabel('Time (s)','FontSize', 32)
ylabel('au','FontSize', 32)
plottools('on')
figure;plot(toriginal,dffilt,'DisplayName',name);
t = title('Baselined');
t.FontSize = 48;
xlabel('Time (s)','FontSize', 32)
ylabel('%\DeltaF/F_{0}','FontSize', 32)
plottools('on')
figure;plot(toriginal,traces,'DisplayName','denoised traces');
t= title('Denoised traces');
t.FontSize = 48;
xlabel('Time (s)','FontSize', 32)
ylabel('%\DeltaF/F_{0}','FontSize', 32)
plottools('on')
figure;plot(toriginal, spikes,'DisplayName',name);
t= title('Events');
t.FontSize = 48;
xlabel('Time (s)', 'FontSize', 32)
ylabel('%\DeltaF/F_{0}','FontSize', 32)
plottools('on')
figurepalette('toff')
propertyeditor('off')
%close all;
%}
%% averaging  per fov

spikes2(spikes2<=0) = NaN; % % removing zeors and negatives before averaging
mspikes = mean (spikes2,'omitnan')'; % average spikes after removing silent cells
spike_rate0 = n_spikes; % icluding inactive cells
spike_rate = n_spikes';  % removingg inactive cell

%% writing tables
T = array2table (spike_rate);
table = sprintf('%s/analyzed/%s.mat', Project_path);
filename = [name, '_event_rate.xls'];
table = fullfile (table, filename);
writetable(T, table)
%% removing cells below thr
spike_rate(spike_rate<= thr) = NaN; % removing inactive cells
avspike_rate0 = mean(spike_rate0,'omitnan'); % icluding inactive cells
avspike_rate = mean(spike_rate,'omitnan'); % without inactive cells%% averaging 


%
save(sprintf('%s/analyzed/%s.mat', Project_path, name), 'name','avspike_rate0', 'avspike_rate','n_spikes','spikes', 'traces', 'dffilt', 'bdata', 'fs', 'center');
save(sprintf('%s/options/%s.mat', Project_path, name),  'name','spikes', 'traces', 'options2');
save(sprintf('%s/analyzeddf/%s.mat', Project_path, name),  'name', 'spikes','mspikes');
%}
