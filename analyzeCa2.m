% for analysis of CamKII expressing gcamp6s, malgamal@mgh.harvard.edu
% no neuropil correction for interneurons


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
%
if exist('Plaque') ==1
    Plaque = Plaque;
else
    Plaque = NaN;
end
%}
%% for traces with neuropil Rois and S2P files
if exist('Fneu') ==1
    iscell2=iscell(:,1);
    Fb= transpose(mean(F')./mean(Fneu'));% brightness relative to neruopil
    Remove2=find(Fb<=1.04);% chacnge to 1.03 to remove dim cells
    Remove1=find(iscell2==0);
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
    data = data';
    data_original = data';
end
%% ROI centers for distance 
if exist('Fneu') ==1
    for i = 1:size (stat,2);
        dim= stat{1,i};
        dim2 = dim.med;
        dim3= [dim2, i];
        center(i, :) = dim3(:)';
    end
    center2 = center;% ROI centers and theri new order (columns 3)
else 
center2 = NaN;
end
%}
%% high pass filter
fs = 1/ dt;% original sample rate

%% neuropil correction and scaling for negative
if exist('Fneu') ==1
   data = F;
    for n = 1:size (data,1);
    ndata3 = data (n,:);
    mindata2 = min (ndata3);
    dataorginal= data_original (n,:);
        if mindata2 <0 % adding an offset for negative data
           ndata2 = ndata3+abs(mindata2);% to avoid negatives
        else
           ndata2 = ndata3;
        end
    ndata (n,:) = ndata2 (:)';
    end
else
    ndata= data;
end


t_time2 = (length (ndata))/fs;  

%% Moving averge baselining, adpted fron feng paper , (chen 2020) without high passfilter

% now take a moving average of a few frames for a first-pass low-pass
% filter
%framesToMovAvg = 6;

% Now let's estimate the baselines. The quantile cut-off is calculated
% based on the statistics of the F distributions.
baselineWindow = 250;
blCutOffs = computeQunatileCutoffs(ndata);
somaticF_BLs=slidingBaseline2(ndata,baselineWindow,blCutOffs);

% now the DF is the (F-BL)/BL
bdata = (ndata - somaticF_BLs)./somaticF_BLs;
bdata = bdata'*100;% %df/f
for n = 1:size (bdata,2);
    y2 = bdata(:,n);
    miny =abs (min(y2));
    y = y2 + miny;
    bmin= prctile (y2, 5);
    [b1, sn1] = estimate_baseline_noise(y, bmin);
    if b1 > prctile (y, 70)
       b3 = 0;
       Possible_Baseline_Error = {name, n}
    else
       b3 = b1-miny;
    end
    b (:,n)= b1(:)';
    b2(:,n)= b3(:)';
    sn (:,n)= sn1(:)';
    bmin2 (:,n)= bmin(:)';
end
%% deconvolution using Oasis
bdata = bdata-b2;
n_traces = size (bdata,2);
t_time = (length (bdata))/fs;
bdata = double(bdata);
for j =  1:n_traces;
    y = bdata(:, j);
    g = exp(-(1/(1.25*fs)));         % AR coefficient 
 
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
traces2 = traces; % to be used to remove zero elements
spikes2 = spikes; % to be used to remove zero elements
%% number of spikes and spike rate
n_traces2 = size(spikes, 2);
for j =  1:n_traces2
    spikes3 = spikes(:, j);
    Nspikes = nnz(spikes3);
    Nspikes= Nspikes/ t_time;
    n_spikes(j, :) = Nspikes(:)';
end
n_spikes = transpose(n_spikes); % spikes per min

%% ploting
toriginal= 0:dt:t_time2-dt;
%
figure;plot(toriginal,ndata,'DisplayName','original data');
title('original')
xlabel('Time (s)')
ylabel('au')
plottools('on')
figure;plot(toriginal,dffilt,'DisplayName',name);
title('baselined df')
xlabel('Time (s)')
ylabel('%\DeltaF/F_{0}')
plottools('on')
figure;plot(toriginal,traces,'DisplayName','denoised traces');
title('denoised traces')
xlabel('Time (s)')
ylabel('%\DeltaF/F_{0}')
plottools('on')
figure;plot(toriginal, spikes,'DisplayName',name);
title('events')
xlabel('Time (s)')
ylabel('%\DeltaF/F_{0}')
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
