function [corrpair] = corrfix_final3(Project_path, name)
% correlation analysis based on deconvolved events
load(sprintf('%s/analyzed/%s.mat', Project_path, name));
%% converting traces below a particular threshould to 0
thr = evalin('base','thr_inactive');
for j =  1:size(spikes,2)
    spikesnewCL2 = spikes(:, j);
    bdata2= dffilt(:, j);
    n_spikes2 = n_spikes(:, j);
    center3 = center(j,:);
    if  n_spikes2<= thr; 
        spikesnewCL2 (:,:)= 0;
        bdata2 (:,:)= 0;
        center3 (1,:)=0;
    else
        spikesnewCL2 = spikesnewCL2;
        center3 = center3;
        bdata2 = bdata2;
    end
    spikesnewCL(j, :) = spikesnewCL2(:)';
    center4(j, :) = center3(:)';
    bdata3(j, :) = bdata2(:)';
end
spikesnewCL= spikesnewCL';
bdata3 = bdata3';
n_spikes3 = n_spikes;
n_spikes3(n_spikes3<=thr) = [];
center5 = center4(any(center4,2),:); 
zeros_idx = find(all(spikesnewCL==0));
bdata3(:,zeros_idx) = [];
spikesnewCL(:,zeros_idx) = [];%deletes columns with all zero values

numtraces = size (spikesnewCL,2);
%% to deal with traces with all below threshould spiking

TF = isempty(spikesnewCL);% if all cells are zereo and the traces are empty output nan
if numtraces ==1 ; % convert noisy taces to 0 to eliminate them and avoid Oaisis touble
   corrpair= NaN;  corrpairAV=NaN; pvpair= NaN; corrP= NaN; corrpair2= NaN;  corrpairAV2=NaN; pvpair2= NaN; corrP2= NaN; n_spikes3=NaN;distances=NaN;  
   save(sprintf('%s/correlations/%s.mat', Project_path, name), 'name','corrP', 'corrpair', 'corrpairAV','pvpair', 'corrP2', 'corrpair2', 'corrpairAV2','pvpair2', 'n_spikes3', 'distances');
elseif TF ==1 
   corrpair= NaN;  corrpairAV=NaN; pvpair= NaN; corrP= NaN; corrpair2= NaN;  corrpairAV2=NaN; pvpair2= NaN; corrP2= NaN; n_spikes3=NaN;distances=NaN;  
   save(sprintf('%s/correlations/%s.mat', Project_path, name), 'name','corrP', 'corrpair', 'corrpairAV','pvpair', 'corrP2', 'corrpair2', 'corrpairAV2','pvpair2', 'n_spikes3', 'distances');
else
spikesnewCL =spikesnewCL;
distances = pdist (center5,'euclidean');
distances2 = squareform(distances);
distances = distances';

spikesnewCL=~~spikesnewCL; % binarize
spikesnewCL= double (spikesnewCL);



%%
dt =1/fs;

%% pearson no lag
[corrP, pvalue] = corr(spikesnewCL, 'type','Pearson');
corrP2 =corrP;
corrP2(eye(size(corrP2))==1) = nan;
corrpairAV= nanmean (corrP2);
corrP2 = triu(corrP,1); % removing repeated correlations 
ii=ones(size(corrP2));
idx= triu(ii,1);
corrP2(~idx)=nan; 
corrpair = reshape(corrP2',[],1); % changing the matrix toa acloumn
corrpair = corrpair(~isnan(corrpair))'; % changing the matrix to a acloumn
corrpair = corrpair'; %Pearson correlation witout lag
pvaluep =pvalue;
pvaluep(eye(size(pvaluep))==1) = nan;
%corrpairAV= nanmean (pvaluep)
pvaluep = triu(pvalue,1); % removing repeated correlations 
ii=ones(size(pvaluep));
idx= triu(ii,1);
pvaluep(~idx)=nan; 
pvpair = reshape(pvaluep',[],1); % changing the matrix toa acloumn
pvpair = pvpair(~isnan(pvpair)); % changing the matrix to a acloumn
pvpair = round (pvpair', 4)'; %Pearson correlation witout lag
%% binning 
binsize= 5;
r = rem(size(spikesnewCL,1),binsize);
if r==0;
spikesnewCL2 = squeeze(sum(reshape(spikesnewCL',numtraces,binsize,[]),2))';
else
extra =zeros(5-r, numtraces);
spikesnewCL = [spikesnewCL;extra];
spikesnewCL2 = squeeze(sum(reshape(spikesnewCL',numtraces,binsize,[]),2))';
end
    
    
spikesnewCL2=~~spikesnewCL2;
spikesnewCL2= double (spikesnewCL2);
%% %% pearson no lag binnned

[corrP3, pvalue3] = corr(spikesnewCL2, 'type','Pearson');
corrP3(eye(size(corrP3))==1) = nan;
corrpairAV2= nanmean (corrP3);
corrP3 = triu(corrP3,1); % removing repeated correlations 
ii=ones(size(corrP3));
idx= triu(ii,1);
corrP3(~idx)=nan; 
corrpair2 = reshape(corrP3',[],1); % changing the matrix toa acloumn
corrpair2 = corrpair2(~isnan(corrpair2))'; % changing the matrix to a acloumn
corrpair2 = corrpair2'; %Pearson correlation witout lag
pvaluep2 =pvalue3;
pvaluep2(eye(size(pvaluep2))==1) = nan;
%corrpairAV= nanmean (pvaluep)
pvaluep2 = triu(pvalue3,1); % removing repeated correlations 
ii=ones(size(pvaluep2));
idx= triu(ii,1);
pvaluep2(~idx)=nan; 
pvpair2 = reshape(pvaluep2',[],1); % changing the matrix toa acloumn
pvpair2 = pvpair2(~isnan(pvpair2)); % changing the matrix to a acloumn
pvpair2 = round (pvpair2', 4)'; %Pearson correlation witout lag
n = size (corrpair,1);
names    = cell(1, n);
names(:) = {name};
names = names';
save(sprintf('%s/correlations/%s.mat', Project_path, name), 'name','corrP', 'corrpair', 'corrpairAV','pvpair', 'corrP2', 'corrpair2', 'corrpairAV2','pvpair2', 'n_spikes3', 'distances');
end