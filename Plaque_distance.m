%distance _estimation
function [avspike_rate] = Plaque_distance(Project_path, name)
load(sprintf('%s/data/%s.mat', Project_path, name));
load(sprintf('%s/analyzed/%s.mat', Project_path, name));
%Plaque = [100,356];
%% for traces with neuropil Rois bu no correction
%dt =0.2;
iscell2=iscell(:,1);
Fb= transpose(mean(F')./mean(Fneu'));% brightness relative to neruopil
Remove2=find(Fb<0.8);% only cells brighter than neuropil are included
Remove1=find(iscell2==0);% Remove non-cell ROIs
Remove = sortrows([Remove2;Remove1]);
Remove = unique(Remove);
F(Remove,:) = [];
Fneu(Remove,:) = [];;
for u= [Remove]; % removing iscell 0 from stat
    stat(u) = [];
end
%% sorting the vectors based on roi center (Y, then X)
for i = 1:size (stat,2);
    dim= stat{1,i};
    dim2 = dim.med;
    dim3= [dim2, i];
    center(i, :) = dim3(:)';
    if exist('Plaque') ==1
        pair =[dim2;Plaque];% 
        distance = pdist(pair,'euclidean'); % calculate the dsitance between cocordiantes
        distance2 = [distance,i]; % distance to plaque bewtween a given roi and plaque and roin numbe ins2p 
        Plaque_dist (i, :) = distance2(:)';
    else
        Plaque_dist2 = NaN;
end
    
end
%center2 = sortrows(center);% ROI centers and theri new order (columns 3)
center2 = center;
order= center2(:,3); % order of cells  in output files relative to S2P ROI#
F2=F(order,:);
Fneu2=Fneu(order,:);
if exist('Plaque') ==1
    pixtoum = 0.994; % 0.497for 512x512 0.994 for 256x256
    Plaque_dist(:, 1)=pixtoum.*Plaque_dist(:, 1);
    Plaque_dist2 = Plaque_dist(order,:);
else
    Plaque_dist2 = NaN;
end
data_original = F2;
save(sprintf('%s/plaque/%s.mat', Project_path, name), 'name', 'Plaque_dist2', 'center2', 'n_spikes');