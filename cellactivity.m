
FileList = dir('*.mat');
DataC    = cell(1, numel(FileList));
for iFile = 1:numel(FileList);
  FileData     = load(FileList(iFile).name);
  %DataC{iFile} = FileData.Plaque_dist2(:,1);
  DataC{iFile} = FileData.data;
  test = struct2cell(FileList);
  DataD{iFile} = test{1,iFile};
end
for i = 1:size (DataC,2);
    temp = DataC{i};
    temp = reshape(temp,[],1);% for cells
    temp = temp(all(~isnan(temp),2),:);
    temp0 = temp;
    temp0(temp0(:,1) < 0.0005, :) = [];
    mouseMean2= mean (temp);
    mouseMean3 = mean (temp0);
    hypo = 100*(sum(temp(:)<0.0001)/(size(temp,1)));
    hyperCAMK = 100*(sum(temp(:)>0.19)/(size(temp,1)));%7 per min %0.09 SST
    hyperPV = 100*(sum(temp(:)>0.33)/(size(temp,1)));%7 per min
    hyperSOM = 100*(sum(temp(:)>0.11)/(size(temp,1)));%7 per min
    %hyper = 100*(sum(temp(:)>0.066)/(size(temp,1)));%7 per min
    %save(DataD{i},'hypo','hyper','-append');
    hypo2(i,:) = hypo(:);
    hyperCAMK2(i,:) = hyperCAMK(:);
    hyperPV2 (i,:) = hyperPV(:);
    hyperSOM2 (i,:) = hyperSOM(:);
    mouseMean(i,:)= mouseMean2 (:);
    mouseMean0(i,:)= mouseMean3 (:);
end
save

