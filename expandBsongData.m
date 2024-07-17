%%%% Expland Series

PathSeries = 'D:\Motif_Results\Datasets\BirdSong\data_SHORT52\';
PathDestination = 'D:\Motif_Results\Datasets\BirdSong\data\';

for i = 1:154
    load([PathSeries,num2str(i),'.mat']);%csvread([PathSeries,num2str(i),'.csv']);
    tempdata = data';
   data=[];
   for varid = 1:size(tempdata,2)
       t = tempdata(:,varid);
        data = [data; interp1(1:numel(t),t,1:.5:numel(t))];
   end
  save([PathDestination,num2str(i),'.mat'],'data');
end