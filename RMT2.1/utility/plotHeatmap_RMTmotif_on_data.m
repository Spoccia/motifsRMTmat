function figure1 = plotHeatmap_RMTmotif_on_data(data, motif_idx, motif_dim,motif_lenght)

% Normalize in 0-255
variate =30;
time= 3;
AllMin = min(data');
AllMax = max(data');
div=AllMax'-AllMin';
y=zeros(size(data));
data =data';
%% plot the data
for i = 1:size(data, 2)
    data(:, i) = data(:, i) - min(data(:, i));
    data(:, i) = data(:, i) / max(data(:, i));
    data(:, i) = data(:, i) + (i - 1) * 1.1;
   % plot(data(:, i), 'color', 'k');
end
data=data';
AllMin = min(data');
AllMax = max(data');
div=AllMax'-AllMin';
for i = 1: size(AllMin,2)
y(i,:)=((data(i,:)-AllMin(i)).*255)/(AllMax(i)-AllMin(i));
end
figure
imshow(y)
Mask = ones(variate,time);
heatmapnoborder=[];

for i=0:size(AllMin,2)-1
 heatmapnoborder = [heatmapnoborder;zeros(variate,size(data,2)*time)];
 for j=0:size(data,2)-1
     heatmapnoborder(i*variate+1:i*variate+variate,j*time+1:j*time+time)= Mask*y(i+1,j+1);
 end
end
figure
imshow(heatmapnoborder)
HeatmapRGB(:,:,1) = heatmapnoborder;
motifHeatmap= zeros(size(heatmapnoborder));
for i = 1:length(motif_idx)
    for k = 1:length(motif_dim{i})
        sub_len=motif_lenght{i}*time; 
        motif_location = motif_idx(i)*time :motif_idx(i)*time + sub_len - 1;
        motifHeatmap( motif_dim{i}(k)*variate : motif_dim{i}(k)*variate+variate,motif_location) = heatmapnoborder(motif_dim{i}(k)*variate : motif_dim{i}(k)*variate+variate,motif_location);
    end
end

HeatmapRGB(:,:,2) = motifHeatmap;
HeatmapRGB(:,:,3) = zeros(size(motifHeatmap))*255;
figure
imshow(HeatmapRGB);