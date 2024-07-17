

%% Paramenters
datasetPath= 'D:\Motif_ResultsCleaning\Datasets\Image\data\Features\';
Experiment='5_RME\';
TESTcase='5_FeaturesPruned Time X Dep Entropy\';
TestID='1';

afterpruning=0;
singleCluster=1;
clusterofInterest=1;
figinimage=1;
steppause=0;
sizeofCluster='5';
DO='2';
TO='2';

AllFeatureLocation= [datasetPath,Experiment,TESTcase,'feature_',TestID,'.mat'];
load(AllFeatureLocation);
centroid= csvread(strcat(datasetPath,Experiment,TESTcase,'DistancesDescriptors_C',sizeofCluster,'\Cluster_Fixed\Descriptor\ClusterMatlab\Centroids_IM_',TestID,'_DepO_',DO,'_TimeO_',TO,'.csv'));
cluster= csvread(strcat(datasetPath,Experiment,TESTcase,'DistancesDescriptors_C',sizeofCluster,'\Cluster_Fixed\Descriptor\ClusterMatlab\Cluster_IM_',TestID,'_DepO_',DO,'_TimeO_',TO,'.csv'));
dependency= csvread(strcat(datasetPath,Experiment,TESTcase,'DistancesDescriptors_C',sizeofCluster,'\DepdScale_IM_',TestID,'_DepO_',DO,'_TimeO_',TO,'.csv'));


indexfeatureGroup = (frame1(6,:)==str2num( TO) & frame1(5,:)==str2num(DO));
X=frame1(:,indexfeatureGroup);
if (afterpruning==1)
    centroid= csvread('D:\Motif_ResultsCleaning\Datasets\Image\data\Features\5_RME\5_FeaturesPruned Time X Dep Entropy\DistancesDescriptors_C10\Cluster_Fixed\Descriptor\ClusterMatlab\Centroids_IM_1_DepO_2_TimeO_2.csv');
    cluster= csvread('D:\Motif_ResultsCleaning\Datasets\Image\data\Features\5_RME\5_FeaturesPruned Time X Dep Entropy\DistancesDescriptors_C10\Cluster_Fixed\Descriptor\ClusterMatlab\Cluster_IM_1_DepO_2_TimeO_2.csv');
    dependency= csvread('D:\Motif_ResultsCleaning\Datasets\Image\data\Features\5_RME\5_FeaturesPruned Time X Dep Entropy\DistancesDescriptors\DepdScale_IM_1_DepO_2_TimeO_2.csv');
    X=csvread('featuresfrom the pruning');
end

  clusterLabel = unique(cluster);
  nCluster     = length(clusterLabel);      
  if (singleCluster ==1)
      nCluster=clusterofInterest;
  end
for i=clusterofInterest : nCluster
    idxClusterfeatures= cluster==i;
    SpecificCentroid= centroid(i,:);
    ClusterFeatures = X(:,idxClusterfeatures);
    ClusterDependency= dependency(:,idxClusterfeatures);
    [~,tempidx]=sort(ClusterFeatures(2,:));
    ClusterFeatures = ClusterFeatures(:,tempidx);
    ClusterDependency= ClusterDependency(:,tempidx);
    DistBFeatures= pdist2(ClusterFeatures(11:138,:)',ClusterFeatures(11:138,:)','cosine')/2;
    for ii=1:size(ClusterFeatures,2)
        DistBFeatures(ii,ii)=inf;
        [minval,nearestFeature]= min(DistBFeatures(ii,:));
        DistBFeatures(ii,ii)=0;
        [maxval,farestFeature]= max(DistBFeatures(ii,:));
        avgFi = mean(DistBFeatures(ii,:));
        fig = zeros(size(data));
        if(figinimage==1)
            fig(ClusterDependency(ClusterDependency(:,ii)>0,ii), round(ClusterFeatures(2,ii)-3*ClusterFeatures(4,ii)):round(ClusterFeatures(2,ii)+3*ClusterFeatures(4,ii)))=...
                data(ClusterDependency(ClusterDependency(:,ii)>0,ii), round(ClusterFeatures(2,ii)-3*ClusterFeatures(4,ii)):round(ClusterFeatures(2,ii)+3*ClusterFeatures(4,ii)));
        else
            fig=data(ClusterDependency(ClusterDependency(:,ii)>0,ii), round(ClusterFeatures(2,ii)-3*ClusterFeatures(4,ii)):round(ClusterFeatures(2,ii)+3*ClusterFeatures(4,ii)));
        end
        figure
            imshow(uint8(fig));
            title({strcat('OT_',TO, 'OD_',DO, 'CL_',num2str(i),'ID',num2str(ii)),...%,'\n',...
                   strcat('NDF=ID',num2str(nearestFeature),'val=', num2str(minval)),...%,'\n',... 
                   strcat('MDF=ID',num2str(farestFeature),'val=', num2str(maxval)),...%,'\n',... 
                   strcat('AVG-Distance=', num2str(avgFi))});
        if(steppause==1)
            pause
        end
    end
end
  