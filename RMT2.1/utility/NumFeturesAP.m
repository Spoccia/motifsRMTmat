path='D:\Motif_ResultsCleaning\Datasets\Image\data\Features\2_RME\';%T\';
filtering='2_FeaturesPruned Time X Dep Entropy\';%'';%'2_FeaturesPruned Time and Dep Entropy\';%
postpath='\Cluster_Fixed\Descriptor\afterPruning\ClusterMatlab\PrunedFeatures_IM_';
imname='1';
k=[5,10,15,20];
OD=2;
OT=2;
AllFeatures=[];
for i = 1:size(k,2);
    SizeofFeatures= [];
    for ii=1:OT
        for iii= 1:OD
            A=csvread([path,filtering,'DistancesDescriptors_C',num2str(k(i)),postpath,imname,'_DepO_',num2str(iii),'_DepT_',num2str(ii),'.csv'] );            
            SizeofFeatures= [SizeofFeatures;[k(i),ii,iii,size(A,2)]];            
        end
    end
    AllFeatures=[AllFeatures;SizeofFeatures];
  end
