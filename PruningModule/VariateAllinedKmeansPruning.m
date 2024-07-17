function [TimeforPruningSubClustering ] = VariateAllinedKmeansPruning(saveFeaturesPath, TEST, imagepath,specificimagepath,imagename,typeofCluster,strategy,prunewith,distanceUsed ,FeaturesRM,USER_OT_targhet,USER_OD_targhet,saveMotifImages,PathOldFeatures )
%(TEST, imagepath,specificimagepath,imagename,typeofCluster,K_valuesCalc,prunewith,distanceUsed ,DictionarySize,histdataimage,FeaturesRM,cleanfeatures,saveMotifImages )
%KMEANSPRUNING Summary of this function goes here
% imagepath=path of the original timeseries
% specificimagepath= ppath to go inside  where the dataset is located
% imagename= the timeseries name
% typeofCluster= ClusterKmedoids or ClusterMatlab it identify the subfolder where  the cluster should be saved and  then it can be readed
% K_valuesCalc =  the value of k can be 'Fixed';%'Threshould';% 'Computed'; a different name  will  help us to go in the specific cluster folder created for the test
% prunewith= Descriptor, Amplitude_Descriptor, Amplitude_Descriptor_overlapping
% distanceUsed= 'Descriptor';%'Amplitude_Descriptor';% it is used inside the cluster to define if use the amplitude or not in the cluster
% DictionarySize= vector containing the k used for each cluster.
% TEST= test folder path
%   Detailed explanation goes here
%% load features
%     if (strcmp(typeofCluster,'ClusterMatlab') ~= 1)
%     'wrong cluster !!!!'
%     typeofCluster
%     pause;
%         return;
%     end
% try
Matlab=1;
TimeforPruningSubClustering=0;

%saveFeaturesPath=[imagepath,specificimagepath,'Features_',FeaturesRM,'\',TEST,'\'];
% if(isempty(PathOldFeatures)==false)
%     saveFeaturesPath=[PathOldFeatures,TEST,'\'];
% end
savepath1 = [saveFeaturesPath,'feature_',imagename,'.mat'];
savepath2 = [saveFeaturesPath,'idm_',imagename,'.mat'];
savepath3 = [saveFeaturesPath,'MetaData_',imagename,'.mat'];
ClusterPath = [saveFeaturesPath,'Distances',distanceUsed,'\ClusterStrategy_',strategy,'\SplitVariate'];
ImageSavingPath=[ClusterPath,'\AP_VA'];
PrunedClusterPath = [ClusterPath,'\AP_VA\',typeofCluster];
load(savepath1);
load(savepath2);
load(savepath3);
TS =data;
clustindfix=0;
k = USER_OT_targhet;
j = USER_OD_targhet;

MotifBag=[];
clustindfix=clustindfix+1;
X=csvread(strcat(ClusterPath,'\Features_IM_',imagename,'_OT_',num2str(k),'_OD_',num2str(j),'.csv'));%frame1(:,indexfeatureGroup);
prunedFeaturesCluster=[];
prunedCluster=[];
prunedDepScale=[];

dpscale = csvread(strcat(ClusterPath,'\DepdScale_IM_',imagename,'_OT_',num2str(k),'_OD_',num2str(j),'.csv'));

C = csvread(strcat(ClusterPath,'\Cluster_IM_',imagename,'_OT_',num2str(k),'_OD_',num2str(j),'.csv'));
mu = csvread(strcat(ClusterPath,'\Centroids_IM_',imagename,'_OT_',num2str(k),'_OD_',num2str(j),'.csv'));
if(Matlab==0)
    mu=mu';
end



clusterLabel = unique(C);
nCluster     = length(clusterLabel);
%                  %% for each cluster

TimeforPruningSubClustering_0=[];
SubclusterDescriptors= [];
for i=1:nCluster
    tic;
    %% features and depscale of each feature in cluster i
    A = X(:, C == clusterLabel(i));
    B =dpscale(:,C == clusterLabel(i));
    [featsize,numfeatures]= size(A);
    descr = A(11:featsize,:);
    
    
    single_std_cluster= std(descr')';
    single_avg_cluster = mean(descr')';
    SubclusterDescriptors= [SubclusterDescriptors;single_avg_cluster];
    CentroidCalc= repmat(single_avg_cluster,1,numfeatures);%Centroid_desriptor
    
    std_cluster = std2(descr');
    avg_cluster = mean2(descr');
    
    %%calculate  descriptor distance between  the  centroid descriptor
    distancecentroid =abs(CentroidCalc-descr);% pdist2(Centroid_desriptor',descr');%,KmeansDescmetric);
    %                     for k1=1:numfeatures
    %                         centroid_distDescriptors(1,k1)= pdist([act_centroid';descr(:,k1)'],KmeansDescmetric);
    %                     end
    IdxClFeat=zeros(1,numfeatures);
    if(strcmp(prunewith,'Descriptor')==1 || strcmp(prunewith,'Amplitude_Descriptor')==1 ||strcmp(prunewith,'Amplitude_Descriptor_overlapping')==1)
        %                        'Prune using just Descriptors'
        for k1=1:numfeatures
            IdxClFeat(1,k1)= sum(distancecentroid(:,k1) <= single_std_cluster*3)==size(distancecentroid,1);
        end
        A1= A(:,IdxClFeat==1);
        B1= B(:,IdxClFeat==1);
        A=A1;
        B=B1;
        %FinalScore= descr_Score(distancecentroid,numfeatures);
    end
    if(strcmp(prunewith,'Amplitude_Descriptor')==1 || strcmp(prunewith,'Amplitude_Descriptor_overlapping')==1)
        'Prune using  Descriptors + Amplitude'
        X_Amp1=amplitudediff(TS,A,gss1,idm1);
        [~,numfeaturesAMP]= size(A);
        AVG_AMP=sum(X_Amp1)/size(X_Amp1,1);
        %% compute amplitude similarity of the features
        AMP_Score= zeros(1,numfeaturesAMP);
        for iii=1:numfeaturesAMP
            Dist2 = abs(X_Amp1(iii)-AVG_AMP)/(X_Amp1(iii)+AVG_AMP);
            AMP_Score(1,iii) = 1/(1+Dist2);
        end
        A1=A(:,AMP_Score > 0.9);
        B1=B(:,AMP_Score > 0.9);
        A=A1;
        B=B1;
        %FinalScore= Amp_descr_Score(distancecentroid,numfeatures,X_Amp1,AVG_AMP); % 1-combinedscore
    end
    if(strcmp(prunewith,'Amplitude_Descriptor_overlapping')==1)
        'Prune using  Descriptors + Amplitude + overlappingvariates'
        %% use overlapping to prune  the feature selected
        %B = depd scale of the cluster
        %[A1,B1]=pruneoverlap(A,B);
        
    end
    
    A1=A;
    B1=B;
    
    
    if (size(A1,2)>0)
        timescope= A1(4,:)*3;
    end
    TimeforPruningSubClustering_0=[TimeforPruningSubClustering_0,toc];
    if(size(A1,2)>1)
        MotifBag{i}.features=A1;
        StartID= round(A1(2,:)-timescope);
        StartID(StartID <1)=1;
        MotifBag{i}.startIdx = StartID';%round(A1(2,:)-timescope)';
        %                         MotifBag{i}.depd=B1;
        %                         MotifBag{i}.Tscope= 2* timescope(:);
        
        for iterator=1:size(MotifBag{i}.startIdx,1)
            MotifBag{i}.depd{iterator}=B1(B1(:,iterator)>0,iterator);
            intervaltime=(round((A1(2,iterator)-timescope(iterator))) : (round((A1(2,iterator)+timescope(iterator)))));
            MotifBag{i}.Tscope{iterator}= size(intervaltime(intervaltime>0 & intervaltime<=size(data,2)),2);%2* timescope(:);
        end
        if(saveMotifImages==1)
            if(exist([ImageSavingPath,'\octaveT_',num2str(k),'_octaveD_',num2str(j)],'dir')==0)
                mkdir([ImageSavingPath,'\octaveT_',num2str(k),'_octaveD_',num2str(j),'\']);
            end
            %                         figure1=figure;
            figure1 = plot_RMTmotif_on_data(data, MotifBag{i}.startIdx, MotifBag{i}.depd,MotifBag{i}.Tscope);
            filename=[ImageSavingPath,'\octaveT_',num2str(k),'_octaveD_',num2str(j),'\TS_',imagename,'_octT_',num2str(k),'_octD_',num2str(j),'_M_',num2str(i),'.eps'];%'.jpg'];
            saveas(figure1,filename,'epsc');
        end
        prunedFeaturesCluster=[prunedFeaturesCluster,A1];
        prunedDepScale = [prunedDepScale,B1];
        prunedsymbols = ones(1,size(A1,2))*i;
        prunedCluster=[prunedCluster,prunedsymbols];
    end
    
    
end
TimeforPruningSubClustering=sum(TimeforPruningSubClustering_0);
if(exist(PrunedClusterPath,'dir')==0)
    mkdir(PrunedClusterPath);
end
save(strcat(PrunedClusterPath,'\Motif_',imagename,'_DepO_',num2str(j),'_DepT_',num2str(k),'.mat'),'MotifBag');
close all;
csvwrite(strcat(PrunedClusterPath,'\PrunedCluster_IM_',imagename,'_DepO_',num2str(j),'_DepT_',num2str(k),'.csv'),prunedCluster);
csvwrite(strcat(PrunedClusterPath,'\Centroids_IM_',imagename,'_DepO_',num2str(j),'_DepT_',num2str(k),'.csv'),mu);
csvwrite(strcat(PrunedClusterPath,'\Centroids_subCluster_IM_',imagename,'_DepO_',num2str(j),'_DepT_',num2str(k),'.csv'),SubclusterDescriptors);
csvwrite(strcat(PrunedClusterPath,'\PrunedFeatures_IM_',imagename,'_DepO_',num2str(j),'_DepT_',num2str(k),'.csv'),prunedFeaturesCluster);
csvwrite(strcat(PrunedClusterPath,'\PrunedDepScaleFeatures_IM_',imagename,'_DepO_',num2str(j),'_DepT_',num2str(k),'.csv'),prunedDepScale);
% catch
% end
end

