function [ TimeforPruningClustering ] = KmeansPruning(saveFeaturesPath, TEST, imagepath,specificimagepath,imagename,typeofCluster,strategy,prunewith,distanceUsed ,FeaturesRM,USER_OT_targhet,USER_OD_targhet,saveMotifImages,PathOldFeatures )
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
TimeforPruningClustering=0;
Matlab=1;
%saveFeaturesPath=[imagepath,specificimagepath,'Features_',FeaturesRM,'\',TEST,'\'];
% if(isempty(PathOldFeatures)==false)
%     saveFeaturesPath=[PathOldFeatures,TEST,'\'];
% end
savepath1 = [saveFeaturesPath,'feature_',imagename,'.mat'];
savepath2 = [saveFeaturesPath,'idm_',imagename,'.mat'];
savepath3 = [saveFeaturesPath,'MetaData_',imagename,'.mat'];
ImageSavingPath=[saveFeaturesPath,'Distances',distanceUsed,'\ClusterStrategy_',strategy,'\AP\'];%,prunewith,'\imageMotifs\'];
PrunedClusterPath=[saveFeaturesPath,'Distances',distanceUsed,'\ClusterStrategy_',strategy,'\AP\',typeofCluster];
ClusterPath = [saveFeaturesPath,'Distances',distanceUsed,'\ClusterStrategy_',strategy,'\'];
load(savepath1);
load(savepath2);
load(savepath3);
TS =data;
clustindfix=0;

k = USER_OT_targhet;
j = USER_OD_targhet;

MotifBag=[];
prunedFeaturesCluster=[];
prunedCluster=[];
prunedDepScale=[];
%            clustindfix=clustindfix+1;
indexfeatureGroup = (frame1(6,:)==k & frame1(5,:)==j);
X=frame1(:,indexfeatureGroup);
dpscale = csvread(strcat(saveFeaturesPath,'Distances',distanceUsed,'\DepdScale_IM_',imagename,'_DepO_',num2str(j),'_TimeO_',num2str(k),'.csv'));
C = csvread(strcat(ClusterPath,'\Cluster_IM_',imagename,'_DepO_',num2str(j),'_TimeO_',num2str(k),'.csv'));%_',imagename,'_',num2str(p),'.csv'));
mu = csvread(strcat(ClusterPath,'\Centroids_IM_',imagename,'_DepO_',num2str(j),'_TimeO_',num2str(k),'.csv'));%_',imagename,'_',num2str(p),'.csv'));
if(Matlab==1)
    mu=mu';
end
clusterLabel = unique(C);
nCluster     = length(clusterLabel);

TimeforPruningClustering_0=[];
    for i=1:nCluster
        tic;
        Centroid_desriptor = mu(:, i);
        %% features and depscale of each feature in cluster i
        A = X(:, C == clusterLabel(i));
        if (size(A,2)>1)
            B =dpscale(:,C == clusterLabel(i));
            [featsize,numfeatures]= size(A);
            descr = A(11:featsize,:);

            single_std_cluster=std(descr')';
            single_avg_cluster = mean(descr')';

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
                    IdxClFeat(1,k1)= sum(distancecentroid(:,k1) <= single_std_cluster*3)==128;
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
            TimeforPruningClustering_0=[TimeforPruningClustering_0,toc];
            %                      tempMinScope1 = min(3*A1(4, :)); % temporal scope of the features in cluster i
            if (size(A1,2)>0)
                timescope= A1(4,:)*3;
            end
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
                    if(exist([ImageSavingPath,'octaveT_',num2str(k),'_octaveD_',num2str(j)],'dir')==0)
                        mkdir([ImageSavingPath,'octaveT_',num2str(k),'_octaveD_',num2str(j),'\']);
                    end
                    %                         figure1=figure;
                    figure1 = plot_RMTmotif_on_data(data, MotifBag{i}.startIdx, MotifBag{i}.depd,MotifBag{i}.Tscope);
                    filename=[ImageSavingPath,'octaveT_',num2str(k),'_octaveD_',num2str(j),'\TS_',imagename,'_octT_',num2str(k),'_octD_',num2str(j),'_M_',num2str(i),'.eps'];
                    saveas(figure1,filename,'epsc');
                end
                prunedFeaturesCluster=[prunedFeaturesCluster,A1];
                prunedDepScale = [prunedDepScale,B1];
                prunedsymbols = ones(1,size(A1,2))*i;
                prunedCluster=[prunedCluster,prunedsymbols];
            end
        end
    end
    TimeforPruningClustering = sum(TimeforPruningClustering_0(:));
    if(exist(PrunedClusterPath,'dir')==0)
        mkdir(PrunedClusterPath);
    end
    save(strcat(PrunedClusterPath,'\Motif_',imagename,'_DepO_',num2str(j),'_DepT_',num2str(k),'.mat'),'MotifBag');
    close all;
    csvwrite(strcat(PrunedClusterPath,'\PrunedCluster_IM_',imagename,'_DepO_',num2str(j),'_DepT_',num2str(k),'.csv'),prunedCluster);
    csvwrite(strcat(PrunedClusterPath,'\Centroids_IM_',imagename,'_DepO_',num2str(j),'_DepT_',num2str(k),'.csv'),mu);
    csvwrite(strcat(PrunedClusterPath,'\PrunedFeatures_IM_',imagename,'_DepO_',num2str(j),'_DepT_',num2str(k),'.csv'),prunedFeaturesCluster);
    csvwrite(strcat(PrunedClusterPath,'\PrunedDepScaleFeatures_IM_',imagename,'_DepO_',num2str(j),'_DepT_',num2str(k),'.csv'),prunedDepScale);

end

