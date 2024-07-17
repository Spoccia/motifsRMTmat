function [Time4Clustering,timeforSubclustering] = cluster_RMT_Start(saveFeaturesPath,justSubCluster, subclusterflag, StrategyClustering, kindOfClustring, distanceUsed, TS_name, USER_OT_targhet, USER_OD_targhet )
    % StrategyClustering= strategy(strategyIDentifier)%2;%1;%3;%
    % 1 - create cluster of feature for the very same  varaites then  in each cluster do  adaptive kmeans on descriptors
    % 2 - create cluster of feature  on similar variates using Adaptive Kmeans then  for each cluster use adaptive kmeans on descriptors
    % 3 - old approach do clustering  then subclustering
    Time4Clustering = 0;
    timeforSubclustering = 0;

    %% read the  features
    savepath1 = [saveFeaturesPath,'feature_',TS_name,'.mat'];
    savepath2 = [saveFeaturesPath,'idm_',TS_name,'.mat'];
    savepath3 = [saveFeaturesPath,'MetaData_',TS_name,'.mat'];

    load(savepath1);
    load(savepath2);
    load(savepath3);
    indexfeatureGroup = (frame1(6,:)==USER_OT_targhet & frame1(5,:)==USER_OD_targhet);
    X=frame1(:,indexfeatureGroup);
    DepdScopeVector=csvread(strcat(saveFeaturesPath,'Distances',distanceUsed,'\DepdScopeVector_IM_',TS_name,'_DepO_',num2str(USER_OD_targhet),'_TimeO_',num2str(USER_OT_targhet),'.csv'));
    tic
    if (StrategyClustering == 1 | StrategyClustering == 4 | StrategyClustering == 7 | StrategyClustering == 10) %% we are interested into  same dependency scope
        possibleset= unique(X(1,:));
        AlltheCluster=[];
        Allthefeatures=[];
        Centroids =[];
        allclusterid=0;
        for classidlabel= 1:size(possibleset,2) % for each set of varaites  create a cluster
            idactfeatures= frame1(1,:)==possibleset(classidlabel);
            ActFeatures = X(:,idactfeatures);
            if(strcmp(kindOfClustring,'AKmeans')==1)
                if(size(ActFeatures,2)<=2)
                    C=ones(size(ActFeatures,2),1)+allclusterid;
                    mu=ActFeatures(11:end,1)';
                else
                    if StrategyClustering==10
                        [C,mu,inertia,tryK,startK]= adaptiveKmeansNormDistances(ActFeatures,2,0.02,2,'sqeuclidean');
                    else
                        [C,mu,inertia,tryK,startK]= adaptiveKmeans(ActFeatures,2,0.02,2,'sqeuclidean');
                    end
                    %                                             [C,mu,inertia,tryK,startK]= adaptiveKmeans(ActFeatures,2,0.02,2,'sqeuclidean');
                end
            elseif(strcmp(kindOfClustring,'DBScan')==1) % strategy==4
                if(size(ActFeatures,2)<=2)
                    C=ones(size(ActFeatures,2),1)+allclusterid;
                    mu=ActFeatures(11:end,1)';
                else
                    C=[];
                    varType=[];
                    if StrategyClustering == 4
                        [C, varType] = dbscan(ActFeatures(11:end,:)',2,'euclidean',0.5);
                    elseif StrategyClustering ==7
                        [C, varType] = dbscan(ActFeatures(11:end,:)',2,'euclidean');
                    end
                    C=C+1;
                    labels = unique(C);
                    mu=zeros(size(labels,1),128);
                    for clusterlabels=1:size(labels,1);
                        mu(clusterlabels,:)= mean(ActFeatures(11:end,C==labels(clusterlabels))');
                    end
                end
            end
            Allthefeatures=[Allthefeatures,ActFeatures];
            Centroids=[Centroids;mu];
            AlltheCluster=[AlltheCluster;C+allclusterid];
            allclusterid=max(AlltheCluster);
        end
        C=AlltheCluster;
        mu=Centroids;
        X=Allthefeatures;
    elseif(StrategyClustering == 2 | StrategyClustering == 5 | StrategyClustering == 8) %% we are interested  into croup of feature on similar variates then we apply  a clustering to get  this groups
        % we first use the depdscopevector to cluster  the features
        % with similar depepndency scope  then we use the
        % descriptor to cluster on the base of the time property
        [depd_Cluster,mu,inertia,tryK,startK]= adaptiveKmeansDependency(DepdScopeVector,2,0.02,1,'hamming');
        possibleset= unique(depd_Cluster);
        AlltheCluster=[];
        Allthefeatures=[];
        Centroids =[];
        allclusterid=0;
        for classidlabel= 1:size(possibleset,1) % for each set of varaites  create a cluster on descriptors
            idactfeatures= depd_Cluster == possibleset(classidlabel);
            ActFeatures = X(:,idactfeatures);
            if(strcmp(kindOfClustring,'AKmeans')==1)
                if(size(ActFeatures,2)<=2)
                    C=ones(size(ActFeatures,2),1)+allclusterid;
                    mu=ActFeatures(11:end,1)';
                else
                    [C,mu,inertia,tryK,startK]= adaptiveKmeans(ActFeatures,2,0.02,1,'sqeuclidean');
                end
            elseif(strcmp(kindOfClustring,'DBScan')==1) % strategy==5
                if(size(ActFeatures,2)<=2)
                    C=ones(size(ActFeatures,2),1)+allclusterid;
                    mu=ActFeatures(11:end,1)';
                else
                    C=[];
                    varType=[];
                    if StrategyClustering == 5
                        [C, varType] = dbscan(ActFeatures(11:end,:)',2,'euclidean',0.5);
                    elseif StrategyClustering == 8
                        [C, varType] = dbscan(ActFeatures(11:end,:)',2,'euclidean');
                    end
                    C=C+1;
                    labels = unique(C);
                    mu=zeros(size(labels,1),128);
                    for clusterlabels=1:size(labels,1);
                        mu(clusterlabels,:)= mean(ActFeatures(11:end,C==labels(clusterlabels))');
                    end
                end
            end
            Allthefeatures=[Allthefeatures,ActFeatures];
            Centroids=[Centroids;mu];
            AlltheCluster=[AlltheCluster;C+allclusterid];
            allclusterid=max(AlltheCluster);
        end
        C=AlltheCluster;
        mu=Centroids;
        X=Allthefeatures;

    elseif((StrategyClustering == 3 | StrategyClustering == 6 | StrategyClustering == 9| StrategyClustering ==11|StrategyClustering ==20) & justSubCluster==0)% classic strategy  we cluster all the features
        if(strcmp(kindOfClustring,'AKmeans')==1)
            if StrategyClustering ==20
                Method   = 'UPGMA';
                Metric   = @(X,Y)pdist([X;Y],'euclidean');%@(X,Y)abs(X-Y);
                Limit    = 0.2;
                Colormap = 'Cool';
                %[Figure,
                [Tree, Clusters, Roots] = hierarchicalcluster(X(11:end,:)',Method,Metric,'Limit',Limit,Colormap);
                allclusterid=0;
                Centroids=[];
                AlltheCluster=[];
                IDXfeatures= [];
                for clusi=1:size(Clusters,1)
                    FeaturesIDX=Clusters{clusi,1};

                    if size(Clusters{clusi,1},2)>1
                        C=ones(size(Clusters{clusi,1},2),1)+allclusterid;
                        mu = mean(X(11:end,FeaturesIDX)')';

                    else
                        C=ones(size(Clusters{clusi,1},2),1)+allclusterid;
                        mu= X(11:end,FeaturesIDX);
                    end
                    IDXfeatures=[IDXfeatures,Clusters{clusi,1}];
                    Centroids=[Centroids,mu];
                    AlltheCluster=[AlltheCluster,C'];
                    allclusterid=allclusterid+1;
                end
                C=AlltheCluster(IDXfeatures)';
                mu=Centroids';
            elseif StrategyClustering ==11
                [C,mu,inertia,tryK,startK]= adaptiveKmeansNormDistances(X,3,0.15,2,'sqeuclidean');
            else
                [C,mu,inertia,tryK,startK]= adaptiveKmeans(X,3,0.02,2,'sqeuclidean');%'cosine');%4th parameter will fix the step to 2 as default 0.02
            end
        elseif(strcmp(kindOfClustring,'DBScan')==1) % strategy==6
            C=[];
            varType=[];
            if StrategyClustering==6
                [C, varType] = dbscan(X(11:end,:)', 2,'euclidean',0.5);
            elseif StrategyClustering==9
                [C, varType] = dbscan(X(11:end,:)', 2,'euclidean');
            end
            C=C+1;
            labels = unique(C);
            mu=zeros(size(labels,1),128);
            for clusterlabels=1:size(labels,1);
                mu(clusterlabels,:)= mean(X(11:end,C==labels(clusterlabels))');
            end
        end
    end
    Time4Clustering=toc;
    % if (StrategyClustering ~= 3 | (StrategyClustering == 3 & justSubCluster==0))
    if(exist(strcat(saveFeaturesPath,'Distances',distanceUsed,'\ClusterStrategy_',num2str(StrategyClustering),'\'),'dir')==0)
        mkdir(strcat(saveFeaturesPath,'Distances',distanceUsed,'\ClusterStrategy_',num2str(StrategyClustering),'\'));
    end
    csvwrite(strcat(saveFeaturesPath,'Distances',distanceUsed,'\ClusterStrategy_',num2str(StrategyClustering),'\Cluster_IM_',TS_name,'_DepO_',num2str(USER_OD_targhet),'_TimeO_',num2str(USER_OT_targhet),'.csv'),C);
    csvwrite(strcat(saveFeaturesPath,'Distances',distanceUsed,'\ClusterStrategy_',num2str(StrategyClustering),'\Centroids_IM_',TS_name,'_DepO_',num2str(USER_OD_targhet),'_TimeO_',num2str(USER_OT_targhet),'.csv'),mu);
    csvwrite(strcat(saveFeaturesPath,'Distances',distanceUsed,'\ClusterStrategy_',num2str(StrategyClustering),'\Features_IM_',TS_name,'_DepO_',num2str(USER_OD_targhet),'_TimeO_',num2str(USER_OT_targhet),'.csv'),X);
    csvwrite(strcat(saveFeaturesPath,'Distances',distanceUsed,'\ClusterStrategy_',num2str(StrategyClustering),'\Time4Clustering',TS_name,'_DepO_',num2str(USER_OD_targhet),'_TimeO_',num2str(USER_OT_targhet),'.csv'),Time4Clustering);
    % end
    if subclusterflag == 1
        if (StrategyClustering == 3 | StrategyClustering == 6| StrategyClustering == 9 | StrategyClustering == 11 | StrategyClustering == 20)
            saveFeaturesPath=[datasetPath,subfolderPath,'Features_',FeaturesRM,'\',TS_name,'\'];
            depdOverLapThreshold = 1;
            timeforSubclustering = subCluster_Varaites(saveFeaturesPath,TS_name,num2str(StrategyClustering),distanceUsed,depdOverLapThreshold,USER_OT_targhet,USER_OD_targhet);
            csvwrite(strcat(saveFeaturesPath,'Distances',distanceUsed,'\ClusterStrategy_',num2str(StrategyClustering),'\timeforSubclustering',TS_name,'_DepO_',num2str(USER_OD_targhet),'_TimeO_',num2str(USER_OT_targhet),'.csv'),timeforSubclustering);
        end
    end

end

