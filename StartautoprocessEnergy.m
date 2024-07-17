close all;
clc;
clear;

%% Location Configuration
PATH_dataset ='E:\RMT\Dataset\Energy\data\';
saveFeaturesPath='E:\RMT\Dataset\Energy\MotifsFeatures\';
subfolderPath= '';%'Z_A_Temp_C\';%
FeaturesRM ='RMT';

TS_num = 709

%% Action to perform
FeatureExtractionFlag = 1; % 1 do it others  skip
createDependencyScale = 1; % 1 do it others  skip
Cluster=1 ; % 1 do it others  skip
    subclusterflag = 0;%1 do it others  skip
    justSubCluster=0; % in the case of strategy 3  we can do just  subclusteringt
pruneCluster = 1;
savecaracteristics = 1;
%% paramters
USER_OT_targhet=2;
USER_OD_targhet=2;
distanceUsed='Descriptor';% use just descriptors to  cluster
strategy=[1,3,10,11,4,20,6,7,9]; % use strategy 1,3,10,11
%                 if StrategyClustering >3  & StrategyClustering < 10
%                     kindOfClustring= 'DBScan';%
kindOfClustring= 'AKmeans';

prunewith='Descriptor';% use this strategy to prune  the outbound features ina  cluster
%% printing functionality
    saveMotifBP = 0; % show the clusters before pruning
    saveMotifAP = 0; % show the clusters after  pruning


for TSname =1:TS_num
    
    TS_name = num2str(TSname);
    
    disp(['EXECUTING JOB ON TIMESERIES: ', PATH_dataset,TS_name,'.mat']);
    %% TIMER
    TIMEFOROCTAVE = 0;
    TimeComputationDepdScale = 0;
    Time4Clustering =0;
    TimeforPruningClustering =0;
    TimeforPruningSubClustering=0;
    timeforSubclustering=0;

    
    
    if FeatureExtractionFlag == 1
        disp('---- EXECUTING FeatureExtractioN...');
        %TIMEFOROCTAVE = Energy_FE_start(PATH_dataset,saveFeaturesPath,TS_name, USER_OT_targhet, USER_OD_targhet);
        TIMEFOROCTAVE = Energy_FE_old_start(PATH_dataset,saveFeaturesPath,TS_name, USER_OT_targhet, USER_OD_targhet);
    end
    
    %% create dependency
    if(createDependencyScale==1)
        %% read the features
        savepath1 = [saveFeaturesPath,'feature_',TS_name,'.mat'];
        savepath2 = [saveFeaturesPath,'idm_',TS_name,'.mat'];
        savepath3 = [saveFeaturesPath,'MetaData_',TS_name,'.mat'];
        saveCSVDepd= strcat(saveFeaturesPath,'Distances',distanceUsed,'\DepdScale_IM_',TS_name,'_DepO_',num2str(USER_OD_targhet),'_TimeO_',num2str(USER_OT_targhet),'.csv');
        savevectorDepd = strcat(saveFeaturesPath,'Distances',distanceUsed,'\DepdScopeVector_IM_',TS_name,'_DepO_',num2str(USER_OD_targhet),'_TimeO_',num2str(USER_OT_targhet),'.csv');
        disp('---- EXECUTING createDependencyScale...');
        TimeComputationDepdScale = Crete_saveDepdScale(savepath1,savepath2,savepath3,USER_OT_targhet,USER_OD_targhet,saveCSVDepd,savevectorDepd,strcat(saveFeaturesPath,'Distances',distanceUsed,'\'));
    end
    
    %% to search for pattern ove different variates groups use subclusterflag = 0;
    for strategyIDentifier = 1:1 %2:2:4
        %% to search over same variates:
        %             for strategyIDentifier = 1:4 % 2 and 4 are slow if
        %             subclaster is active
        %                 subclusterflag = 1;
        StrategyClustering= strategy(strategyIDentifier);%2;%1;%3;%
        % 1 - create cluster of feature for the very same  varaites then  in each cluster do  adaptive kmeans on descriptors
        % 2 - create cluster of feature  on similar variates using Adaptive Kmeans then  for each cluster use adaptive kmeans on descriptors
        % 3 - old approach do clustering  then subclustering
        kindOfClustring= 'AKmeans';
        if StrategyClustering >3  & StrategyClustering < 10
            kindOfClustring= 'DBScan';%
        end
        
        % reset time
        Time4Clustering =0;
        TimeforPruningClustering =0;
        TimeforPruningSubClustering=0;
        timeforSubclustering=0;
        
        %% DO CLUSTERING
        if (Cluster==1 | justSubCluster==1)
            
            
            [Time4Clustering,timeforSubclustering] = cluster_RMT_Start(saveFeaturesPath,justSubCluster, subclusterflag, StrategyClustering, kindOfClustring, distanceUsed, TS_name, USER_OT_targhet, USER_OD_targhet);
        end
        
        %% DO PRUNING        
        if(pruneCluster==1)
            if (StrategyClustering == 3 | StrategyClustering == 6| StrategyClustering == 9  | StrategyClustering ==11 |StrategyClustering ==20 )
                TimeforPruningClustering = KmeansPruning(TS_name,PATH_dataset,subfolderPath,TS_name,kindOfClustring,num2str(StrategyClustering),prunewith,distanceUsed ,FeaturesRM,USER_OT_targhet,USER_OD_targhet,saveMotifAP);%1);
                if subclusterflag == 1
                    TimeforPruningSubClustering = VariateAllinedKmeansPruning(saveFeaturesPath,TS_name,PATH_dataset,subfolderPath,TS_name,kindOfClustring,num2str(StrategyClustering),prunewith,distanceUsed ,FeaturesRM,USER_OT_targhet,USER_OD_targhet,saveMotifAP);
                end
            else
                TimeforPruningClustering = KmeansPruning(saveFeaturesPath,TS_name,PATH_dataset,subfolderPath,TS_name,kindOfClustring,num2str(StrategyClustering),prunewith,distanceUsed ,FeaturesRM,USER_OT_targhet,USER_OD_targhet,saveMotifAP);
            end
        end
        %% DO SAVE RESULTS                
        if(savecaracteristics==1)
            savepath1 = [saveFeaturesPath,'feature_',TS_name,'.mat'];
            savepath2 = [saveFeaturesPath,'idm_',TS_name,'.mat'];
            savepath3 = [saveFeaturesPath,'MetaData_',TS_name,'.mat'];
            load(savepath1);
            load(savepath2);
            load(savepath3);
            a=[];
            k=USER_OT_targhet;
            j=USER_OD_targhet;
            indexfeatureGroup = (frame1(6,:)==k & frame1(5,:)==j);
            X=frame1(:,indexfeatureGroup);
            SizeFeaturesforImages=[k,j,size(X,2)];
            csvwrite(strcat(saveFeaturesPath,'NumFeatures.csv'),SizeFeaturesforImages);
            %                         xlswrite(strcat(saveFeaturesPath,'NumFeatures.xls'),SizeFeaturesforImages);
            
            %                         col_header={char(strcat('OT',num2str(USER_OT_targhet),'_OD',num2str(USER_OD_targhet)))};%{'OT1_OD1','OT1_OD2','OT2_OD1','OT2_OD2'};
            %                         rowHeader ={'FeatureEstraction';'ComputationDepdScale';'Clustering';'VaraiteAllineament';'PruningStandarDev_V_allined';'PruningStandarDev_Clusters'};
            csvwrite(strcat(saveFeaturesPath,'Strategy_',num2str(StrategyClustering),'_TIME1.csv'),[TIMEFOROCTAVE;TimeComputationDepdScale;Time4Clustering;timeforSubclustering;TimeforPruningSubClustering;TimeforPruningClustering]);
            %xlswrite(strcat(saveFeaturesPath,'Strategy_',num2str(StrategyClustering),'_TIME1.xls'),[TIMEFOROCTAVE;TimeComputationDepdScale;Time4Clustering;timeforSubclustering;TimeforPruningSubClustering;TimeforPruningClustering],'TIME','B2');
            %                         xlswrite(strcat(saveFeaturesPath,'Strategy_',num2str(StrategyClustering),'_TIME1.xls'),rowHeader,'TIME','A2');
            %                         xlswrite(strcat(saveFeaturesPath,'Strategy_',num2str(StrategyClustering),'_TIME1.xls'),col_header,'TIME','B1');
        end
    end
end


