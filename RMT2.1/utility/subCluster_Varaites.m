function  timeforSubclustering = subCluster_Varaites(saveFeaturesPath,TS_name,SizeofK,distanceUsed,depdOverLapThreshold,USER_OT_targhet,USER_OD_targhet)
%  % design as a row vector
%             reOrgFeatures       = [];
%             currentFeatureLabel = [];
%             reorgCentroids      = [];
%             currentFeatDepd     = [];
%             prevFeaturesLabel   = [];
savepath1 = [saveFeaturesPath, 'feature_', TS_name,'.mat'];
savepath2 = [saveFeaturesPath, 'idm_', TS_name,'.mat'];
savepath3 = [saveFeaturesPath, 'MetaData_', TS_name,'.mat'];
ClusterPath = [saveFeaturesPath,'Distances',distanceUsed,'\ClusterStrategy_',SizeofK,'\'];%,'\Cluster_',SizeofK,'\',distanceUsed,'\',typeofCluster];

load(savepath1);
load(savepath2);
load(savepath3);

timeforSubclustering=0;%zeros(1,4);
timeOctave = USER_OT_targhet;
depdOctave = USER_OD_targhet;

if(exist(strcat(ClusterPath, '\Cluster_IM_', TS_name, '_DepO_', num2str(depdOctave), '_TimeO_', num2str(timeOctave), '.csv'), 'file')~=0)
    
    % for features from each octave
    indexfeatureGroup = (frame1(6,:) == timeOctave & frame1(5,:) == depdOctave);
    X = frame1(:,indexfeatureGroup);
    depdScaleRead = csvread(strcat(saveFeaturesPath, 'Distances',distanceUsed,'\DepdScale_IM_', TS_name, '_DepO_', num2str(depdOctave),'_TimeO_',num2str(timeOctave),'.csv'));
    
    % read output from preivous k-means
    C = csvread(strcat(ClusterPath,'\Cluster_IM_',TS_name,'_DepO_',num2str(depdOctave),'_TimeO_',num2str(timeOctave),'.csv'));
    mu = csvread(strcat(ClusterPath,'\Centroids_IM_',TS_name,'_DepO_',num2str(depdOctave),'_TimeO_',num2str(timeOctave),'.csv'));
    clusterLabel = unique(C);
    numberOfClusters = length(clusterLabel);
    
    % subCluster label for interesting octave monotonically increases
    % for features from different octaves, subClusterLabel reset to 1
    subClusterLabel = 1;
    
    reOrgFeatures       = [];
    currentFeatureLabel = [];
    reorgCentroids      = [];
    currentFeatDepd     = [];
    prevFeaturesLabel   = [];
    
    tic;
    for clusterIndex = 1 : numberOfClusters
        % get the features from current clutser index
        clusterFeatures = X(:, C == clusterLabel(clusterIndex));
        
        % create sub-cluster below
        depdScale = depdScaleRead(:, C == clusterLabel(clusterIndex));
        
        % should use cell structure to index cluster features
        Depd_OverLapping = zeros(size(clusterFeatures, 2));
        
        for queryFeatureIndex = 1 : size(clusterFeatures, 2)
            % dataFeatureIndex starts from queryFeatureIndex since our goal is to group feature from the same simulation data
            for dataFeatureIndex = queryFeatureIndex : size(clusterFeatures, 2)
                % use symmetric approach to define the similiar motifs below
                %                         if(dataFeatureIndex==17)
                %                             'check'
                %                         end
                overlap = overlapJaccardSimilarity(depdScale(:, queryFeatureIndex), depdScale(:, dataFeatureIndex));
                Depd_OverLapping(queryFeatureIndex, dataFeatureIndex) = overlap;
                Depd_OverLapping(dataFeatureIndex, queryFeatureIndex) = overlap;
            end
        end
        Depd_OverLapping = Depd_OverLapping >= depdOverLapThreshold;
        stopFlag= 0;
        while(stopFlag == 0 )
            columnSum = sum(Depd_OverLapping); % row vector
            [maxSum, index] = max(columnSum);
            if(maxSum<=1)%(maxSum == 0) % check maxSum condition at the beginning Silv Did <=1 because 1 means single cluster
                stopFlag = 1;
            else
                subClusterFeatureIndex = Depd_OverLapping(:, index); % column vector
                %nonZeroEntryIndex = find(subClusterFeatureIndex ~= 0); % nonZeroEntryIndex is the actual feature index from current cluster
                subClusterFeatures =  clusterFeatures(:, subClusterFeatureIndex);
                
                % re-org features and the feathre labels
                reOrgFeatures       = [reOrgFeatures, subClusterFeatures];
                currentFeatureLabel = [currentFeatureLabel, repmat(subClusterLabel, 1, size(subClusterFeatures, 2))]; % design as a row vector
                reorgCentroids      = [reorgCentroids,mean(subClusterFeatures(11:end,:)')'];
                currentFeatDepd     = [currentFeatDepd,depdScale(:,subClusterFeatureIndex)];
                
                prevFeaturesLabel   = [prevFeaturesLabel,repmat(clusterIndex, 1, size(subClusterFeatures, 2))];
                % save as index: subClusterLabel
                subClusterLabel = subClusterLabel + 1;
                
                % now update Depd_OverLapping matrix
                %                         zeroColumn = zeros(size(subClusterFeatureIndex));
                Depd_OverLapping(:, subClusterFeatureIndex) = Depd_OverLapping(:, subClusterFeatureIndex)*0;%zeroColumn;
                Depd_OverLapping(subClusterFeatureIndex,:) = Depd_OverLapping(subClusterFeatureIndex,:)*0;%zeroColumn;
            end
        end
    end % end marking current cluster at interesting octave
    timeforSubclustering=toc;%(timeOctave*depdOctave)
    if(exist([ClusterPath,'\SplitVariate\'],'dir')==0)
        mkdir([ClusterPath,'\SplitVariate\']);
    end
    csvwrite([ClusterPath,'\SplitVariate\Features_IM_',TS_name,'_OT_',num2str(timeOctave),'_OD_',num2str(depdOctave),'.csv'],reOrgFeatures);
    csvwrite([ClusterPath,'\SplitVariate\Cluster_IM_',TS_name,'_OT_',num2str(timeOctave),'_OD_',num2str(depdOctave),'.csv'],currentFeatureLabel);
    csvwrite([ClusterPath,'\SplitVariate\Centroids_IM_',TS_name,'_OT_',num2str(timeOctave),'_OD_',num2str(depdOctave),'.csv'],reorgCentroids);
    csvwrite([ClusterPath,'\SplitVariate\DepdScale_IM_',TS_name,'_OT_',num2str(timeOctave),'_OD_',num2str(depdOctave),'.csv'],currentFeatDepd);
    csvwrite([ClusterPath,'\SplitVariate\ParentCluster_IM_',TS_name,'_OT_',num2str(timeOctave),'_OD_',num2str(depdOctave),'.csv'],prevFeaturesLabel);
    
end