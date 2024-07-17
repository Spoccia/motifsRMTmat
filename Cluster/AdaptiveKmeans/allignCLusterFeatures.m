% This matlab script reads a cluster of a features and saves features that have similiar variates as one cluster
% function returned= allignCLusterFeatures()
clc;
clear;

SubDSPath = 'data\';%'FlatTS_MultiFeatureDiffClusters\';%'CosineTS_MultiFeatureDiffClusters\';%'MultiFeatureDiffClusters\';
datasetPath = 'D:\Motif_Results\Datasets\Mocap\';
subfolderPath = '';%'Z_A_Temp_C\';%
FeaturesRM ='RMT';%'RME';%

% Path Parameters
TEST ='1';%
typeofCluster='ClusterMatlab';
distanceUsed ='Descriptor';%'Amplitude_Descriptor';%
SizeofK = 'Fixed';
TS_name = TEST;
imagename = TS_name;
saveFeaturesPath = [datasetPath, subfolderPath, 'Features_', FeaturesRM, '\', TS_name, '\'];

savepath1 = [saveFeaturesPath, 'feature_', TS_name,'.mat'];
savepath2 = [saveFeaturesPath, 'idm_', TS_name,'.mat'];
savepath3 = [saveFeaturesPath, 'MetaData_', TS_name,'.mat'];
ClusterPath = [saveFeaturesPath,'DistancesDescriptors\Cluster_',SizeofK,'\',distanceUsed,'\',typeofCluster];
%         ImageSavingPath=[saveFeaturesPath,'DistancesDescriptors\Cluster_',K_valuesCalc,'\',distanceUsed,'\BP_Kmeans_CosineDescriptor'];%\imageMotifs\',imagename];
%         RebSeriesPath = [saveFeaturesPath,'DistancesDescriptors\Cluster_',K_valuesCalc,'\',distanceUsed,'\BP_Kmeans_CosineDescriptor\rebClusters\'];
load(savepath1);
load(savepath2);
load(savepath3);
depdOverLapThreshold = 0.5;

            reOrgFeatures       = [];
            currentFeatureLabel = [];
            reorgCentroids      = [];
            currentFeatDepd     = [];
            prevFeaturesLabel   = [];

 % design as a row vector
for timeOctave = 2:DeOctTime
    for depdOctave = 2:DeOctDepd
        if(exist(strcat(ClusterPath, '\Cluster_IM_', imagename, '_DepO_', num2str(depdOctave), '_TimeO_', num2str(timeOctave), '.csv'), 'file')~=0)
            
            % for features from each octave
            indexfeatureGroup = (frame1(6,:) == timeOctave & frame1(5,:) == depdOctave);
            X = frame1(:,indexfeatureGroup);
            depdScaleRead = csvread(strcat(saveFeaturesPath, 'DistancesDescriptors\DepdScale_IM_', imagename, '_DepO_', num2str(depdOctave),'_TimeO_',num2str(timeOctave),'.csv'));
            
            % read output from preivous k-means
            C = csvread(strcat(ClusterPath,'\Cluster_IM_',imagename,'_DepO_',num2str(depdOctave),'_TimeO_',num2str(timeOctave),'.csv'));
            mu = csvread(strcat(ClusterPath,'\Centroids_IM_',imagename,'_DepO_',num2str(depdOctave),'_TimeO_',num2str(timeOctave),'.csv'));
            clusterLabel = unique(C);
            numberOfClusters = length(clusterLabel);
            
            % subCluster label for interesting octave monotonically increases
            % for features from different octaves, subClusterLabel reset to 1
            subClusterLabel = 1;             
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
                        overlap = overlapJaccardSimilarity(depdScale(:, queryFeatureIndex), depdScale(:, dataFeatureIndex));
                        Depd_OverLapping(queryFeatureIndex, dataFeatureIndex) = overlap;
                        Depd_OverLapping(dataFeatureIndex, queryFeatureIndex) = overlap;
                    end
                end
                Depd_OverLapping = Depd_OverLapping > depdOverLapThreshold;
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
                        zeroColumn = zeros(size(subClusterFeatureIndex));
                        %use matrix product to set the feature identified
                        %in a cluster to 0 
                        Depd_OverLapping(:, subClusterFeatureIndex) = Depd_OverLapping(:, subClusterFeatureIndex)*0;%zeroColumn;
                        Depd_OverLapping(subClusterFeatureIndex,:) = Depd_OverLapping(subClusterFeatureIndex,:)*0;%zeroColumn;
                    end
                end
            end % end marking current cluster at interesting octave
        end
    end
 end


% end

