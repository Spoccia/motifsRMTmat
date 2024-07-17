function [ output_args ] = ShowVaraiteallineadCluster(TEST,imagepath,specificimagepath,imagename,K_valuesCalc,distanceUsed,typeofCluster,histdataimage,FeaturesRM,saveMotifImages )
%SHOWKMEANSCLUSTER Summary of this function goes here
%   Detailed explanation goes here
% if (strcmp(typeofCluster,'ClusterMatlab')~=1)
%     'wrong cluster !!!!'
%     typeofCluster
%     pause;
%     return;
% end
saveFeaturesPath=[imagepath,specificimagepath,'Features_',FeaturesRM,'\',TEST,'\'];
savepath1 = [saveFeaturesPath,'feature_',imagename,'.mat'];
savepath2 = [saveFeaturesPath,'idm_',imagename,'.mat'];
savepath3 = [saveFeaturesPath,'MetaData_',imagename,'.mat'];

ClusterPath = [saveFeaturesPath,'Distances',distanceUsed,'\Cluster_',K_valuesCalc,'\SplitVariate'];%'\Cluster_',K_valuesCalc,'\',distanceUsed,'\',typeofCluster,'\SplitVariate'];
ImageSavingPath=[ClusterPath,'\BP_VA'];%\imageMotifs\',imagename];
RebSeriesPath = [ClusterPath,'\BP_VA\info\'];

load(savepath1);
load(savepath2);
load(savepath3);
colours= ['b';'g';'r';'c';'m';'y';'k';'w'];
symbols= ['.';'o';'x';'+';'*';'s';'d';'v'];
symbComb=[];
for i =1: length(symbols)
    for j = 1:length(colours)
        symbComb = [symbComb;strcat(colours(j,1),symbols(i,1))];
    end
end

for k=1:DeOctTime
    for j=1:DeOctDepd
        if(exist(strcat(ClusterPath,'\Cluster_IM_',imagename,'_OT_',num2str(k),'_OD_',num2str(j),'.csv'),'file')~=0)
            % indexfeatureGroup = (frame1(6,:)==k & frame1(5,:)==j);
            X=csvread(strcat(ClusterPath,'\Features_IM_',imagename,'_OT_',num2str(k),'_OD_',num2str(j),'.csv'));%frame1(:,indexfeatureGroup);
            %             DictionarySizeApplied= floor(abs(size(X,2))/10);
            %              if(abs(size(X,2))>=DictionarySizeApplied)
            dpscale = csvread(strcat(ClusterPath,'\DepdScale_IM_',imagename,'_OT_',num2str(k),'_OD_',num2str(j),'.csv'));
            %                 csvread(strcat(saveFeaturesPath,'DistancesDescriptors\DepdScale_IM_',imagename,'_DepO_',num2str(j),'_TimeO_',num2str(k),'.csv'));
            
            C = csvread(strcat(ClusterPath,'\Cluster_IM_',imagename,'_OT_',num2str(k),'_OD_',num2str(j),'.csv'));
            mu = csvread(strcat(ClusterPath,'\Centroids_IM_',imagename,'_OT_',num2str(k),'_OD_',num2str(j),'.csv'));
            clusterLabel = unique(C);
            nCluster     = length(clusterLabel);
%             dataid=zeros(size(data,1),size(data,2),nCluster);
%             histdataid=zeros(size(data,1),size(data,2),nCluster);
%             FeatureLocation=zeros(size(data,1),size(data,2),3,nCluster);
            MotifBag=[];
            for ii=1:nCluster
                A = X(:, C == clusterLabel(ii));
                B =dpscale(:,C == clusterLabel(ii));
                timescope= A(4,:)*3;
                
                %% use overlapping to prune  the feature selected
                %                     [A1,B1]=pruneoverlap(A,B);
                %                     A=A1;
                %                     B=B1;
                if (size(A,2)>0)
                    timescope= A(4,:)*3;
                end
                MotifBag{ii}.features=A;
                StartID= round(A(2,:)-timescope);
                StartID(StartID <1)=1;
                MotifBag{ii}.startIdx = StartID';%round(A1(2,:)-timescope)';
                for iterator=1:size(MotifBag{ii}.startIdx,1)
                    MotifBag{ii}.depd{iterator}=B(B(:,iterator)>0,iterator);
                    intervaltime=(round((A(2,iterator)-timescope(iterator))) : (round((A(2,iterator)+timescope(iterator)))));
                    MotifBag{ii}.Tscope{iterator}= size(intervaltime(intervaltime>0 & intervaltime<=size(data,2)),2);%2* timescope(:);
                end
                if(saveMotifImages==1)
                    if(exist([ImageSavingPath,'\octaveT_',num2str(k),'_octaveD_',num2str(j),],'dir')==0)
                        mkdir ([ImageSavingPath,'\octaveT_',num2str(k),'_octaveD_',num2str(j),'\']);
                    end
                    figure1 = plot_RMTmotif_on_data(data, MotifBag{ii}.startIdx, MotifBag{ii}.depd,MotifBag{ii}.Tscope);
                    filename=[ImageSavingPath,'\octaveT_',num2str(k),'_octaveD_',num2str(j),'\HistIm_',imagename,'_octT_',num2str(k),'_octD_',num2str(j),'_M_',num2str(ii),'.eps'];%'.jpg'];
                    saveas(figure1,filename,'epsc');
                end
                %                     %% close  this portion of code
                %                     % try to add square under the feature
                %                     Starting = zeros(size(data,1),size(data,2));
                %                     Ending   = zeros(size(data,1),size(data,2));
                %
                %                     for iii=1: size(A,2)
                %                         intervaltime=(round((A(2,iii)-timescope(iii))) : (round((A(2,iii)+timescope(iii)))));
                %                         dataid(B((B(:,iii)>0),iii),intervaltime((intervaltime>0 & intervaltime<=size(data,2))),ii)= data(B(B(:,iii)>0,iii),intervaltime((intervaltime>0 & intervaltime<=size(data,2))));
                %                         histdataid(B((B(:,iii)>0),iii),intervaltime((intervaltime>0 & intervaltime<=size(data,2))),ii)= histdataimage(B(B(:,iii)>0,iii),intervaltime((intervaltime>0 & intervaltime<=size(data,2))));
                %
                %                         Xs=min(intervaltime((intervaltime>0 & intervaltime<=size(data,2))));
                %                         Xe=max(intervaltime((intervaltime>0 & intervaltime<=size(data,2))));
                %
                %                         FeatureLocation(B((B(:,iii)>0),iii),Xs,1,ii)= Starting(B((B(:,iii)>0),iii),Xs)+1;
                %                         FeatureLocation(B((B(:,iii)>0),iii),Xe,2,ii)= Ending(B((B(:,iii)>0),iii),Xe)+1;
                %                         FeatureLocation(B((B(:,iii)>0),iii),intervaltime((intervaltime>0 & intervaltime<=size(data,2))),3,ii)= data(B(B(:,iii)>0,iii),intervaltime((intervaltime>0 & intervaltime<=size(data,2))));
                %
                %                         Starting(B((B(:,iii)>0),iii),Xs) = Starting(B((B(:,iii)>0),iii),Xs)+1;
                %                         Starting(min(B((B(:,iii)>0),iii)),(Xs):(Xe)) = Starting(min(B((B(:,iii)>0),iii)),(Xs):(Xe))+1;
                %                         Ending(B((B(:,iii)>0),iii),Xe) = Ending(B((B(:,iii)>0),iii),Xe)+1;
                %                         Ending(max(B((B(:,iii)>0),iii)),(Xs):(Xe)) = Ending(max(B((B(:,iii)>0),iii)),(Xs):(Xe))+1;
                %                     end
                %
                %                     if(size(A,2)>1)
                %                         if(exist([ImageSavingPath,'\octaveT_',num2str(k),'_octaveD_',num2str(j),],'dir')==0)
                %                             mkdir ([ImageSavingPath,'\octaveT_',num2str(k),'_octaveD_',num2str(j),'\']);
                %                         end
                %
                %                         Starting= Starting*255;
                %                         Ending= Ending*255;
                %                         Combined(:,:,1)= uint8(Ending);%red
                %                         Combined(:,:,2)= uint8(Starting);%red
                %                         Combined(:,:,3)= (histdataid(:,:,ii));%blue
                %                         %                         figure
                %                         %                         imshow(Combined)
                %
                %                         imwrite(uint8(histdataid(:,:,ii)),[ImageSavingPath,'\octaveT_',num2str(k),'_octaveD_',num2str(j),'\HistIm_',imagename,'_octT_',num2str(k),'_octD_',num2str(j),'_Cl_',num2str(ii),'.jpg']);
                %                         imwrite(uint8(dataid(:,:,ii)),[ImageSavingPath,'\octaveT_',num2str(k),'_octaveD_',num2str(j),'\Im_',imagename,'_octT_',num2str(k),'_octD_',num2str(j),'_Cl_',num2str(ii),'.jpg']);
                %                         imwrite(Combined,[ImageSavingPath,'\octaveT_',num2str(k),'_octaveD_',num2str(j),'\CHistIm_',imagename,'_octT_',num2str(k),'_octD_',num2str(j),'_Cl_',num2str(ii),'.jpg']);
                %                     end
            end
            if(exist(RebSeriesPath,'dir')==0)
                mkdir (RebSeriesPath);
            end
            save(strcat(RebSeriesPath,'MotifBag_',imagename,'_Toctave_',num2str(k),'_Doctave_',num2str(j),'_KC_',num2str(nCluster),'.mat'),'MotifBag');
            close all;
            %                 save(strcat(RebSeriesPath,'Series_Feature_',imagename,'_Toctave_',num2str(k),'_Doctave_',num2str(j),'_KC_',num2str(nCluster),'.mat'),'FeatureLocation');
            %                 save(strcat(RebSeriesPath,'RebSeries_',imagename,'_Toctave_',num2str(k),'_Doctave_',num2str(j),'_dic_',num2str(nCluster),'.mat'),'dataid');
            %             end
        end
    end
end
end

