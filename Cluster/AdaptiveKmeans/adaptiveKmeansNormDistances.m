function  [C1,mu1,inertia,tryK,startK]=adaptiveKmeansNormDistances(features1,K_start,saturation,Step,distance)
   isFound=false;
   inertia=[];
 %  inertiatest=[];
   AllFeaturesDistances = pdist2(features1(11:end,:)',features1(11:end,:)');
   maxdistance2instances= max(AllFeaturesDistances(:));
   if Step==0
       Step=2;
   end
   itr=1;
    startK=K_start;
    tryK=[];
    meanSilhouette=[];
    evaluation=[];
    C1=[];
    mu1=[];
    inp=[];
   while(~isFound)
      ProposedC1=C1;
      Proposedmu=mu1;
      try
%           if(size(features1(11:size(features1,1),:),2)-1<=(startK-1))
%               [C1,mu1,SUMD, D]=kmeans(features1(11:size(features1,1),:)',startK-2,'Distance',distance,'Replicates',5);%'Display','final'
%           else
            [C1,mu1,SUMD, D]=kmeans(features1(11:size(features1,1),:)',startK,'Distance',distance,'Replicates',5);%'Display','final'
%           end
      catch
            [C1,mu1,SUMD, D]=kmeans(features1(11:size(features1,1),:)',startK-2,'Distance',distance,'Replicates',5);%'Display','final'
       %     isFound=true;
      end
%       figure
%       [silh5,h] = silhouette(features1(11:end,:)',C1,'Euclidean');
%       meanSilhouette=[meanSilhouette,mean(silh5)];
%       evaluation = [evaluation,C1];
      
      % normalize the distances
      % get the number of instances in each cluster
      % sum the D rows
      % averageOf the SUM is SUMD
      % inertia should be teh average of the SUMD
      
      labels= unique(C1);
      D1 =zeros (size(D));
      D2=D1;
%       D3=D1;
     % maxdistance2instances1=max(max(D1));
      NumofIntancesforClusters=zeros(1,K_start);
      for i=1:length(labels)
          D1(C1==i,i)=D(C1==i,i);    % save just the distances between  instances and each centroid
          NumofIntancesforClusters(1,i)=sum(C1==i);
           D2(:,i)=D1(:,i)/maxdistance2instances;%/Maximum_D_C; %normalize the distances  with the max distance
            
      end
%       maxdistance2instances1=max(D1(:));
%       D3=D1/maxdistance2instances1;
%       NumofIntancesforClusters = sum(D1>0);
      
      SUMD1 = sum(D1);
      SUMD2 = sum(D2);
%       SUMD3 = sum(D3);
      MeanD1 = SUMD1./NumofIntancesforClusters;
      mean2 = sum(SUMD1)/sum(NumofIntancesforClusters);
%        MeanD2= mean(D2>0);
      MeanD2 = SUMD2./NumofIntancesforClusters;
      MD2=max(MeanD2);
      mean3  = sum(SUMD2)/sum(NumofIntancesforClusters);
   %   meanD3= sum(SUMD3)/sum(NumofIntancesforClusters);
       MeasureToUse=sum(SUMD2);%mean3;%mean2;%MD2;%sum(SUMD2)%sum(SUMD1);%;%mean(mean3);%MeanD2);% sum(SUMD2);%
%        MeasureToUse=sum(MeanD2);
%         saturation = /mean(NumofIntancesforClusters);

      % [C1,mu1,SUMD] = kmeans(features1(11:size(features1,1),:)',startK,'Distance',distance);%'sqeuclidean');%'cosine');%
       %% computeinertia
       
       if(itr==1)
%        inertia=[inertia,sum(SUMD)];
       inertia=[inertia,MeasureToUse];%mean(SUMD1)];
%       inertiatest=[inertiatest,meanD3];
%         inertia=[inertia,sum(SUMD1)/size(D,1)];
%        inp=[inp,mean3];

       tryK=[tryK,startK];
       startK=startK+Step;
       end
       if(itr>1)
%         inertia=[inertia,sum(SUMD)];
        inertia=[inertia,MeasureToUse];%mean(SUMD1)];
 %       inertiatest=[inertiatest,meanD3];
%         inp=[inp,mean3];
%         inertia=[inertia,sum(SUMD1)/size(D,1)];
        tryK=[tryK,startK];
        error = abs(inertia(itr) - inertia(itr-1));
        if error<=saturation  || startK >= (size(features1,2)-1) || inertia(itr)<=saturation %& inertia(itr)<=0.05
%             eva = evalclusters(features1(11:end,:)',evaluation,'CalinskiHarabasz');
%               C1=ProposedC1;
%               mu1=Proposedmu;
%         GetbadClusters =   MeanD2<1;
       
            isFound=true;
        end
        startK=startK+Step;        
       end
       if(~isFound)
       itr=itr+1;
       end
   end
end