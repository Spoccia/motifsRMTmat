function  [C1,mu1,inertia,tryK,startK]=adaptiveKmeansDependency(DepdScopeVector,K_start,saturation,Step,distance)%adaptiveKmeans(features1,K_start,saturation,Step,distance)
   isFound=false;
   inertia=[];
   AllFeaturesDistances = pdist2(DepdScopeVector',DepdScopeVector',distance);
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
   while(~isFound)
      ProposedC1=C1;
      Proposedmu=mu1;
      [C1,mu1,SUMD, D]=kmeans(DepdScopeVector',startK,'Distance',distance,'Replicates',5);%'Display','final'

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
      NumofIntancesforClusters=zeros(1,K_start);
      for i=1:length(labels)
          D1(C1==i,i)=D(C1==i,i);    % save just the instances of the distances from  each centroid
          NumofIntancesforClusters(1,i)=sum(C1==i);
           D2(:,i)=D1(:,i)/maxdistance2instances;%/Maximum_D_C; %normalize the distances  with the max distance
      end
      
%       NumofIntancesforClusters = sum(D1>0);
      
      SUMD1 = sum(D1);
      SUMD2 = sum(D2);
      MeanD1 = SUMD1./NumofIntancesforClusters;
      mean2 = SUMD1/sum(NumofIntancesforClusters);
%       MeanD2= mean(D2>0);
      MeanD2 = SUMD2./NumofIntancesforClusters;
      mean3  = sum(SUMD2)/sum(NumofIntancesforClusters);
       MeasureToUse=mean3;%sum(SUMD2)%sum(SUMD1);%;%mean(mean3);%MeanD2);% sum(SUMD2);%
%        MeasureToUse=sum(MeanD2);
%         saturation = /mean(NumofIntancesforClusters);

      % [C1,mu1,SUMD] = kmeans(features1(11:size(features1,1),:)',startK,'Distance',distance);%'sqeuclidean');%'cosine');%
       %% computeinertia
       if(itr==1)
%        inertia=[inertia,sum(SUMD)];
       inertia=[inertia,MeasureToUse];%mean(SUMD1)];
%         inertia=[inertia,sum(SUMD1)/size(D,1)];

       tryK=[tryK,startK];
       startK=startK+Step;
       end
       if(itr>1)
%         inertia=[inertia,sum(SUMD)];
        inertia=[inertia,MeasureToUse];%mean(SUMD1)];
%         inertia=[inertia,sum(SUMD1)/size(D,1)];
        tryK=[tryK,startK];
        error = abs(inertia(itr) - inertia(itr-1));
        if error<=saturation || startK >= (size(DepdScopeVector,2)-1) || inertia(itr)<=saturation
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