function [features,dpscale,tostore] = pruningEntropyThresh(features,dpscale,EntropyTHR,data)
    timescope= features(4,:)*3;
    iii=1;
    EachFentropy=zeros(2,size(features,2));
    for iii =1 : size(features,2)
%     while iii <= size(features,2)
        CandidateEntropicFeature=features(:,iii);
        CandidateDepScale=dpscale(:,iii);
        CandidateTimerange=(round((CandidateEntropicFeature(2,1)-timescope(1,iii))) : (round((CandidateEntropicFeature(2,1)+timescope(1,iii)))));
        dataF1 = globalQuantization(data(CandidateDepScale(CandidateDepScale>0,1),CandidateTimerange((CandidateTimerange>0 & CandidateTimerange<=size(data,2)))));
        EntropyF1 = entropy(dataF1/max(dataF1(:)));%EntropySingVariate_mex(dataF1',-Inf);
        EachFentropy(1,iii)=EntropyF1;
        EntropyF1 = EntropySingVariate_mex(dataF1',-Inf);
        EachFentropy(2,iii)=EntropyF1;
%         if(EntropyF1<EntropyTHR)
%             features(:,iii)=[];
%             dpscale(:,iii)=[];
%             timescope(:,iii)=[];
%             iii=iii-1;
%         end
    end
    maximum = max(EachFentropy(2,:));
    thr = maximum*EntropyTHR;
    toremove = EachFentropy(2,:)>=thr;
    tostore = ones(1,size(features,2))-toremove;
    features = features(:,EachFentropy(2,:)<=thr);
    dpscale = dpscale(:,EachFentropy(2,:)<=thr);
%     descriptors = features(11:end,:);
%     Qdescr= globalQuantization(descriptors);
%     Edescriptors=zeros(1,size(features,2));
%     for i = 1: size(features,2)
%         Edescriptors(i) = EntropySingVariate_mex(Qdescr(:,i),-Inf);
%     end
    
end