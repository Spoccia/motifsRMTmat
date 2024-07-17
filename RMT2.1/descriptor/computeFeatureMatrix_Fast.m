function [fgss,VariatePosScale] = computeFeatureMatrix_Fast(octave,sminT,sminD, sigmaDep0, St,Sd,H,H1,directed,uniqueComposedVector)

%uniqueComposedVector contains vector of scale and variate for the features extracted
% Index offset
soT = -sminT+1 ;
soD = -sminD+1 ;

[M_Time, N_variate,NSD,NST]= size(octave(:,:,:,:) );
NumFeatureMetrix = size(uniqueComposedVector,2);
VariatePosScale=zeros(1,NSD);
StartPos =1;
for idScaleCreated = 1: VariatePosScale
    Fscale   = uniqueComposedVector(2,idScaleCreated);
    Fvaraite = uniqueComposedVector(1,idScaleCreated);
    
    ScaleSigma = sigmaDep0 * 2^(Fscale/Sd);
    ScopeDependency = ceil(3* ScaleSigma);
    
    V_i_rearrange= zeros(M_Time,2*ScopeDependency+1);
    StartPos= floor(ScopeDependency+1);
    VariatePosScale(1,Fscale+soD)=StartPos;
    V_i_rearrange(:,StartPos)= octave(:,Fvaraite,Fscale+soD,Fscale+soT);
    
    fgss{Fscale+soD,Fscale+soT,Fvaraite} =zeros(M_Time, N_variate);
    
    for i= 1:ScopeDependency
        % construct matrix on right side
        temp = (H^(i)*octave(:,:,Fscale+soD,Fscale+soT)')';
        V_i_rearrange(:,StartPos+i)= temp(:,Fvaraite);
        if(directed==0) % undirected graph
            temp1= (H1^(i)*octave(:,:,Fscale+soD,Fscale+soT)')';
            V_i_rearrange(:,StartPos-i)= temp1(:,Fvaraite);
        else% directed graph
            V_i_rearrange(:,StartPos-i)= temp(:,Fvaraite);
        end
    end
    fgss{Fscale+soD,Fscale+soT,Fvaraite}=V_i_rearrange;
end