function [fgss,VariatePosScale] = computeFeatureMatrix(octave,sminT,sminD, sigmaDep0, St,Sd,H,H1,directed)
% For each combination of scale  dependency and Time in the octave create a
% structure where  the matrix  represent the neighboorhood of the  specific
% variate V_i
% octave: all the scale combination of the secific octave combination for 
%         higher octave the scale are  reduced in size then we can consider
%         without loss the same sigma scale because the sigma 2 between
%         octave is considered in the halfsize.
% sminT: minimum scale over time
% sminD: maimum scale over time 
% sigmaDep0: here we are precomputing the  matrices needed for each variate for each scale
% St: number of scale time
% Sd: number of scale dependency
% H: forward connectivity row normalized F in the paper
% H1:backward connectivity wor normalized
% directed flag to determine if the  graph used is directed or undirected:
% 1. graph directed the  backward and forward matrices are  equals
% 0. graph undirected the  backward and forward matrices are  different


% Index offset
soT = -sminT+1 ;
soD = -sminD+1 ;

[M_Time, N_variate,NSD,NST]= size(octave(:,:,:,:) );
VariatePosScale=zeros(1,NSD);
StartPos =1;
for si_T= sminT:St+1
    for si_D= sminD:Sd+1
        ScaleSigma = sigmaDep0 * 2^(si_D/Sd);
        ScopeDependency = ceil(3* ScaleSigma);
        for V_i=1:N_variate
            V_i_rearrange= zeros(M_Time,2*ScopeDependency+1);
            StartPos= floor(ScopeDependency+1);
            VariatePosScale(1,si_D+soD)=StartPos;
            V_i_rearrange(:,StartPos)= octave(:,V_i,si_D+soD,si_T+soT);
             fgss{si_D+soD,si_T+soT,V_i} =zeros(M_Time, N_variate);
%             fgss(:,:,si_D+soD,si_T+soT,V_i) = zeros(M_Time, N_variate);
            for i= 1:ScopeDependency
                % construct matrix on right side
                temp = (H^(i)*octave(:,:,si_D+soD,si_T+soT)')';
                V_i_rearrange(:,StartPos+i)= temp(:,V_i);
                if(directed==0) % undirected graph
                  temp1= (H1^(i)*octave(:,:,si_D+soD,si_T+soT)')';
                  V_i_rearrange(:,StartPos-i)= temp1(:,V_i);
                else% directed graph
                  V_i_rearrange(:,StartPos-i)= temp(:,V_i);
                end
            end
           fgss{si_D+soD,si_T+soT,V_i}=V_i_rearrange;
           
        end

    end
end
    