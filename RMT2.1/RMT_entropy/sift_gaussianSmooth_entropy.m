function [tempFrames,descriptor_silv,gss,dogss,depd,idm, TIMESCALE, TIMEKEYPOINTS, TIMEDESCRIPTORS]=sift_gaussianSmooth_entropy(I, LocM, Ot, Od, St, Sd, sigmaTime ,sigmaDepd, NBP, gthresh, r,sBoundary, eBoundary)
% [M,N,C] = size(I) ;
% O = floor(log2(min(M,N)))-2 ; % up to 8x8 images
% time  = zeros(1, Ot*Od);
time = zeros(Ot, Od);
timeDescr = zeros(Ot, Od);
timee = zeros(1,2);

featureTimeScale = [];
featureDepdScale = [];
% thresh = 0.04 / St / 2 ;
thresh = 0.04 / 3 / 2 ; % why???
NBO    = 8;
magnif = 3.0;
NBP_Time = 4;
NBP_Depd = 4;
% frames      = [] ;
tempFrames = [];
descriptors = [] ;
descriptor_silv=[];
ktime = 2^(1/(St-3));
kdepd = 2^(1/(Sd-3));

stmin=0;%-1;
sdmin=0;%-1;
otmin=0;
odmin=0;

p = tic;
% Compute scale spaces
% [gss, depd, idm] = gaussianss_asynchronousMote(I, LocM, Ot, Od,St, Sd, sigmaTime, sigmaDepd, gthresh);
% work original Sicong with  graph normalized
%[gss, depd, idm] = gaussianss_asynchronousMote_Silv_1(I, LocM, Ot, Od,St, Sd, sigmaTime, sigmaDepd, gthresh);
% work scale created just with entropy function on the original image
% [gss, depd, idm] = gaussianss_asynchronousMote_Silv_Entropy(I, LocM, Ot, Od,St, Sd, sigmaTime, sigmaDepd, gthresh);
%[gss, depd, idm] = gaussianss_asynchronousMote_Silv_Entropy_smooth(I, LocM, Ot, Od,St, Sd, sigmaTime, sigmaDepd, gthresh);
[gss, depd, idm,TIMESCALE] = EntropyScale_Silv_fromorg_justUserOctave        (I, LocM, Ot,Od,St,Sd,otmin,odmin,stmin,sdmin,St+1,Sd+1, sigmaTime, sigmaDepd,gthresh,-1,-1);
% [gss, depd, idm,TIMESCALE] = gaussianss_Silv_fromorg_justUserOctave(I, LocM, Ot,Od,St,Sd,otmin,odmin,stmin,sdmin,St+1,Sd+1, sigmaTime, sigmaDepd,gthresh,-1,-1);
TIMEKEYPOINTS=zeros(1,4);


dogss = diffss_asynchronous_justUserOctave(gss,Od,Ot);%dogss = diffss_asynchronous(gss); % difference of gaussians
doess = diffss_asynchronous_Entropy_justUserOctave(gss,Od,Ot);%doess = diffss_asynchronous_Entropy(gss);
% dogss = diffss_asynchronousTest(gss);
otime = Ot;
odepd = Od;
% % oframes_Sicong=[];
% % oframes_likeSicong=[];
% % oframes_EntropyScale=[];
tic;
% Local maxima of the DOG octave
scaleDiff = St - 1;

forwardIdx_Entropy = siftlocalmax_directed_100(doess.Entropyoctave{otime, odepd}{3},...
    doess.Entropyoctave{otime, odepd}{2},...
    doess.Entropyoctave{otime, odepd}{1},...
    0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);


[i,j, s1] = ind2sub( size( gss.octave{otime, odepd} ), forwardIdx_Entropy ) ;
y=i-1;
x=j-1;
s1=s1-1+gss.sminT;
s2 = s1;
forwardIdx_Entropy = [x(:)';y(:)';s1(:)'] ;


dogss.octave{otime, odepd}{1} = -dogss.octave{otime, odepd}{1};
dogss.octave{otime, odepd}{2} = -dogss.octave{otime, odepd}{2};
dogss.octave{otime, odepd}{3} = -dogss.octave{otime, odepd}{3};
%         backwardIdx = siftlocalmax_directed_100(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
%
doess.Entropyoctave{otime, odepd}{3}= -doess.Entropyoctave{otime, odepd}{3};
doess.Entropyoctave{otime, odepd}{2}= -doess.Entropyoctave{otime, odepd}{2};
doess.Entropyoctave{otime, odepd}{1}= -doess.Entropyoctave{otime, odepd}{1};

backwardIdx_Entropy = siftlocalmax_directed_100(...
    doess.Entropyoctave{otime, odepd}{3},...
    doess.Entropyoctave{otime, odepd}{2},...
    doess.Entropyoctave{otime, odepd}{1},...
    0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);% Entropy scale
% Silv Entropy Sclale detected features
[i,j, s1] = ind2sub( size( gss.octave{otime, odepd} ), backwardIdx_Entropy ) ;
y=i-1;
x=j-1;
s1=s1-1+gss.sminT;
s2 = s1;
backwardIdx_Entropy = [x(:)';y(:)';s1(:)'] ;


oframes = [forwardIdx_Entropy, backwardIdx_Entropy];

[C, ia, ic] = unique(oframes', 'rows');
oframes = C';
for hsize = 1:size(dogss.octave{otime, odepd}{3},3)% iterate on the scale in the specific octave
    HY(:,:,hsize) = (NormalizeByRow(depd{odepd})*(-dogss.octave{otime, odepd}{3}(:,:,hsize)'))';
    HY2(:,:,hsize) = (NormalizeByRow(depd{odepd}')*(-dogss.octave{otime, odepd}{3}(:,:,hsize)'))';
end
rad = 0;
if(size(oframes,2)~=0)
    sel= ... % select feature that are centered in the timeseries
        oframes(2,:)-rad >= 1  & ...
        oframes(2,:)+rad <= size(gss.octave{otime, odepd},1)      ;
    oframes=oframes(:,sel) ;
    % Silv concern: If wedo not do this we should use 0 curvature (I do not rememebr )
    % adding a 0 vector to the oframes.
    %pricurRatio = zeros(1,size(oframes,2));  % Silv setted the pricur ratio = 0 ;
    
    % pruning Features Extrcted on the base of curvature
    oframes = siftrefinemx_directed(oframes, -dogss.octave{otime, odepd}{3},HY,HY2,gss.sminT,thresh,r,0) ;
    % Silv commented the pruning step here
    
end
TIMEKEYPOINTS(otime + odepd)= toc;
clear HY HY2
TIMEDESCRIPTORS=zeros(1,otime + odepd);
tic;
if size(oframes,2) >0
    pricurRatio = zeros(1,size(oframes,2));%oframes(4,:);
else
    pricurRatio = zeros(1,0);
end

if(size(oframes, 2) ~=0)
    oframes(4,:) = odepd; % Save the Octave Dependency
    oframes(5,:) = otime; % Save the Octave Time
    % Store frames
    x = oframes(1, :); % the original code report to the variate of the  specific scale 2^(o-1+gss.omin) * oframes(1,:) ;
    y  = 2^(gss.ds{otime, odepd}(1)+gss.otmin-1) * oframes(2,:) ; % otmin starts from 0 for timeseries
    %report tievalue to the original scale
    
    tempDepd = oframes(1,:) ;
    tempTime = oframes(2,:) ;
    dependencyScale = oframes(3,:);
    %dependencyScale = oframes(3,:)+1; %SICONG
    % timeScale = oframes(4,:)+1;
    timeScale = oframes(3,:);% SICONG
    featureDepdScale = [featureDepdScale dependencyScale];
    featureTimeScale = [featureTimeScale timeScale];
    
    %timeScale = oframes(3,:)+1; % Silv SAY: ASK Sicong why this line of code is like this... It seems correct this are not the octave
    % sigma = 2^(o-1+gss.omin) * gss.sigma0 * 2.^(oframes(3,:)/gss.S) ;
    sigmad =  2^(odepd-1+gss.odmin) * gss.sigmad * 2.^(oframes(3,:)/gss.Sd);
    %sigmad =  sigmaDepd*2^(odepd-1)*kdepd.^(dependencyScale-1) ;%gss.sigmat = sigmatimezero  SICONG
    sigmat =  2^(otime-1+gss.odmin) * gss.sigmat * 2.^(oframes(3,:)/gss.Sd);
    %sigmat =  (sigmaTime*2^(otime-1))*ktime.^(timeScale-1);%gss.sigmat = sigmatimezero SICONG
    
    
    % append difference-of-Gausssian values to output
    TimescaleSicongNormalized = timeScale+1;
    [timeDoGs, depdDoGs, bothDoGs] = appendDogs(dogss.octave{otime, odepd}, tempDepd, tempTime, dependencyScale, TimescaleSicongNormalized);
    % unique goes here
    tempFrames = [tempFrames, [x(:)'+ones(1,size(x,1)) ; y(:)' ; sigmad(:)' ;sigmat(:)' ; oframes(4,:); oframes(5,:); oframes(3,:);pricurRatio; timeDoGs(:)'; depdDoGs(:)'; bothDoGs(:)']];%[x(:)'
    %           tempFrames = [tempFrames, [x(:)'+ones(1,size(x,1)) ; y(:)' ; sigmad(:)' ;sigmat(:)' ; oframes(4,:); oframes(5,:); pricurRatio; timeDoGs(:)'; depdDoGs(:)'; bothDoGs(:)']];%[x(:)'
end

dogss.octave{otime, odepd}{1} = -dogss.octave{otime, odepd}{1};
dogss.octave{otime, odepd}{2} = -dogss.octave{otime, odepd}{2};
dogss.octave{otime, odepd}{3} = -dogss.octave{otime, odepd}{3};
% 1 means directed graph
[fgss_silv,Pseudo_centerVaraite]= computeFeatureMatrix_Silv(gss.octave{otime, odepd},gss.sminT,gss.sminD,gss.sigmad,gss.St,gss.Sd,NormalizeByRow(depd{odepd}),NormalizeByRow(depd{odepd}'),1);
% Descriptors
if(size(oframes, 2) > 0)
    p = tic;
    % for f=1:size(oframes,2)
    for f=1:size(oframes,2)
        
        oframeSilv=oframes(:,f);
        centerV= Pseudo_centerVaraite(1,oframeSilv(3,1)-gss.sminD+1);
        oframeSilv(1,1)=centerV;
        % oframe (3, 1); scale Dependency
        % oframe (4, 1); scale Time
        % oframe (5, 1); eventual orientation value
        oframeSilv(4,1)=oframeSilv(3,1);
        oframeSilv(5,1)=0;
      
        sh_silv=siftdescriptor_Silv(...
            fgss_silv{oframeSilv(3,1)-gss.sminD+1,oframeSilv(4,1)-gss.sminT+1,oframes(1,f)+1},...
            oframeSilv(:,1),...%oframes(:,f), ...
            gss.sigmat, ...
            gss.sigmad	,...
            gss.St, ...
            gss.Sd, ...
            gss.sminT	, ...
            gss.sminD, ...
            magnif, ...
            NBP_Time, ...
            NBP_Depd, ...
            NBO) ;
        
        descriptor_silv=[descriptor_silv,sh_silv];
    end
end
TIMEDESCRIPTORS(otime + odepd)= toc;
clear fOframes bOframes fgss
end
