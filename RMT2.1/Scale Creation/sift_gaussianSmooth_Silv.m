function [tempFrames,descriptor_silv,gss,dogss,depd,idm, TIMESCALE, TIMEKEYPOINTS, TIMEDESCRIPTORS]=sift_gaussianSmooth_Silv(I, LocM, Ot, Od, St, Sd, sigmaTime ,sigmaDepd, NBP, gthresh, r,sBoundary, eBoundary)
%, featureDepdScale, featureTimeScale]
% [M,N,C] = size(I) ;
% O = floor(log2(min(M,N)))-2 ; % up to 8x8 images
% time  = zeros(1, Ot*Od);
time = zeros(Ot, Od);
timeDescr = zeros(Ot, Od);
timee = zeros(1,2);

featureTimeScale = [];
featureDepdScale = [];
% thresh = 0.04 / St / 2 ;
thresh = 0.04 / 3 / 2 ; %value picked from the vidaldi code
NBO    = 8;
NBP_Time = 4;
NBP_Depd = 4;
magnif = 3.0;
% frames      = [] ;
tempFrames = [];
descriptors = [] ;
descriptor_silv=[];
ktime = 2^(1/(St-3));
kdepd = 2^(1/(Sd-3));

%p = tic;
% Compute scale spaces
% [gss, depd, idm] = gaussianss_asynchronousMote(I, LocM, Ot, Od,St, Sd, sigmaTime, sigmaDepd, gthresh);
%[gss, depd, idm] = gaussianss_asynchronousMote_Silv_1(I, LocM, Ot,Od,St,Sd, sigmaTime, sigmaDepd, gthresh);%this was in use

% Try this function
%[gss, depd, idm] = gaussianss_asynchronous_SilvRewrite(I, LocM, Ot, Od,St,Sd, sigmaTime, sigmaDepd, gthresh);
stmin=0;%-1;
sdmin=0;%-1;
otmin=0;
odmin=0;

%   [gss, depd, idm,TIMESCALE] = gaussianss_Silv_fromorg(I, LocM, Ot,Od,St,Sd,otmin,odmin,stmin,sdmin,St+1,Sd+1, sigmaTime, sigmaDepd,gthresh,-1,-1);
  [gss, depd, idm,TIMESCALE] = gaussianss_Silv_fromorg_justUserOctave(I, LocM, Ot,Od,St,Sd,otmin,odmin,stmin,sdmin,St+1,Sd+1, sigmaTime, sigmaDepd,gthresh,-1,-1);

%[gss, depd, idm] = gaussianss_asynchronousMote_Silv_2(I, LocM, Ot, Od,St, Sd, sigmaTime, sigmaDepd, gthresh);
% timee(1) = 0;%timee(1)+ toc(p);
TIMEKEYPOINTS=zeros(1,4);
% counter=1;

%p = tic;
% dogss = diffss_asynchronous(gss);
dogss = diffss_asynchronous_justUserOctave(gss,Od,Ot);

% dogss = diffss_asynchronousTest(gss);
% timee(2) = timee(2)+ toc(p);

otime = Ot;
odepd = Od;
% for otime = 1: size(gss.octave,1)
%     for  odepd = 1: size(gss.octave,2)
%         p = tic;
        tic;
        % Local maxima of the DOG octave
        scaleDiff = St - 1;
        % forwardIdx = siftlocalmax_directed_bak12142015(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
        % forwardIdx = siftlocalmax_directed_95(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
        %% used
             forwardIdx = siftlocalmax_directed_100(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
        %% optimized
            % forwardIdx = siftlocalmax_directed_sicong_optimize(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
        % forwardIdx = siftlocalmax_directed_999(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
        % forwardIdx = siftlocalmax_directed_998(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
        [i,j, s1] = ind2sub( size( dogss.octave{otime, odepd}{3}), forwardIdx ) ;
        y=i-1;
        x=j-1;
        s1=s1-1+gss.sminT;
        % s1=s1-1-1 ; SICONG
        s2 = s1;
        forwardIdx = [x(:)';y(:)';s1(:)'] ;
        
        dogss.octave{otime, odepd}{1} = -dogss.octave{otime, odepd}{1};
        dogss.octave{otime, odepd}{2} = -dogss.octave{otime, odepd}{2};
        dogss.octave{otime, odepd}{3} = -dogss.octave{otime, odepd}{3};
        
        % backwardIdx = siftlocalmax_directed_bak12142015(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
        % backwardIdx = siftlocalmax_directed_95(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
        %% used
             backwardIdx = siftlocalmax_directed_100(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
        %% optimized
            % backwardIdx = siftlocalmax_directed_sicong_optimize(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
        % backwardIdx = siftlocalmax_directed_999(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
        % backwardIdx = siftlocalmax_directed_998(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);
%         time(otime, odepd) = 0;% time(otime, odepd) + toc(p);
        [i,j, s1] = ind2sub( size( gss.octave{otime, odepd} ), backwardIdx ) ;
        y=i-1;
        x=j-1;
        s1=s1-1+gss.sminT;
        % s1=s1-1-1;% for the new gaussian : +dogss.sminT SICONG
        s2 = s1;
        backwardIdx = [x(:)';y(:)';s1(:)'] ;
        % backwardIdx = [x(:)';y(:)';s1(:)';s2(:)'] ;
        
        oframes = [forwardIdx, backwardIdx];
        [C, ia, ic] = unique(oframes', 'rows');
        oframes = C';
        for hsize = 1:size(dogss.octave{otime, odepd}{3},3)% iterate on the scale in the specific octave
            HY(:,:,hsize) = (NormalizeByRow(depd{odepd})*(-dogss.octave{otime, odepd}{3}(:,:,hsize)'))';
            HY2(:,:,hsize) = (NormalizeByRow(depd{odepd}')*(-dogss.octave{otime, odepd}{3}(:,:,hsize)'))';
        end
        rad = 0;
        sel= ... % select feature that are centered in the timeseries
            oframes(2,:)-rad >= 1  & ...
            oframes(2,:)+rad <= size(gss.octave{otime, odepd},1)      ;
        oframes=oframes(:,sel) ;
        % pruning Features Extrcted on the base of curvature
        %   ofra =oframes;
        oframes = siftrefinemx_directed(oframes, -dogss.octave{otime, odepd}{3},HY,HY2,gss.sminT,thresh,r,0) ;
        % this is my function more allined with the original code
        %  oframes = siftrefinemx_directed_Silv(oframes, -dogss.octave{otime, odepd}{3},HY,HY2,gss.sminT,thresh,r,0) ;
        % The last term 0 I think refer to the starting octave in hte original
        % code there is no that term
        % Silv concern: If wedo not do this we should use 0 curvature
        % adding a 0 vector to the oframes.
        TIMEKEYPOINTS(otime + odepd)= toc;
        clear HY HY2
        TIMEDESCRIPTORS=zeros(1,otime + odepd);
        tic;
        if size(oframes,2) >0
            pricurRatio = oframes(4,:);%zeros(1,size(oframes,2));%
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

% %             featureDepdScale = [featureDepdScale dependencyScale];
% %             featureTimeScale = [featureTimeScale timeScale];
            
            %timeScale = oframes(3,:)+1; % Silv SAY: ASK Sicong why this line of code is like this... It seems correct this are not the octave
            % sigma = 2^(o-1+gss.omin) * gss.sigma0 * 2.^(oframes(3,:)/gss.S) ;
            sigmad =  2^(odepd-1+gss.odmin) * gss.sigmad * 2.^(oframes(3,:)/gss.Sd);
            %sigmad =  sigmaDepd*2^(odepd-1)*kdepd.^(dependencyScale-1) ;%gss.sigmat = sigmatimezero  SICONG
            sigmat =  2^(otime-1+gss.odmin) * gss.sigmat * 2.^(oframes(3,:)/gss.Sd);
            %sigmat =  (sigmaTime*2^(otime-1))*ktime.^(timeScale-1);%gss.sigmat = sigmatimezero SICONG
            
            % append difference-of-Gausssian values to output
            TimescaleSicongNormalized = timeScale+1;
            [timeDoGs, depdDoGs, bothDoGs] = appendDogs(dogss.octave{otime, odepd}, tempDepd, tempTime, dependencyScale, TimescaleSicongNormalized);

            %% unique goes here
%           tempFrames = [tempFrames, [x(:)'+ones(1,size(x,1)) ; y(:)' ; sigmad(:)' ;sigmat(:)' ; oframes(4,:); oframes(5,:); pricurRatio; timeDoGs(:)'; depdDoGs(:)'; bothDoGs(:)']];%[x(:)'
            tempFrames = [tempFrames, [x(:)'+ones(1,size(x,1)) ; y(:)' ; sigmad(:)' ;sigmat(:)' ; oframes(4,:); oframes(5,:); oframes(3,:);pricurRatio; timeDoGs(:)'; depdDoGs(:)'; bothDoGs(:)']];%[x(:)'

        end
        
        dogss.octave{otime, odepd}{1} = -dogss.octave{otime, odepd}{1};
        dogss.octave{otime, odepd}{2} = -dogss.octave{otime, odepd}{2};
        dogss.octave{otime, odepd}{3} = -dogss.octave{otime, odepd}{3};
        % 1 means directed graph
        [fgss_silv,Pseudo_centerVaraite]= computeFeatureMatrix_Silv(gss.octave{otime, odepd},gss.sminT,gss.sminD,gss.sigmad,gss.St,gss.Sd,NormalizeByRow(depd{odepd}),NormalizeByRow(depd{odepd}'),1);
        % Descriptors
        if(size(oframes, 2) > 0)
%             p = tic;

%            fgss = computeFeatureMatrix_directed_bak12112015(gss.octave{otime, odepd}, oframes, NormalizeByRow(depd{odepd}));
            % for f=1:size(oframes,2)
            for f=1:size(oframes,2)
%                 oframe_old = oframes(:,f);
%                 oframe_old(4,:)= [];%remove octave deendency   
%                 oframe_old(4,:)=oframe_old(4,:)*0;
                oframeSilv=oframes(:,f);
                centerV= Pseudo_centerVaraite(1,oframeSilv(3,1)-gss.sminD+1);
                oframeSilv(1,1)=centerV;
                % oframe (3, 1); scale Dependency
                % oframe (4, 1); scale Time
                % oframe (5, 1); eventual orientation value
                oframeSilv(4,1)=oframeSilv(3,1);
                oframeSilv(5,1)=0;
                %                 fgss_silv{oframes(3,1)-gss.sminD+1,oframes(4,1)-gss.sminT+1,oframes(1,f)},...% gss.octave{otime, odepd}(:,:,oframes(3,1)-gss.sminD+1,oframes(4,1)-gss.sminT+1),...%THIS IS THE SCALE OF THE FEATURE F (TIME,VARIATE,SD,ST)
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
                
%                 sh = siftdescriptor_original(...
%                     fgss{f}, ...
%                     oframe_old(:,1),...
%                     gss.sigmat, ...
%                     gss.St, ...
%                     0, ...
%                     'Magnif', magnif, ...
%                     'NumSpatialBins', NBP, ...
%                     'NumOrientBins', NBO) ;
%                 descriptors = [descriptors, sh] ;
            end
          
        end
        TIMEDESCRIPTORS(otime + odepd)= toc;
%             counter=counter+1;
        clear fOframes bOframes fgss
        
%        timeDescr(%otime, odepd) = timeDescr(otime, odepd)+toc(p);
%     end
    
% end
%time = TIMEKEYPOINTSDESCRIPTORS+TIMESCALE;
end
