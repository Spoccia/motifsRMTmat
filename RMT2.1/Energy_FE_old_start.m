

function[TIMEFOROCTAVE] = Energy_FE_old_start(PATH_dataset,saveFeaturesPath,TS_name, USER_OT_targhet, USER_OD_targhet)


    %% Parameters for  feature extractions
    DeOctTime = USER_OT_targhet;
    DeOctDepd = USER_OD_targhet;
    DeLevelTime = 4;%6;%
    DeLevelDepd = 4;%6;%
    DeSigmaDepd = 0.5;%0.6;%0.5;%0.3;%
    DeSigmaTime = 4*sqrt(2)/2;%1.6*2^(1/DeLevelTime);%*2;%4*sqrt(2);%1.6*2^(1/DeLevelTime);%4*sqrt(2);%2*1.6*2^(1/DeLevelTime);%  8;%4*sqrt(2);%1.2*2^(1/DeLevelTime);%
    thresh = 0.04 / DeLevelTime / 2 ;%0.04;%
    DeGaussianThres = 0.1;%0.1;%0.4;%1;%0.6;%2;%6; % TRESHOLD with the normalization of hte distance matrix should be  between 0 and 1
    DeSpatialBins = 4; %NUMBER OF BINs
    r= 10; %5 threshould variates

    LocM1 =[0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
        ];
    LocM2 = zeros(9,9);
    % LocM2=[0	1	1	0	0	0	0	0	0	0	0	0
    %        1	0	0	0	0	0	0	0	0	0	0	0
    %        1	0	0	0	0	0	0	0	0	0	0	0
    %        0	0	0	0	1	1	0	0	0	0	0	0
    %        0	0	0	1	0	0	0	0	0	0	0	0
    %        0	0	0	1	0	0	0	0	0	0	0	0
    %        0	0	0	0	0	0	0	1	1	0	0	0
    %        0	0	0	0	0	0	1	0	0	0	0	0
    %        0	0	0	0	0	0	1	0	0	0	0	0
    %        0	0	0	0	0	0	0	0	0	0	1	1
    %        0	0	0	0	0	0	0	0	0	1	0	0
    %        0	0	0	0	0	0	0	0	0	1	0	0];
    LocM3=[];
    IDM1=1:27;
    % IDM2= [1,2,3,1,2,3,4,5,6,7,8,9,4,5,6,4,5,6,7,8,9,10,11,12,10,11,12];
    IDM2 =[1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9];
    IDM3=[];
    
    idm2{1} = IDM1;
    idm2{2} = IDM2;
    idm2{3} = IDM3;
    
    %% START CODE FOR FEATURE EXTRACTION
    data = csvread([PATH_dataset,TS_name,'.csv']);
    data(isnan(data))=0;
    
    sBoundary=1;
    eBoundary=size(data',1);
    frames1=[];
    descr1=[];
    gss1=[];
    dogss1=[];
    depd1=[];
    idm1=[];
    time=[];
    timee=[];
    timeDescr=[];
    p=tic;
    
    %     [frames1,descr1,gss1,dogss1,depd1,idm1, time, timee, timeDescr] = ...
    %         sift_gaussianSmooth_BirdSong_orgDataDescr(data',...
    %         LocM1 ,LocM2,LocM3,IDM1, IDM2, IDM3, DeOctTime, DeOctDepd,...
    %         DeLevelTime, DeLevelDepd, DeSigmaTime ,DeSigmaDepd,...
    %         DeSpatialBins, DeGaussianThres, r, sBoundary, eBoundary);
    
    %% this commented code will run the RMT extraction from all octaves
    p=tic;
    [frames1,descr1,gss1,dogss1,depd1,idm1, time, timee, timeDescr] = ...
        sift_gaussianSmooth_Proximity(data',...
        LocM1 ,LocM2,LocM3,IDM1, IDM2, IDM3, DeOctTime, DeOctDepd,...
        DeLevelTime, DeLevelDepd, DeSigmaTime ,DeSigmaDepd,...
        DeSpatialBins, DeGaussianThres, r, sBoundary, eBoundary);
    
    
    while(size(frames1,2)==0)
        frames1 = zeros(4,1);
        descr2 = zeros(128,1);
    end
    frame1 = [frames1;descr1];
    if( isnan(sum(descr1(:))))
        TS_name
        nanIDX=  isnan(sum(descr1));
        frame1(:,nanIDX)  = [];
        descr1(:,nanIDX)  = [];
        frames1(:,nanIDX) = [];
    end
    if size(frame1,1)>10
        frame1(7,:) = [];
    end
    feature = frame1;
    TIMEFOROCTAVE=toc(p);
    savepath1 = [saveFeaturesPath,'feature_',TS_name,'.mat'];
    savepath2 = [saveFeaturesPath,'idm_',TS_name,'.mat'];
    savepath3 = [saveFeaturesPath,'MetaData_',TS_name,'.mat'];
    savepath5 = [saveFeaturesPath,'GaussianSmoothing/DepdMatrix_',TS_name,'.mat'];
    savepath6 = [saveFeaturesPath,'/ComparisonTime_',TS_name,'.csv'];
    savepath7 = [saveFeaturesPath,'/ScaleTime_',TS_name,'.csv'];
    savepath8 = [saveFeaturesPath,'/DescrTime_',TS_name,'.csv'];
    
    save(savepath1, 'data','gss1', 'frame1','depd1');%
    save(savepath2,'idm1');
    save(savepath3,'DeOctTime', 'DeOctDepd', 'DeSigmaTime','DeSigmaDepd', 'DeLevelTime','DeLevelDepd', 'DeGaussianThres', 'DeSpatialBins', 'r', 'descr1' );
    save(savepath5, 'depd1');
end