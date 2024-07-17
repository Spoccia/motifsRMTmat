

function[TIMEFOROCTAVE] = MoCap_FE_start(PATH_dataset,saveFeaturesPath,TS_name, USER_OT_targhet, USER_OD_targhet)

    %% Parameters for  feature extractions
    DeOctTime = USER_OT_targhet;
    DeOctDepd = USER_OD_targhet;
    DeLevelTime = 4;%6;
    DeLevelDepd = 4;%6;
    DeSigmaDepd = 0.5;%0.6;%0.5;%0.4;%
    DeSigmaTime = 4*sqrt(2)/2;%1.6*2^(1/DeLevelTime)*2;%4*sqrt(2);%1.6*2^(1/DeLevelTime);%4*sqrt(2);%2*1.6*2^(1/DeLevelTime);%  8;%4*sqrt(2);%1.2*2^(1/DeLevelTime);%
    DeGaussianThres = 0.3;%0.3;%0.1;%0.4;%1;%0.6;%2;%6; % TRESHOLD with the normalization of hte distance matrix should be  between 0 and 1
    thresh = 0.04 / DeLevelTime / 2 ;%0.04;%
    DeSpatialBins = 4; %NUMBER OF BINs
    r= 10; %5 threshould variates

    %% Location sensory data
    RELATION=csvread(strcat(PATH_dataset,'location\LocationSensor_aggregate.csv'));;%

    %% START CODE FOR FEATURE EXTRACTION
    %load([PATH_dataset,TS_name,'.mat']);
    data = csvread([PATH_dataset,TS_name,'.csv']);
    data(isnan(data))=0;
    if(exist(saveFeaturesPath,'dir')==0)
        mkdir(saveFeaturesPath);
        mkdir([saveFeaturesPath,'GaussianSmoothing\']);
    end

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

    [frames1,descr1,gss1,dogss1,depd1,idm1, time, timee, timeDescr] = sift_gaussianSmooth_Silv(data',RELATION, DeOctTime, DeOctDepd,...
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
    frame1(7,:) = [];
    feature = frame1;
    TIMEFOROCTAVE=toc(p);
    savepath1 = [saveFeaturesPath,'feature_',TS_name,'.mat'];
    savepath2 = [saveFeaturesPath,'idm_',TS_name,'.mat'];
    savepath3 = [saveFeaturesPath,'MetaData_',TS_name,'.mat'];
    savepath5 = [saveFeaturesPath,'GaussianSmoothing/DepdMatrix_',TS_name,'.mat'];
    
    savepath6 = [saveFeaturesPath,'/ComparisonTime_',TS_name,'.csv'];
    savepath7 = [saveFeaturesPath,'/ScaleTime_',TS_name,'.csv'];
    savepath8 = [saveFeaturesPath,'/DescrTime_',TS_name,'.csv'];
    
    save(savepath1,'data', 'gss1', 'frame1','depd1');
    save(savepath2,'idm1');
    save(savepath3,'DeOctTime', 'DeOctDepd', 'DeSigmaTime','DeSigmaDepd', 'DeLevelTime','DeLevelDepd', 'DeGaussianThres', 'DeSpatialBins', 'r', 'descr1' );
end