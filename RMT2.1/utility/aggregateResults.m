%% AGGREGATE RESULTS FORM SICONG

% Path= 'D:\Motif_Results\Datasets\SynteticDataset\Features_RMT_new\ComputedAccuracy\Results_Motif_Save_NewlyAdded_Excluding_3_and_6\Motif1\Strategy_4\';
Path= 'D:\Motif_Results\Datasets\SynteticDataset\Features_RMT_new\ComputedAccuracy\Results_Motif_Save_NewlyAdded\Motif1\Strategy_3\';
percentage=[0,0.1,0.5,0.75,1];
timeoverlap=[0.1,0.25,0.5,0.75];
AllPrecision=[];
AllRecall=[];
Allfscore=[];
for i= 1:length(percentage)
    Precision=[];
    Recall=[];
    fscore=[]
    for j= 1:length(timeoverlap)
        foldername =['amp_scale_',num2str(percentage(i)),'_TO_',num2str(timeoverlap(j)),'/'];
        Pname='RMTPrecision_aggregated.csv';
        Rname ='RMTRecall_aggregated.csv';
        Fname= 'RMTFScore_aggregated.csv';
        Precisionid = csvread([Path,foldername,Pname]);
        Recallid = csvread([Path,foldername,Rname]);
        Fscoreid = csvread([Path,foldername,Fname]);
        
        if j==1
            Precision=[Precision,Precisionid(:,1:2)];
            Recall=[Recall,Recallid(:,1:2)];
            fscore=[fscore,Fscoreid(:,1:2)];
        else
            Precision=[Precision,Precisionid(:,2)];
            Recall=[Recall,Recallid(:,2)];
            fscore=[fscore,Fscoreid(:,2)];
        end
    end
    AllPrecision{i}=Precision;
    AllRecall{i}=Recall;
    Allfscore{i}=fscore;
end