
function all = StructureForFiles(path,datasetname,listsize);
   all =[];
    for i=1:listsize 
    %  strcat(path,'\',disrurb,'\',datasetname,'\',num2str(i),'_Percent_',num2str(percent),'_instance_',num2str(instanceN),'.csv')
        SAX_TS_Q = csvread(strcat(path,'\',num2str(i),'.csv'));
        %SAX_TS_Q = SAX_TS_Q(5:12,:);
        all= [all,struct(datasetname,SAX_TS_Q)];
    end
%m= all(2).matix
end