function [timeDoGs, depdDoGs, bothDoGs] = appendDogs(dogssInput, depdIndex, timeIndex, dependencyScale, timeScale)
%% append DoGs values for the frames computed
% input parameters timeIndex, depdIndex, timeScale and dependencyScale have
% the same size
featureNum = size(depdIndex, 2);
timeDoGs = zeros(1, featureNum);
depdDoGs = zeros(1, featureNum);
bothDoGs = zeros(1, featureNum);
for i = 1 : featureNum
%     timeDoGs(1, i) = dogssInput{2}(timeIndex(1, i), depdIndex(1, i), timeScale(1, i), dependencyScale(1, i));
%     depdDoGs(1, i) = dogssInput{2}(timeIndex(1, i), depdIndex(1, i), timeScale(1, i), dependencyScale(1, i));
%     bothDoGs(1, i) = dogssInput{2}(timeIndex(1, i), depdIndex(1, i), timeScale(1, i), dependencyScale(1, i));


%     timeDoGs(1, i) = dogssInput{2}(timeIndex(1, i), depdIndex(1, i), timeScale(1, i));
%     depdDoGs(1, i) = dogssInput{2}(timeIndex(1, i), depdIndex(1, i), timeScale(1, i));
%     bothDoGs(1, i) = dogssInput{2}(timeIndex(1, i), depdIndex(1, i), timeScale(1, i));
    timeDoGs(1, i) = dogssInput{2}(timeIndex(1, i), depdIndex(1, i), timeScale(1, i));
    depdDoGs(1, i) = dogssInput{1}(timeIndex(1, i), depdIndex(1, i), timeScale(1, i));
    bothDoGs(1, i) = dogssInput{3}(timeIndex(1, i), depdIndex(1, i), timeScale(1, i));

end
