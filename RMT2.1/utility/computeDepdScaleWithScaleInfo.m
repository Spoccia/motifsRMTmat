function [depdScale1] = computeDepdScaleWithScaleInfo(feature1, gss1, idm1, featureDepdScale, featureTimeScale)
% compute match score using matlab
% feature1, feature2, descr1, descr2 -- same size for column
depdSize = size(idm1{1}, 1);
featureSize = size(feature1, 2);
depdScale1 = zeros(depdSize, featureSize);

for i = 1 : size(feature1, 2)
    frame1 = feature1(1:7, i);
    depdin1 = depdscale(gss1, frame1(:,1), idm1, featureDepdScale(i), featureTimeScale(i));
    depdScale1(1:size(depdin1,2), i) = depdin1(1,:)';
end

function depdin = depdscale(gss1, frame, idm1, DepdScale, TimeScale)
% [x(:)' ; y(:)' ; sigmad(:)' ;sigmat(:)' ; octave_Depd--oframes(4,:); octave_Time--oframes(5,:); pricurRatio
% gss.ds{toctave, doctave}
% gss.octave{toctave, doctave}
DOctave = frame(5, 1);
TOctave = frame(6, 1);
% TOctave = frame(3, 1);
OriginalIDM = idm1{1,1};
ActIDM =idm1{1,DOctave};
S1 = gss1.smoothmatrix{DOctave, TOctave}(:,:, DepdScale);
% ignore small number
tempDepdIndex = frame(1, 1);%+1;% Silv Change +1

% if DOctave > 1
%     tempDepdIndex = floor (tempDepdIndex * (1/2)^DOctave);
% end
tempVector = S1(max(1, tempDepdIndex),:);
tempOnes = ones(1, size(tempVector, 2));
tempCombine = tempVector + tempOnes;
if(DOctave == 1)
    temp(1,1) = frame(1, 1);%+1;% Silv add +1
    % temp(2,1) = frame(2, 1); % mark the feature center time 
else
   partialcombination = find(tempCombine ~= 1.000);
   indexoriginalpear=ismember(ActIDM,partialcombination);
     temp(1,:) = OriginalIDM(indexoriginalpear);
    % temp(2,:) = max(1,floor(frame(2,1))); % mark the feature center time
end
depdin = temp;
