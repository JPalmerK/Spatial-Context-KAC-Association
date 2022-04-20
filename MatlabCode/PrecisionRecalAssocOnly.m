function [Detector, RavenTable] = PrecisionRecalAssocOnly(simStruct, truth)



% Input -
% Examp: simulation class object pre-processed for similarity matrix
% Truth: Validated Raven selection table



%% Create a raven table from the detections

% Convert arrival time to seconds
start_times =simStruct.arrivalArray(:,1);


scores =simStruct.RefScores;

% Validation info 1 is target, others are not. In graywhale dataset call
% type 2 is fin and call type 0 is noise. 
Pruned = simStruct.pruned;
dex = simStruct.dex(:,10);
RavenTable = table(...
    start_times,...
   scores,...
    Pruned,...
    dex,...
    'VariableNames',{ ...
    'BeginS','Scores','pruned','dex'});



% Avoid division by 0
RavenTable.voting= zeros([height(RavenTable),1])/0;


RavenTable.voting(RavenTable.voting==Inf) = ...
    max(RavenTable.voting(isfinite(RavenTable.voting)));




% Detection scores greater than the thresh and represent a true positive
thresh= unique((RavenTable.Scores))';

TP = sum(bsxfun(@gt,RavenTable.Scores,thresh).*...
    repmat(RavenTable.pruned,[1, length(thresh)]))
% Detection scores greater than the threshold but not associated with a
% validated call
FP = sum(bsxfun(@gt,RavenTable.Scores,thresh).*...
    repmat(~RavenTable.pruned,[1, length(thresh)]));

% Calls in the truth table that were not assigned a score
FN  =sum(bsxfun(@lt,truth.DetectorScore,thresh)); 

Prec = TP./(TP+FP);
Recall =TP./(TP+FN);
% % 
%  figure;plot(Recall,Prec)
% xlabel('Recall')


Detector= struct();
Detector.Precision = Prec;
Detector.Recall =Recall;
Detector.Thresh = thresh';
Detector.FPR = FP;
Detector.FN = FN;

end