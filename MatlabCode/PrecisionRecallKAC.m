function [Clustered, Detector, RavenTable] = PrecisionRecallKAC(simStruct, truth)



% Input -
% Examp: simulation class object pre-processed for similarity matrix
% Truth: Validated Raven selection table



cluster_ids = unique(simStruct.Cluster_id);

%% Create a raven table from the detections

% Convert arrival time to seconds
start_times =simStruct.arrivalArray(:,1);


scores =simStruct.RefScores;
ClusterId = simStruct.Cluster_id;

% Validation info 1 is target, others are not. In graywhale dataset call
% type 2 is fin and call type 0 is noise. 
Pruned = simStruct.pruned;
dex = simStruct.dex(:,10);
RavenTable = table(...
    start_times,...
   scores,...
    ClusterId,...
    Pruned,...
    dex,...
    'VariableNames',{ ...
    'BeginS','Scores','ClusterId', 'pruned','dex'});



% Avoid division by 0
RavenTable.voting= zeros([height(RavenTable),1])/0;


for jj = 1:length(cluster_ids)
    
    clus = cluster_ids(jj);
    
    
    
    corrScores = RavenTable.Scores(RavenTable.ClusterId==clus);
    corrScores = corrScores(~isnan(corrScores));
%     
%              LR =  log(prod(corrScores./(1.00001-corrScores)))/length(corrScores);
%              LR =prod(corrScores)^(1/length(corrScores));
%              LR = max(corrScores);
%              LR = log((corrScores./(1.00001-corrScores)))
%              
%              log10(prod(corrScores))-log10(prod(1.00001-corrScores))
%     %
     LR = mean(corrScores);
    
    RavenTable.voting(RavenTable.ClusterId==clus)=LR;
    
    %         RavenTable.LRscores(RavenTable.ClusterId==clus) =...
    %             log((corrScores./(1.00001-corrScores)));
%     
    RavenTable.LRscores(RavenTable.ClusterId==clus) =LR;
    
    
    
    
    % Calculate the entropy of class labels within each cluster
    prunedvals= ( RavenTable.pruned(RavenTable.ClusterId==clus));
    prunedIds = unique(prunedvals);
    
    
    
    
end

% Update the truth scores
newRef = RavenTable.LRscores
%newRef =  mean([keepScores,detectionsRef.Score],2, 'omitnan')

simStruct.RefScores=  newRef
% Updatet the 'truth' scores
[simStruct.pruned newTruth]= validateDetection(truth, simStruct, 10);



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
FN  =sum(bsxfun(@lt,newTruth.DetectorScore,thresh)); 

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


% Likelihood ratio
thresh= unique((RavenTable.LRscores))';
TP = sum(bsxfun(@gt,RavenTable.LRscores,thresh).*...
    repmat(RavenTable.pruned,[1, length(thresh)]))
% Detection scores greater than the threshold but not associated with a
% validated call
FP = sum(bsxfun(@gt,RavenTable.LRscores,thresh).*...
    repmat(~RavenTable.pruned,[1, length(thresh)]));

% Calls in the truth table that were not assigned a score
FN  =sum(bsxfun(@lt,truth.DetectorScore,thresh)) 


%%


prec = TP./(TP+FP);
Recall =TP./(TP+FN);
% % 
%   hold on;plot(Recall,prec)
%   xlabel('Recall');ylabel('Precision')


Clustered= struct();
Clustered.Precision = prec;
Clustered.Recall =Recall;
Clustered.Thresh = thresh';
Clustered.FPR = FP;
Clustered.FN = FN;



%% Figures

% 
% 
% figure
% labelSpp=[{'False Positive'}, {'Right Whale'}];
% idxvals =[0 1];
% for ii=1:length(idxvals)
%     
%     idx = find(RavenTable.pruned == idxvals(ii));
%     scores = RavenTable.Scores(idx);
%     
%     subplot(1,2,ii); hist((scores),[0:.05:1]);
%     title([labelSpp(ii)]);
%     xlabel('Neural Network Detection Scores');
%     median(scores)
%     xlim([0.4 1])
% end
% 
% figure
% labelSpp=[{'False Positive'}, {'Right Whale'}];
% idxvals =[0 1];
% for ii=1:length(idxvals)
% 
%     idx = find(RavenTable.pruned == idxvals(ii));
%     scorescChanted = RavenTable.Scores(idx)-RavenTable.LRscores(idx);
% 
%     subplot(1,2,ii); hist((scorescChanted(scorescChanted~=0)));
%     title([labelSpp(ii)]);
%     xlabel('Neural Network Detection Scores');
%     mean(scorescChanted(scorescChanted~=0))
% 
% end
% 
% aa =[ RavenTable.Scores(idx) RavenTable.LRscores(idx)]
% aa(:,3)=aa(:,1)>aa(:,2)
% aa(:,4)=aa(:,1)<aa(:,2);

end