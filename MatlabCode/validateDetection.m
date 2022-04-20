function [pruned, truth]= validateDetection(truth, simStruct, ref_chan)


truth.Dex = zeros(height(truth), 1);

% Link validation spreadsheet to the detector output


% Create a validation column (T.H. code calls this 'pruned' sticking
% with the naming scheme)
pruned = zeros(size(simStruct.arrivalArray,1),1);

% Create temporary subset of the truth labels
truthSub = truth(truth.Channel==ref_chan,:);

detected=0;
% Step through the calls with TDOA info and if they match a validation
% call assign 1 to the pruned value if it's a rw and 0 otherwise.
for jj=1:size(simStruct.arrivalArray,1)
       
    callStart = simStruct.arrivalArray(jj,1);
    
    % See if any calls are within .3 of a seconds from a truth label
    [diffval, diffIdx] = min(abs(truthSub.BeginTime_s_ - callStart));
  
    
    % If a detection correlates with the truth value add the pruned fx
    if diffval<1.5
        %disp('blarg')
        
        truth.DetectorScore(truth.Selection ==...
            truthSub.Selection(diffIdx))= simStruct.RefScores(jj);
        
        
        pruned(jj) = 1;
        
        %             truthSub(diffIdx,:)=[];
        jj
        truth.Dex(truth.Selection ==...
            truthSub.Selection(diffIdx))=simStruct.dex(jj,10);
    end
end




thresh = unique(round(simStruct.RefScores,2))';
% precision = sum(bsxfun(@gt,detectionsOut.Score,thresh).*repmat(detectionsOut.TP,[1, length(thresh)]))/height(detectionsOut)
% recall=sum(bsxfun(@lt,truthSub.Score,thresh))/height(truthSub)


% Detection scores greater than the thresh and represent a true positive
TP = sum(bsxfun(@gt,simStruct.RefScores,thresh).*repmat(pruned,[1, length(thresh)]));
   

% Detection scores greater than the threshold but not associated with a
% validated call
FP = sum(bsxfun(@gt,simStruct.RefScores,thresh).*repmat(~pruned,[1, length(thresh)]));

% Calls in the truth table that were not assigned a score
FN  =sum(bsxfun(@lt,truth.DetectorScore,thresh));

Prec = TP./(TP+FP);
Recall =TP./(TP+FN);
% 
% figure;plot(Recall,Prec)
% xlabel('Recall')
% ylabel('Precision')





end
