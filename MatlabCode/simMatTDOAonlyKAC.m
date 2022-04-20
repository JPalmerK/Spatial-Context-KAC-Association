function Sim_mat = simMatTDOAonlyKAC(simStruct)

% Create the simulation matrix using TDOA values only

Sim_mat = nan(size(simStruct.arrivalArray,1));
arrivalArray= (simStruct.arrivalArray);
diffTimeMat = Sim_mat;

for ii =1:(size(simStruct.TDOA_vals,1))
    
    tdoa_orig =simStruct.TDOA_vals(ii,:);
    
    % Figure out the number of time gaps within the maximum
    % allowed correlation time (time_cut)
    time_gaps = arrivalArray(ii:end, 1)-...
        arrivalArray(ii, 1);
    
    diff_times = diff(arrivalArray(ii:end, 1));
    
    % Find first big gap
    idx_end = find(diff_times>= simStruct.maxEltTime,1)-1;
    
    if isempty(idx_end)
        idx_end = length(time_gaps);
    end
    
    time_gaps = time_gaps(1:idx_end);
    TDOA_next = simStruct.TDOA_vals(ii:ii+idx_end-1,:);
    
    
    deltaTDOA = bsxfun(@minus, tdoa_orig,TDOA_next);
    elapsedTime = time_gaps;
    
%     % Binary mask for cross correlation score
%     % Columns to keep
%     keepCol = 1:10; keepCol([4,parent])=[];
%     BinMask = simStruct.arrivalTable.CrossScore(ii:ii+idx_end-1,keepCol);
%     BinMask(BinMask<.8)=nan;BinMask(BinMask>=0.8)=1;
%     deltaTDOA=deltaTDOA.*BinMask;
    
    deltaTDOALklhd=zeros(size(deltaTDOA));
    
    for jj=1:size(deltaTDOA,2)
        
        mu = zeros(size(TDOA_next,1),1);
        sigmaSwim = sqrt(simStruct.drift^2+ simStruct.PosUncertsigma+(elapsedTime *...
            (simStruct.s)/simStruct.c).^2);
        x =deltaTDOA(:,jj); %values
        
        likelihood = normpdfRep(x,mu,sigmaSwim);
        
        % Normalizing factor
        LikelihoodNormFac= normpdfRep(0,0,sigmaSwim);
        NormLikelihood = likelihood./LikelihoodNormFac;
        
        %%%%
% %         NormLikelihood(isnan(NormLikelihood))=1;
%         NormLikelihood(isnan(NormLikelihood))=0;
        %%%
        
  
        
        % Normalized likelihood
        deltaTDOALklhd(:,jj)= NormLikelihood; % sigma
        % Create normalizing factor
        
        
   
    end
    
    simValues = min(deltaTDOALklhd,[],2, 'omitnan' );
    
    % Take the minimum value and fill in the similarity matrix
    
    Sim_mat(ii, ii:ii+length(simValues)-1) = simValues;
    Sim_mat(ii:ii+length(simValues)-1,ii) = simValues;
    Sim_mat(ii,ii) = 1;
    
   
    
end



end