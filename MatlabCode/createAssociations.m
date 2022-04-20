function [delays, dex, crossScores, xcorVal] = ...
    createAssociations(detections,idx,file_listTable, hydDelay, ref_chan, fs)


noisePadLength =10; % seconds on either side of the even to pad

% convert noise buffer from seconds to samples
noisebuffpts = round(noisePadLength * fs);
% get max delay
max_delay = max(hydDelay);
% to acquire the data needed for all calculations to follow
event_pad = max(noisebuffpts, max_delay);

% duration is always 2 sec for NN output
totalpts = ceil(2 * fs) + (2 * event_pad);

%% use the extent of the event to bandlimit
flo = detections.LowFreq_Hz_(1);
fhi = detections.HighFreq_Hz_(1);

detectionsRef = detections(detections.Channel==10,:);

num_chan=max(detections.Channel)-min(detections.Channel)+1;
delays = nan(1, 10);
dex = delays;
crossScores = dex;
arrivalArray = dex;
detectionsRef =detections(detections.Channel==ref_chan,:);


% convert noise buffer from seconds to samples
noisebuffpts = round(noisePadLength * fs);



callStartSamp = round(detectionsRef.BeginTime_s_(idx)*fs);
callStopSamp=   round(detectionsRef.EndTime_s_(idx)*fs);


aa = (((callStopSamp+event_pad)/fs)>detections.BeginTime_s_ & detections.Channel~=ref_chan);
aa1 = (detections.EndTime_s_>((callStartSamp-event_pad)/fs)& detections.Channel~=ref_chan);

inxCheck=find(aa.*aa1);

% If there are any detections in the region load
if ~isempty(inxCheck)
    
    % Following KAC
    % I. Get RL of reference event
    %   1) get event plus BG noise pad
    %   2) separate event chunk from BG noise chunk
    %   3) calculate event STFT
    %   4) calculate noise STFT to get NSE
    %   5) de-noise event STFT
    %   6) calculate RL from de-noised event STFT
    
    [xcorVal] = KACRevisedXcorr(file_listTable,hydDelay, noisebuffpts, ...
        ref_chan, num_chan, flo, fhi, callStartSamp, callStopSamp);
    
    
    % Get the matching calls
    detectionsSub = detections(inxCheck,:);
    chanDetections = unique(detectionsSub.Channel);
    
    for jj =1:length(chanDetections)
        
        currChan = chanDetections(jj);
        detRel =detectionsSub(detectionsSub.Channel==currChan,:);
        
        [timeDiff, indexDiff] =   min(abs(...
            detRel.BeginTime_s_-xcorVal.event_time(chanDetections(jj))));
        
        if timeDiff<.75
            delays(currChan)=xcorVal.event_time(currChan)-callStartSamp/fs;
            crossScores(currChan)=xcorVal.all_pk_xcorr_norm(currChan);
            dex(currChan)=detRel.Selection(indexDiff);
        end
        
        
    end
    
    
else
    xcorVal=[];
    
end



dex(:,ref_chan)=detectionsRef.Selection(idx);







end