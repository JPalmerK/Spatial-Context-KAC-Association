%% Create a subset of interleaved times

[counts, ids] = groupcounts(aa(:,13));
df =table(counts, ids);
idkeep = ids(counts>30)

D = bsxfun(@eq, (aa(:,13)), idkeep')

aa= [arrivalArray simStruct.pruned];
aa(:,12)= [0;diff(aa(:,10))];
aa(:,13)= simStructAcousticEncounders.Cluster_id;
D = bsxfun(@eq, (aa(:,13)), idkeep')
D = D*1
D(D==0)=nan;
cc = find(prod(D,2, 'omitnan'))


arrayTrimmed = arrivalArray(cc,:);
dexTrimmed = dex(cc,:);
delaysTrimmed = delays(cc,:);


%% Clean out the truth table and the detection table
TruthTrimmed = truth;
TruthTrimmed.Keep = zeros(height(truth),1);
detectionTrimmed =detections;
detectionTrimmed.Keep = zeros(height(detections),1);
for ii= 1:length(idkeep)
    
    keepTime = [min(aa((aa(:,13)==idkeep(ii)),10)) max(aa((aa(:,13)==idkeep(ii)),10))]
    
    TruthTrimmed.Keep(...
        TruthTrimmed.BeginTime_s_>=keepTime(1) &...
        TruthTrimmed.EndTime_s_<=keepTime(2))=1;
    
     detectionTrimmed.Keep(...
        detectionTrimmed.BeginTime_s_>=keepTime(1) &...
        detectionTrimmed.EndTime_s_<=keepTime(2))=1;
    
    
end

TruthTrimmed= TruthTrimmed(TruthTrimmed.Keep==1,:);
detectionTrimmed= detectionTrimmed(detectionTrimmed.Keep==1,:);
%%


simStruct.arrivalArray=arrayTrimmed;
simStruct.TDOA_vals = delaysTrimmed;

simStruct.fs =fs;
simStruct.PosUncertsigma = 0.0004^2 +.1^2 + .3^2; % seconds see EM Nosal 20
simStruct.PosUncertsigma = 0.0067^2 +.003^2 + .25^2; % seconds see EM Nosal 20- R position, ssp, arrival time
simStruct.drift = 0.25; % clock drift uncertainty
simStruct.c =1430;  % speed of sound
simStruct.maxEltTime = 30*60; % Maximum elapsed time between detections to make another encounter
simStruct.s = 3;% maximum swim speed of the animal in m/s
simStruct.RefScores = detections.Score(detectionTrimmed.Channel==ref_chan);


simStruct.dex = dexTrimmed;
[simStruct.pruned truthTrimmed]= validateDetection(TruthTrimmed, simStruct, 10);
% Three different apporaches


%1 - TDOA
simStructTDOA= simStruct;
simStructTDOA.Sim_mat = simMatTDOAonlyKAC(simStruct);

%2 - Time Clustering
simStructAcousticEncounders = simStruct;
simStructAcousticEncounders.Cluster_id = acEnc(simStructAcousticEncounders);
simStructAcousticEncounders.Sim_mat = ones(size(simStructTDOA.Sim_mat));
