%%% Included in Workspace %%%
close all; clear all;clc
cd('D:\Spatial-Context\Scripts Sans GPL')

% Parameters
fs =2000;

% List of sound files
cwd ='D:\Data\DCLDE2013\DCLDE_2013_10Channel\Day1';
cwd = 'D:\Data\DCLDE2013\DCLDE_2013_10Channel\All';
file_list=(  dir(fullfile(cwd, '*.aif')));


detectionsLocation = 'D:\Data\DCLDE2013\NN_output'
NNoutput = dir(fullfile(detectionsLocation, '*.txt'))

detectionsOut =[];
for ii = 1:4%length(NNoutput)

    
    fname = fullfile(NNoutput(ii).folder, NNoutput(ii).name)
    detections = readtable(fname);
    detections.streamSampleStart = round(detections.BeginTime_s_*fs)+(24*(ii-1)*60*60*fs)+1;
    detections.streamSampleStop = detections.streamSampleStart+(fs*2);
    detections = sortrows(detections,'streamSampleStart','ascend');
    detections.date = ones([height(detections), 1])*datenum(NNoutput(ii).name(end-11:end-4), 'yyyymmdd');
    detections.BeginTime_s_= detections.BeginTime_s_+((ii-1)*24*60*60);
    detections.EndTime_s_= detections.EndTime_s_+((ii-1)*24*60*60);
    
    chanCall=ones(1,10);
    
    
    detectionsOut = [detectionsOut; detections];
    
    
end

 detectionsOut = sortrows(detectionsOut,'streamSampleStart','ascend');
% %  
%   detectionsOut(detectionsOut.Score<0.7,:)=[];
%  

detections = detectionsOut;
%%



% Load DCLDE meta data
% Load/Create Hydrophone Structure
dclde_2013_meta = xlsread('D:/DCL2013_NEFSC_SBNMS_200903_metadata.xlsx');

% Create labels for plotting locations (undocumented function)
labels=sprintfc('%d',1:10);

plot(dclde_2013_meta(:,12), dclde_2013_meta(:,11), 'o')
text(dclde_2013_meta(:,12), dclde_2013_meta(:,11),labels,'VerticalAlignment','top','HorizontalAlignment','left')

% Convert the meta data to a structure in the format that the
% GPL/localisation code expects

hydrophone_struct= struct();
for ii=1:size(dclde_2013_meta,1)
    hydrophone_struct(ii).name = num2str(dclde_2013_meta(ii,1));
    hydrophone_struct(ii).location = dclde_2013_meta(ii,[11:12]);
    hydrophone_struct(ii).depth= abs(dclde_2013_meta(ii, 13));
    hydrophone_struct(ii).channel=ii;
    
end



%
% parm.waveform=0;
%% Create a list of the soundfiles



file_listTable = struct2table(file_list)
file_listTable.MatlabDate = cell2mat(vertcat(cellfun(@(x) datenum(x(end-18:end-4), 'yyyymmdd_HHMMSS'), ...
    file_listTable.name, 'UniformOutput', false)));

% Samples into the entire stream
file_listTable.SampleStrt = fs*60*60*24*(file_listTable.MatlabDate-floor(min(file_listTable.MatlabDate)));

% Seconds inot the s
file_listTable.SecStart = (file_listTable.MatlabDate-floor(file_listTable.MatlabDate))*60*60*24;
file_listTable.SampStop =file_listTable.SampleStrt+(15*60*fs);



%get the x corr times between hydrophone 10 (primary) and others

for ii=1:10
    % distance between the primary and secondary
    dist(ii) = vdist2(hydrophone_struct(10).location(1),...
        hydrophone_struct(10).location(2),...
        hydrophone_struct(ii).location(1),...
        hydrophone_struct(ii).location(2));
end
hydDelay = dist/1350;
% convert delays from seconds to samples
hydDelay = round(hydDelay * fs);
% get max delay
max_delay = max(hydDelay);

% Reference channel
ref_chan=10;
num_chan=10;

%% Create the associations!

delays = nan(sum((detections.Channel==10)), 10);
dex = delays;
crossScores = dex;
arrivalArray = dex;
detectionsRef =detections(detections.Channel==ref_chan,:);


% Primary detections
for ii =1:height(detectionsRef)
    
    
    [delays(ii,:), dex(ii,:), crossScores(ii,:),val] = ...
        createAssociations(detections, ii, file_listTable,hydDelay, ref_chan, fs);
    
    if ~isempty(val)
        valueout(ii)=val;
    end
    dex(ii,ref_chan)=detectionsRef.Selection(ii);
    
    ii
    
    
end

%% Remove associations where multiple primary channel detections were linked to secondary detections

%%%%%%%%%% This actually doesn't seem to help at all. 
for jj =1:9
    
    [counts, ids] = groupcounts(dex(~isnan(dex(:,jj)),jj));
    doubleBooked =ids(counts>1);
    
    % Step through the double booked values, find the highest correlationa call
    % and for all othe rprimary calls, remove the association
    
    for ii=1:length(doubleBooked)
        
        callId = doubleBooked(ii)
        primaryCalldexs = find(dex(:,jj)==callId)
        
        % get the cross correlation scores
        corrscores = vertcat(valueout(primaryCalldexs).all_pk_xcorr_raw)
        
        % index of the strongest correlation
        [~, keepIdx] =max(corrscores(:,jj))
        primaryCalldexs(keepIdx)=[]
        
        % Allow that one to stay but clear out the arrival array, dex, and
        % cross scores for the remaining double booked detections
        delays(primaryCalldexs,jj)=nan;
        dex(primaryCalldexs,jj)=nan;
        crossScores(primaryCalldexs,jj)=nan;
        
        
    end
end

%house keeping
clear counts ids doubleBooked callId primaryCalldexs corrscores keepIdx



%% Load the truth data
truthFiles = dir('D:\Data\DCLDE2013\LogsWithNewSNR\RavenST\CorrectChannels\*.csv')

truthStream =[]
for ii=1:length(truthFiles)
    
    truthtemp = readtable(fullfile(truthFiles(ii).folder, truthFiles(ii).name));
    truthtemp=truthtemp(:,1:7);
    
    
    if ii>1
        truthtemp.Selection=truthtemp.Selection+height(truthStream)
    end
    
    truthtemp.mstart = datenum(truthFiles(ii).name(7:14), 'yyyymmdd')+truthtemp.BeginTime_s_/86400;
    truthtemp.mend = datenum(truthFiles(ii).name(7:14), 'yyyymmdd')+truthtemp.EndTime_s_/86400;
    truthtemp.mid = (truthtemp.mend+truthtemp.mstart)/2;
    truthtemp.Sec = (truthtemp.mstart-min(truthtemp.mstart))*60*60*24+(ii-1)*24*60*60;
    truthtemp.DetectorScore = zeros([height(truthtemp) 1 ]);
    truthtemp.BeginTime_s_= truthtemp.BeginTime_s_+((ii-1)*24*60*60);
    truthtemp.EndTime_s_= truthtemp.EndTime_s_+((ii-1)*24*60*60);
    truthStream=[truthStream; truthtemp];
    clear truthtemp
    ii
end

truth = truthStream;
clear truthStream truthtemp truthFiles


%% Remove detections that were not associated at least one secondary hydrophone
% delaysOrig = delays;
% dexOrig= dex;
% crossScoresOrig = crossScores;
% 
% 
% % Which detections are not catured at least one secondary channel
% delaysDitch = sum(~isnan(delays),2)==0
% crossScores(delaysDitch,:)=[];
% delays(delaysDitch,:)=[];
% dex(delaysDitch,:)=[];
% detectionsRef(delaysDitch,:)=[];
dexOrig =dex;
delaysOrig =delays;
crossScoresOrig = crossScores;
detectionsRefOrig =detectionsRef;
idxKeep = find(sum(~isnan(delays),2)>0);

% Which detections to remove from the truth
dexRem=dexOrig(find(sum(~isnan(delaysOrig),2)==0),10)

dex=dex(idxKeep,:);
delays=delays(idxKeep,:);
crossScores=crossScores(idxKeep,:);
detectionsRef=detectionsRef(idxKeep,:)


%% Set up the simulation structure
arrivalArray = detectionsRef.BeginTime_s_*ones(1,10)+delays;
arrivalArray(:,10)=detectionsRef.BeginTime_s_;

simStruct = struct();
simStruct.arrivalArray = [arrivalArray(:,ref_chan) arrivalArray(:,1:end-1)];
simStruct.TDOA_vals = delays;


simStruct.fs =fs;
simStruct.PosUncertsigma = 0.0004^2 +.1^2 + .3^2; % seconds see EM Nosal 20
simStruct.PosUncertsigma = 0.0067^2 +.003^2 + .25^2; % seconds see EM Nosal 20- R position, ssp, arrival time
simStruct.drift = 0.25; % clock drift uncertainty
simStruct.c =1500;  % speed of sound
simStruct.truncateKm=18; % replaced by half normal function, previously maximum detection distance
simStruct.maxEltTime = 30*60; % Maximum elapsed time between detections to make another encounter
simStruct.s = 3;% maximum swim speed of the animal in m/s

refChanScores =detections.Score(detections.Channel==ref_chan);
refScores=refChanScores(idxKeep)
simStruct.RefScores = refScores;


simStruct.dex = dex;
[simStruct.pruned truth]= validateDetection(truth, simStruct, ref_chan);
% Three different apporaches

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean out the truth data where detections have been removed
[tf,idx] = intersect(dexRem, truth.Dex);
truth(idx,:)=[];
[simStruct.pruned truth]= validateDetection(truth, simStruct, ref_chan);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%1 - TDOA
simStructTDOA= simStruct;
simStructTDOA.Sim_mat = simMatTDOAonlyKAC(simStruct);

%2 - Time Clustering
simStructAcousticEncounders = simStruct;
simStructAcousticEncounders.Cluster_id = acEnc(simStructAcousticEncounders);
simStructAcousticEncounders.Sim_mat = ones(size(simStructTDOA.Sim_mat));

%% dBscan clustering
clear clusterer
data = [detectionsRef.BeginTime_s_, fliplr(delays(:,7:9))]%detections.Score(detections.Channel==ref_chan)]
% data(isnan(data))= i;
% data= data(~isnan(data(:,2)),:)
%data(:,2)=[0;diff(data(:,2))] 
% data(:,1)=[1:size(data,1)];

clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(2,10,1,2,.9,1)

dataOut =zeros(size(delays,1),9)

for ii=1:9
    
    data = [detectionsRef.BeginTime_s_, fliplr(delays(:,ii))]%detections.Score(detections.Channel==ref_chan)]
    idxKeep = ~isnan(data(:,2));
    data= data(~isnan(data(:,2)),:);
    try
    clusterer = HDBSCAN(data)
    clusterer.fit_model
    clusterer.run_hdbscan(2,10,1,2,.9,1)
    dataOut(idxKeep,ii)=single(clusterer.labels)+max(dataOut(:));
    catch
        disp('oh well')
    end
    
    
end

% remove channel 4 - also not complete shit
dataOut(:,4)=[]
dataOut(dataOut(:)==0)=max(dataOut(:)):max(dataOut(:))+length(isnan(dataOut(:)))-1
datOut =[detectionsRef.BeginTime_s_, dataOut]
clusterer = HDBSCAN(datOut)
clusterer.fit_model
clusterer.run_hdbscan(2,10,1,2,.9,1)

% Try bool for present/absent - meets baseline higher (.4)
data = [detectionsRef.BeginTime_s_, ~isnan(delays(:,[1:3,5:9]))]
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(2,10,1,2,.9,1)

% keep the delays, skip channel 4- not bad
data = [detectionsRef.BeginTime_s_, fliplr(delays(:,[1:3,5:9]))]
data(isnan(data))= i;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(2,10,1,2,.9,1)


% just the local channels- no
data = [detectionsRef.BeginTime_s_, fliplr(delays(:,[6:9]))]
data(isnan(data))= i;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(2,10,1,2,.9,1)


% keep the delays, skip channel 4 and concatenated silliness meh
data = [detectionsRef.BeginTime_s_, (delays(:,[1:3,5:9])),datOut]
data(isnan(data))= i;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(2,10,1,2,.9,1)

% keep the delays, sset null value to 0
data = [detectionsRef.BeginTime_s_, fliplr(delays(:,[1:3,5:9]))]
data(isnan(data))= 0;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(2,10,1,2,.9,1)


%Create bool for which instruments detected- nope
customFx = @(aa) bin2dec(sprintf('%d', aa));
aa =~isnan(delays(:,[1:3,5:9]))
C = num2cell(aa', 1);             %# Collect the columns into cells
columnSums = cellfun(customFx, C)';  %# A 1-by-10 vector of sums for each cell
data = [detectionsRef.BeginTime_s_, fliplr(delays(:,[1:3,5:9])),columnSums]
data(isnan(data))= i;


% Skip the time element
data = [ fliplr(delays(:,[1:3,5:9]))]
data(isnan(data))= i;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(2,3,1,2,.9,1)

% Skip the time element fuck with some parameters
data = [ fliplr(delays(:,[1:3,5:9]))]
data(isnan(data))= i;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(2,3,1,2,.2,1)

% Skip the time element fuck with some parameters
data = [[0; diff(detectionsRef.BeginTime_s_)],fliplr(delays(:,[1:3,5:9]))]
data(isnan(data))= i;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(3,5,1,2,.2,1)


% Get diff time of the acoustic encounters
% Initialzie some values
timeDiffs = [0; diff(arrivalArray(:,10))]; % Time diffs for acoustic encs.

% Indicies fo the start of each acoustic encounter
acousticEncounterBreaks = [1; (find(timeDiffs>30)); height(dataTable)];

% Initialize a datatable for indexing around and preserving the origional
% information
dataTable = table();

dataTable.timeDiffs = timeDiffs;
dataTable.AcousticEncounters = zeros(size(timeDiffs));% this is just zero
dataTable.encounterTimediff = dataTable.AcousticEncounters; % this is just zero
deltaTDOAprobs = zeros(size(delays))

% Demarkate the acoustic encounters in the datatable
for ii=1:length(acousticEncounterBreaks)-1
    dataTable.AcousticEncounters(acousticEncounterBreaks(ii):end) =ii;
    dataTable.encounterTimediff(acousticEncounterBreaks(ii):end) = ...
        arrivalArray(acousticEncounterBreaks(ii):end,10)-...
        arrivalArray(acousticEncounterBreaks(ii),10);

    
    
end

data = [(dataTable.encounterTimediff) fliplr(delays(:,[1:3,5:9]))]
data(isnan(data))= i;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(3,5,1,2,.2,1)


% Justnormailize the arrival time
data = [normalize(detectionsRef.BeginTime_s_), fliplr(delays(:,[1:3,5:9]))]
data(isnan(data))= 0;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(2,10,1,2,.9,1)

%Create bool for which instruments detected- ok
dataoutTemp = dataOut+1
data = [detectionsRef.BeginTime_s_, fliplr(delays(:,[1:3,5:9])),normalize(prod(dataoutTemp,2))]
data(isnan(data))= i;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(2,10,1,2,.9,1)

%Create bool for which instruments detected- nope
dataoutTemp = dataOut+1
data = [normalize(detectionsRef.BeginTime_s_), fliplr(delays(:,[1:3,5:9])),normalize(prod(dataoutTemp,2))]
data(isnan(data))= i;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(2,10,1,2,.9,1)


%Create bool for which instruments detected ok but worse than above
dataoutTemp = dataOut+1
data = [detectionsRef.BeginTime_s_, fliplr(delays(:,[1:3,5:9])),normalize(prod(dataoutTemp,2))]
data(isnan(data))= i;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(2,5,1,2,.2,1)

%Create bool for which instruments detected- ok
dataoutTemp = dataOut+1
data = [detectionsRef.BeginTime_s_, fliplr(delays(:,[1:3,5:9])),normalize(prod(dataoutTemp,2))]
data(isnan(data))= i;
clusterer = HDBSCAN(data)
clusterer.fit_model
clusterer.run_hdbscan(3,10,1,2,.9,1)



%% dbscan chains
labels=single(clusterer.labels)
max(labels)
idxZeros = find(labels == 0)
labels(idxZeros) = 1+[max(labels):max(labels)+length(idxZeros)-1]
simStructTDOA.Cluster_id= labels


% Evaluate
[ClusteredTDOA, ~,~] = PrecisionRecallKAC(simStructTDOA, truth);

%% Create the results field to look at precision and recall as a function
% of the max elapsed time and similarity threshold

figure;
simThreshs = [ .85 .98];%linspace(.95,.998, 3);
TimeThreshs = [5 10 15];

 simThreshs = [.9];%linspace(.95,.998, 3);
% TimeThreshs = [30 ];
%%%%%%%%%%%%%%%%%% FUCING AROUND
%%simStructTDOA.arrivalArray(:,1)=1:size(simStructTDOA.arrivalArray,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Markers = {'-o','-.','-+','-*', '->', '-o','-.','-+','-*', '->'};
sizeVal = [2, 6,4,2,2,2, 6,4,2,2,2, 6,4,2,2];
markercounter=1;
for ii =1:length(simThreshs)
    
    simThresh = simThreshs(ii);
    simStructTDOA.cutoff = simThresh;
    
    for jj =1:length(TimeThreshs)
        
        
        simStructAcousticEncounders.cutoff=.5;
        disp('next simThresh')
        
        % Create precision recall curves for each similarity and time threshold
        simStructTDOA.maxEltTime = TimeThreshs(jj);
        simStructAcousticEncounders.maxEltTime = TimeThreshs(jj);
        
        % Same similarity threshold each time, only run if the time changes
        if ii==1
            
            % Baseline clustering, time only
            simStructAcousticEncounders.chains =...
                updateChainsEncounterFirstKAC(simStructAcousticEncounders);
            simStructAcousticEncounders.Cluster_id=...
                updateClusterID(simStructAcousticEncounders);
            
            [ClusteredAcousticEncouvers, ~,~] = ...
                PrecisionRecallKAC(simStructAcousticEncounders, truth);
            
            
            Results.AcousticEncounters(ii,jj) = ClusteredAcousticEncouvers;
            Results.AcousticEncountersClusters(ii,jj) = ...
                {simStructAcousticEncounders.Cluster_id};
            
             plot(Results.AcousticEncounters(ii,jj).Recall,...
                Results.AcousticEncounters(ii,jj).Precision,...
                Markers{markercounter},'color',[0/255 158/255 115/255],...
                'MarkerSize',sizeVal(markercounter))
            title('Temporal Clustering')
            hold on
            
            
            
        end
        
        %Run the clustering TDOA
        simStructTDOA.chains = updateChainsEncounterFirstKAC(simStructTDOA);
        simStructTDOA.Cluster_id= updateClusterID(simStructTDOA);
        
        % TDOA precision/recall
        [ClusteredTDOA, ~,~] = PrecisionRecallKAC(simStructTDOA, truth);
        Results.TDOA(ii,jj) = ClusteredTDOA;
        Results.TDOAClusters(ii,jj) = {simStructTDOA.Cluster_id};
        
        hold on
        plot(Results.TDOA(ii,jj).Recall,...
            Results.TDOA(ii,jj).Precision, ...
            Markers{markercounter},'color',[230/255 159/255 0/255],...
            'MarkerSize',sizeVal(markercounter))
        markercounter=markercounter+1
        
%         
%         
%         hold on
%         plot(Results.Baseline(1,1).Recall,...
%             Results.Baseline(1,1).Precision,...
%             'k-','Linewidth',2)
%         xlabel('Recall'); ylabel('Precision')
        
        
    end
    
    %cluster in time
    
end

% Do the detector onnly
[~, detector,~] = PrecisionRecallKAC(simStructTDOA, truth);
Results.Baseline(1,1) = detector;

        hold on
        plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision,...
            'k-','Linewidth',2)
        xlabel('Recall'); ylabel('Precision')
legend({'Temporal Clustering', 'TDOA Method', 'Detector' })
%%

xx = Results.Baseline(1,1).Recall;
yy = Results.Baseline


