
%%% Included in Workspace %%%
close all; clear all;clc
cd('D:\Spatial-Context\Scripts Sans GPL')

% Parameters
fs =2000;

% List of sound files
cwd ='D:\Data\DCLDE2013\DCLDE_2013_10Channel\Day1';
cwd = 'D:\Data\DCLDE2013\DCLDE_2013_10Channel\All';
file_list=(  dir(fullfile(cwd, '*.aif')));

% LOAD THE ARRAY STRUCTURE
load('D:\Spatial-Context\Scripts\DCLDE2013_RWDCLDE_2013_10_Chan201914101_96_array_struct_data.mat')

detectionsLocation = 'D:\Data\DCLDE2013\NN_output'
NNoutput = dir(fullfile(detectionsLocation, '*.txt'))
% Load DCLDE meta data
% Load/Create Hydrophone Structure
dclde_2013_meta = xlsread('D:/DCL2013_NEFSC_SBNMS_200903_metadata.xlsx');

detectionsOut =[];
% load the NN detectior output (set to 1 for debugging)
for ii = 1:1%length(NNoutput)

    
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
truthFiles = dir('D:\Data\DCLDE2013\LogsWithNewSNR\RavenST\CorrectChannels\*.csv')
detections = detectionsOut;
%% Hydrophone data


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




%% Load the truth data


truthStream =[]
for ii=1:1%length(truthFiles)
    
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
%% list of the soundfiles

file_listTable = struct2table(file_list)
file_listTable.MatlabDate = cell2mat(vertcat(cellfun(@(x) datenum(x(end-18:end-4), 'yyyymmdd_HHMMSS'), ...
    file_listTable.name, 'UniformOutput', false)));

% Samples into the entire stream
file_listTable.SampleStrt = fs*60*60*24*(file_listTable.MatlabDate-floor(min(file_listTable.MatlabDate)));

% Seconds inot the s
file_listTable.SecStart = (file_listTable.MatlabDate-floor(file_listTable.MatlabDate))*60*60*24;
file_listTable.SampStop =file_listTable.SampleStrt+(15*60*fs);



%get the x corr times

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


dexOrig =dex;
delaysOrig =delays;
crossScoresOrig = crossScores;
detectionsRefOrig =detectionsRef;

idxKeep = find(sum(~isnan(delays),2)>0);

dex=dex(idxKeep,:);
delays=delays(idxKeep,:);
crossScores=crossScores(idxKeep,:);
detectionsRef=detectionsRef(idxKeep,:);
refScores = detections.Score(detectionsRef.Channel==ref_chan);
%% Set up the simulation structure
arrivalArray = detectionsRef.BeginTime_s_*ones(1,10)+delays;
arrivalArray(:,10)=detectionsRef.BeginTime_s_;

simStruct = struct();
simStruct.arrivalArray = [arrivalArray(:,ref_chan) arrivalArray(:,1:end-1)];

% trim the cross scores


simStruct.parent =10;
simStruct.fs =fs;
simStruct.PosUncertsigma = 0.0004^2 +.1^2 + .3^2; % seconds see EM Nosal 20
simStruct.PosUncertsigma = 0.0067^2 +.003^2 + .25^2; % seconds see EM Nosal 20- R position, ssp, arrival time
simStruct.drift = 0.25; % clock drift uncertainty
simStruct.c =1500;  % speed of sound
simStruct.truncateKm=18; % replaced by half normal function, previously maximum detection distance
simStruct.maxEltTime = 10*60; % Maximum elapsed time between detections to make another encounter
simStruct.s = 3;% maximum swim speed of the animal in m/s
simStruct.RefScores = refScores;
%simStruct.propAmbLoc ='C:\Data\DCLDE2013_ProjAmbSurfs_KAC'
simStruct.array_struct = array_struct_data(ref_chan).array;
simStruct.dex = dex;
simStruct.TDOA_vals=delays;
simStruct.child =1:9;
simStruct.child_idx= 1:9;
% Three different apporaches

% Link the truth detections with the NN output
[simStruct.pruned truth]= validateDetection(truth, simStruct, ref_chan);


%%
% % Create the similarity matrix using pre-computed ambiguity surfaces
% simStruct = simMatMaxofProdPreComputedKAC(simStruct);
% exampSpatial.Sim_mat= round(simMatTemp,2);


%%

%1 - TDOA
simStructTDOA= simStruct;
simStructTDOA.Sim_mat = simMatTDOAonlyKAC(simStruct);



%%
%2 - Time Clustering
simStructAcousticEncounders = simStruct;
simStructAcousticEncounders.Cluster_id = acEnc(simStructAcousticEncounders);
simStructAcousticEncounders.Sim_mat = ones(size(simStructTDOA.Sim_mat));






%% Create the results field to look at precision and recall as a function
% of the max elapsed time and similarity threshold

simThreshs = [.6 .75 .98];
simThresh=linspace(.9,.998, 3);
TimeThreshs = [5 15  30];
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
%%
