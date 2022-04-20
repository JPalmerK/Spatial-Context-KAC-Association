
function [samps] = readSampData(file_listTable, callStartSamp, sampleStop, event_pad)
% Pull the signal and buffer associated with the detection



    sampleStart = callStartSamp-event_pad;
    sampleStop = sampleStop+event_pad;

    % If early in the file pad may overstep
    if sampleStart<1
        beginBuffpts=-sampleStart;
        sampleStart=1;
        
    else
        beginBuffpts=[]
    end
    
 
    

    % the file where the clip is
    fileIdx = find(file_listTable.SampleStrt+1 > sampleStart ,1)-1;
    
    if isempty(fileIdx)
        fileIdx  = height(file_listTable);
    end
    
    % Where in file the sample starts
    fstart = round(sampleStart-file_listTable.SampleStrt(fileIdx));
    fstop = fstart+(sampleStop-sampleStart);
    
   
    
%     if fstop>(fs*60*15 + (fs*3))
%         disp('feeeck')
%     end
    
    % samples overlap files
    if fstop>1800000
        
        % load the audio
        samps = audioread(fullfile(file_listTable.folder{fileIdx},...
            file_listTable.name{fileIdx}),...
            [fstart 1800000]);
        
        try
        samps = [samps;...
            audioread(fullfile(file_listTable.folder{fileIdx+1},...
            file_listTable.name{fileIdx+1}),...
            [1 fstop-1800000-1])];
        catch
            disp('last files')
         samps = [samps; zeros(fstop-1800000, size(samps,2)); ] ;
        
        end
            
    else
        
        samps = audioread(fullfile(file_listTable.folder{fileIdx},...
            file_listTable.name{fileIdx}),...
            [fstart fstop-1]);  
    end
    
    if beginBuffpts>0
    % Add the buffer points to the sample
    samps = [zeros(beginBuffpts, size(samps,2)); samps] ;
    end
    

end