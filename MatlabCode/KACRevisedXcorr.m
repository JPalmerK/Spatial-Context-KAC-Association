function [value] = KACRevisedXcorr( file_listTable,Hyddelay, noisebuffpts, ...
    ref_chan, num_chan, flo, fhi, callStartSamp, callStopSamp)

%% Input
% file_listTable- a table containing the name of each sound file, folder location,
% stream sample start of each file, stream sample stop of each file, folder
% Hyddelay- time delay (in seconds) between the primary hydrophone and all 
% secondary hydrophones
% noisebuffpts- length in samples of the noise buffer. Typically 5 seconds
% worth
% ref_chan- which is the primary channel from which to pull the initial
% detection
% num_chan- number of channels (this could be gleaned from the data
% assuming all channels you want to xcorr)
% flo- low frequency cut-off
% fhigh - high frequency cutoff
% callStartSamp - sample in the datastream of the call start
% callStopSamp - sample in the datastream of the call stop


use_ce = 1;
fs= 2000;
% get max delay
 max_delay = max(Hyddelay);
% noisePadLength =0; % seconds on either side of the even to pad
% 
% % convert noise buffer from seconds to samples
% noisebuffpts = round(noisePadLength * fs);

% to acquire the data needed for all calculations to follow
event_pad = max(noisebuffpts, max_delay);


% duration is always 2 sec for NN output
totalpts = ceil(2 * fs) + (2 * event_pad);

% KJP modification- pull the sample data (previously readsound or similar
% from xbat)
[sigdata] = readSampData(file_listTable, callStartSamp, callStopSamp, event_pad);

% DC offset per channel
sigdata = sigdata - ones(totalpts, 1) * mean(sigdata);

eventinx = (event_pad + 1):(totalpts - event_pad);
noise_start = max(1, eventinx(1) - noisebuffpts);
noise_end = min(totalpts, eventinx(end) + noisebuffpts);
noiseinx = [noise_start:(eventinx(1) - 1), (eventinx(end) + 1):noise_end];

% scaling factor, the window size equals fft size
win_func = feval(get_taper_func('Hann'), 512);
win_scale = sum(win_func.^2) / 512;


% calculate the overlap for this particular
% taper to give correct sum-square calculation
% in Parseval-sense
ovrlap_pts = round(512 * (1 - win_scale));
hop_pts =512 - ovrlap_pts;
df = fs / 512;

% calculate event and noise buffer STFFTs

% zero-pad the begin and end of the signal with ovrlap_pts/2
% number of points to eliminate end-effects coming from use
% of the overlap method for STFT generation
end_pad_pts = round(ovrlap_pts/2);
end_pad = zeros(end_pad_pts, 1);

%% get event data
sig = [end_pad; sigdata(eventinx, ref_chan); end_pad];
[ref_eventspcs] = specgram(sig, 512, fs, win_func, ovrlap_pts);


num_spc = size(ref_eventspcs, 2);

% get total number of points used for spectrogram calculation
Ntotal = ((num_spc - 1) * (512 - ovrlap_pts)) + 512;
% get total number of event data points used (i.e., excluding zero-
% padded section)
Npts = min(length(eventinx), (Ntotal - end_pad_pts));

% convert to correctly scaled mag-sqrd values in Parseval-sense
ref_eventspcs = ref_eventspcs .* conj(ref_eventspcs);
ref_eventspcs(2:end-1, :) = 2 * ref_eventspcs(2:end-1, :);
% don't need this if end zero-pad is used
ref_eventspcs(:, [1, end]) = ref_eventspcs(:, [1, end]) / (win_scale * 2);
ref_eventspcs = ref_eventspcs / 512;


% get background noise
sig = [end_pad; sigdata(noiseinx, ref_chan); end_pad];
[noisespcs] = specgram(sig,512, fs, win_func, ovrlap_pts);

% convert to correctly scaled mag-sqrd values in Parseval-sense
noisespcs = noisespcs .* conj(noisespcs);
noisespcs(2:end-1, :) = 2 * noisespcs(2:end-1, :);

% don't need this if end zero-pad is used
noisespcs(:, [1, end]) = noisespcs(:, [1, end]) / (win_scale * 2);
noisespcs = noisespcs / 512;

% get average noise spectrum
avg_noisespc = sum(noisespcs, 2) / size(noisespcs, 2);

% do some quick error checking and get frequency pointers
maxInx = size(ref_eventspcs, 1);
frange = max(1, min(maxInx, 1+round(flo/df))):max(1, min(maxInx, 1+round(fhi/df)));


% since we're only interested in the sum-squared value
% (and not display of the cleaned spectrogram),
% we can just subtract the contribution of the additive
% noise at the end
noise_contrib = num_spc * sum(avg_noisespc(frange, :));


% get sum-squared pressure
ref_ss_press = sum(sum(ref_eventspcs(frange, :)));

% correct for contribution of additive noise
ref_ss_press = ref_ss_press - noise_contrib;

% get RMS pressure
ref_rms_press = sqrt(ref_ss_press / Npts);
ref_rms_intensity = ref_rms_press^2 / (1023.6 * 1430);



% get ref event spectrogram without end padding
[ref_eventspcs] = specgram(sigdata(eventinx, ref_chan), 512, fs, win_func, ovrlap_pts);
% bandlimit
refeventspcs_blm = ref_eventspcs(frange, :);
% convert to magnitude
refeventspcs_blm = abs(refeventspcs_blm);
% convert to dB (may want to use power instead...)
refeventspcs_blm = 20 * log10(refeventspcs_blm);
% standardize for x-corr - kjp modification for normalizing across time, 
refeventspcs_blm = standardize_mat(refeventspcs_blm, 'row');

% get the length of the reference event (time wavform) in points
refevent_len = length(eventinx);

for inx = 1:num_chan
    
    % correlate event with the appropriate length of signal from cross channel
    % (why correlate across any more than I need to...)
    
    % get index for xchan buff relative the start of the extracted
    % chunk; be sure it doesn't go out of bounds (even though it shouldn't)
    xchaninx{inx} = max(1, eventinx(1) - Hyddelay(inx)):min(totalpts, eventinx(end) + Hyddelay(inx));
    
    
    % calculate the cross channel spectrogram
    xeventspc{inx} = specgram(sigdata(xchaninx{inx}, inx), 512, fs, win_func, ovrlap_pts);
    % bandlimit the cross channel spectrogram
    xeventspc_blm{inx} = xeventspc{inx}(frange, :);
    % convert to magnitude
    xeventspc_blm{inx} = abs(xeventspc_blm{inx});
    % convert to dB (may want to use power instead)
    xeventspc_blm{inx} = 20 * log10(xeventspc_blm{inx});
    % standardize the spectrogram for x-corr
    % (standardize over the entire chunk for now; i.e., don't worry yet about standardizing over
    % the reference event duration only)
    xeventspc_blm{inx} = standardize_mat(xeventspc_blm{inx}, 'row');
    
    % get correlation function over max allowed lag (i.e., the entire length of the x-channel chunk)
    maxlag = size(xeventspc_blm{inx}, 2);
    [corr_func, lags] = specx_time(xeventspc_blm{inx}, refeventspcs_blm, maxlag); %%, 'diagnostics');
    
    
    
    
    if use_ce
        
        % calculate the complex envelope of the xcorr function
        % to make peak picking easier
        complex_env = abs(hilbert(corr_func));
        
        % by convention,  specx_time(a, b) slides b against a, so we can read
        % off lags directly; get best lag
        [junk, lag(inx)] = max(complex_env);
        
    else
        
        % use straight correlation function
        % by convention,  specx_time(a, b) slides b against a, so we can read
        % off lags directly; get best lag
        [junk, lag(inx)] = max(corr_func);
        
    end
    
    % get lag in spec bins
    lag(inx) = lags(lag(inx));
    % convert to lag in time bins
    lag(inx) = lag(inx) * hop_pts;
    
    
    % get index to the cross channel event
    xev_start = (xchaninx{inx}(1) + lag(inx));
    xev_end = (xchaninx{inx}(1) + lag(inx) + refevent_len - 1);
    
    % pad each end of the x-channel event with hop_pts if possible
    % for time xcorr refinement below
    xev_start = max(1, (xev_start - hop_pts));
    xev_end = min(totalpts, (xev_end + hop_pts));
    
    xeventinx{inx} = xev_start:xev_end;
    
    
    
    
end
% code pulled from KA Cortopassi 'calc_Location_v2p2_measurement.m'
% determine the order for the FIR bandlimiting (??or, should this go right into 'bandlimitsnd'??)
data_len = size(sigdata, 1); %% length of data to be filtered
fftsz = 512; %% # of bins for parsing the frequency range
k = 5; %% # of frequency bins to use for transition band
max_ord = 500; %% arbitrary upper limit on order
c = 2; %% constant for FIR order formula M = c * fs / delta_f_trans (determined empirically

%% set the width of the transition band in Hz
delta_f_trans =  k*fs/fftsz;
if (flo ~= 0)
    delta_f_trans = min([flo, delta_f_trans]);
end
if (fhi ~= fs/2)
    delta_f_trans = min([(fs/2-fhi), delta_f_trans]);
end

% set the FIR filter order
orderVal = min([round(c*fs/delta_f_trans), floor(data_len/3), max_ord]);

% bandlimit the signal data (time waveform) for use below
% (filter the sound data in the time domain over the event frequency bounds)
[sigdata_blm] = bandlimsnd(sigdata, fs, flo, fhi, orderVal);

% get bandlimited ref event time waveform
refevent_blm = sigdata_blm(eventinx, ref_chan);

for inx = 1:num_chan
    
    % correlate ref event with cross channel event
    
    % get bandlimited cross channel event time waveform
    xevent_blm{inx} = sigdata_blm(xeventinx{inx}, inx);
    
    % refine the alignment between the ref event and cross channel event with a time xcorr
    % return the raw correlation values
    % set max allowed lag to twice hop points
    maxlag = 2*hop_pts;
    [corr_func, lags] = xcorr(xevent_blm{inx}, refevent_blm, maxlag, 'none');
    
    
    if use_ce
        
        % calculate the complex envelope of the xcorr function
        % to make peak picking easier
        complex_env = abs(hilbert(corr_func));
        
        % by convention,  xcorr(a, b) slides b against a, so we can read
        % off lags directly; get xcorr value and lag
        [xc_val(inx), lag2(inx)] = max(complex_env);
        
    else
        
        % use straight correlation function
        % by convention,  xcorr(a, b) slides b against a, so we can read
        % off lags directly; get xcorr value and lag
        [xc_val(inx), lag2(inx)] = max(corr_func);
        
    end
    
    % get any additional lag in time bins
    lag2(inx) = lags(lag2(inx)) - hop_pts;
    
    % adjust the overall lag in time bins from above
    lag(inx) = lag(inx) + lag2(inx);
    
    % use lag to find event start in seconds releative to sound stream
    % remeber first position in stream is zero (duration is always the same
    % as event duration)
    event_time(inx) = (callStartSamp-event_pad + (xchaninx{inx}(1) - 1) + lag(inx)) / fs;
    
    % also get index to the cross channel event
    xeventinx{inx} = (xchaninx{inx}(1) + lag(inx)):(xchaninx{inx}(1) + lag(inx) + refevent_len - 1);
    
end


% get the tricky little scaling factors for all channels (their xcorr peaks
% divided by the ref auto-corr peak) to use in calculating their pressure values
% from the reference event pressure value
channel_scale_W = xc_val / xc_val(ref_chan);

% Wait!  do this slightly differently, instead of using reference channel
% auto-corr peak (which hasn't been noise corrected, and which will
% show a significant contribution due to the auto-corr of the noise
% with itself) use the noise corrected reference sum-of-squares pressure;
% which should equal the auto-corr of the noise-free reference signal
channel_scale_W = xc_val / ref_ss_press;

% set the scaling of the reference channel to one, as it should be
channel_scale_W(ref_chan) = 1;

% get all channel sum-sqr and RMS pressures
all_ss_press = ref_ss_press * (channel_scale_W .^ 2);
all_rms_press = ref_rms_press * channel_scale_W;

% get all channel rms intensity (mag) in Watts / meter-cubed using density and speed of sound
all_rms_intensity = all_rms_press.^2 / (1023.6 * 1430);


% now make sure that the reference event time is really equal to the reference event time
% i.e. force a match to what we were given
event_time(ref_chan) = callStartSamp/fs;

% make a holder for any 'bad' channels
crap_chan = [];

for i = 1:num_chan
    % find the presumed events on the other channels, see if they occurred in
    % acceptable spots
    
    if lag(i) < 0
        % no part of the cross channel matched the reference event well
        % within the maximum allowed time delays (i.e., ref event fell off
        % the left edge)
        crap_chan = [crap_chan, i];
        % stick in a place holder for this bad channel and deal with it below
        xevent_blm{i} = 1;
        
    elseif (lag(i)+refevent_len) > length(xchaninx{i})
        % no part of the cross channel matched the reference event well within the
        % maximum allowed time delays (i.e., ref event fell off the right edge)
        crap_chan = [crap_chan, i];
        % stick in a place holder for crap channel and deal with it below
        xevent_blm{i} = 1;
        
    else
        % there was a good match within the allowed time delays
        % get the bandlimited cross channel event and calculate its norm
        xevent_blm{i} = norm(sigdata_blm(xeventinx{i}, i));
        
    end
    
end


% convert to an array of norms
xevent_blm = cell2mat(xevent_blm);

% calculate the norm of the reference event
refevent_blm = norm(refevent_blm);

% use these norms to calculate the normalized cross-correlation value
xc_val2 = xc_val ./ (xevent_blm * refevent_blm);

% set the norm xcorr value for crap channels to zero
xc_val2(crap_chan) = 0;

% return the raw and normalized x-corr values
value.all_pk_xcorr_raw = [xc_val]; %% raw peak cross-correlation values for all channels with reference channel
value.all_pk_xcorr_norm = [xc_val2]; %% normalized peak cross-correlation values for all channels with reference channel

%% calculate all the pairwise time lags
index = 1;
for i = 1:num_chan-1
    for j = 2:num_chan
        if (i < j)
            %pairwise_delta_t(index) = abs(event_time(i) - event_time(j));
            %% return signed values
            pairwise_delta_t(index) = event_time(i) - event_time(j);
            index = index + 1;
        end
    end
end

% return the time-delays
value.pairwise_delta_t =  pairwise_delta_t; %% all pairwise time delays between channels in the order 1x2, 1x3, ..., 1xn, 2x3, 2x4, ..., 2xn, ..., (n-1)xn
value.event_time=event_time;

end