function [snd] = bandlimsnd(snd, fs, lo, hi, ord);

%**
%** "bandlimsnd.m"
%**
%** Kathryn A. Cortopassi, 2001-2004
%**
%** Function to bandlimit a sound waveform or cell array of sound waveforms
%**
%** Syntax: [snd] = bandlimsnd(snd, fs, lo, hi, ord);
%**
%** 'snd' is a single sound waveform or cell array of sound waveforms
%** 'fs' is the sampling rate of the sound waveforms, or a cell array of
%** sampling rates as appropriate (even if 'snd' is a cell array,
%** 'fs' can still be a scalar)
%** 'lo' is the low frequency corner in Hz
%** 'hi' is the high frequency corner in Hz
%** 'ord' is the filter order
%**
%** 'snd' returns the bandlimited sound waveform
%**

%**
%** by Kathryn A. Cortopassi
%** 12-November-2004
%** add 'ord' June 2005
%**



%% set flag to show internal diagnostic plots
plot_internal_diag = 0;


% check for correct number of inputs
argLo = 5; argHi = inf;
error(nargchk(argLo, argHi, nargin));


if ~iscell(fs) & length(fs) == 1
  wlo = lo/(fs/2);
  whi = hi/(fs/2);
end


if iscell(snd)

  for i = 1:length(snd)
    %% bandpass filter the sound waveforms based on lo and hi

    %% first get the FIR coefficients
    %% get coefficients for an order 'ord' bandpass filter

    if iscell(fs)
      wlo = lo/(fs{i}/2);
      whi = hi/(fs{i}/2);
    elseif length(fs) > 1
      wlo = lo/(fs(i)/2);
      whi = hi/(fs(i)/2);
    end

    if (wlo == 0) & (0 < whi & whi < 1)
      %% generate a lowpass filter
      b = fir1(ord, whi);
    elseif (0 < wlo & wlo < 1) & (whi == 1)
      %% generate a highpass filter
      b = fir1(ord, wlo, 'high');
    elseif (0 < wlo & wlo < 1) & (0 < whi & whi < 1)
      %% generate a bandpass filter
      b = fir1(ord, [wlo, whi]);
    elseif (wlo == 0) & (whi == 1)
      %% no filtering needed, drop through loop
      continue;
    elseif wlo == whi
      %% no filtering possible, drop through loop
      continue;
    end

    %% then do zero phase filtering; because of forward and reverse
    %% filtering, actual order will be 2*ord
    snd{i} = filtfilt(b, 1, snd{i});

    if plot_internal_diag
      fvtool(b, 1);
    end

  end


else
  %% bandpass filter the sound waveform based on lo and hi

  %% first get the FIR coefficients
  %% get coefficients for an order 'ord' bandpass filter
  if (wlo == 0) & (0 < whi & whi < 1)
    %% generate a lowpass filter
    b = fir1(ord, whi);
  elseif (0 < wlo & wlo < 1) & (whi == 1)
    %% generate a highpass filter
    b = fir1(ord, wlo, 'high');
  elseif (0 < wlo & wlo < 1) & (0 < whi & whi < 1)
    %% generate a bandpass filter
    b = fir1(ord, [wlo, whi]);
  elseif (wlo == 0) & (whi == 1)
    %% no filtering needed, drop out
    return;
  elseif wlo == whi
    %% no filterign possible, drop out
    return;
  end

  %% then do zero phase filtering; because of forward and reverse
  %% filtering, actual order will be 2*ord
  snd = filtfilt(b, 1, snd);

  if plot_internal_diag
    fvtool(b, 1);
  end

end


% end function
return;