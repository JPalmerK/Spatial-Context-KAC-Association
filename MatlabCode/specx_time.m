function [corrfunc, lagvec] = specx_time(spec1, spec2, maxlag, diagnostics);


%%
%% 'specx_time'
%%
%% Kathryn A. Cortopassi, August 2004
%%
%% Perform a matrix (typically spectrogram) cross-correlation
%% over varying time lags; return the cross-correlation function
%% Cross is over time lags only, no frequency lags
%%
%% syntax:
%% -------
%%   corrfunc = specx_time(spec1, spec2, maxlag);
%%
%% input:
%% ------
%%   spec1 == spectrogram matrix
%%   spec2 == spectrogram matrix
%%   maxlag == max allowed time lag in bins
%%   diagnostics == optional flag indicating whether to plot diagnostics
%%
%% output:
%% -------
%%   corrfunc == cross-correlation function
%%   lagvec == lag vector in bins
%%


%% check that number of rows, i.e., frequency bands are equal
if size(spec1, 1) ~=  size(spec2, 1)
  error('number of rows (frequency bands) for matrices (spectrograms) are not equal');
  return;
end


%% calculate the spectrogram cross
crossmat = (spec2') * (spec1);


[m, n] = size(crossmat);

%% make sure the specified lag range does not exceed the
%% available lag range
krange = max(-maxlag, -(m-1)):min(maxlag, n-1);


%% generate the cross-correlation function for the
%% lag range
index = 0;
for k = krange
  
  index = index + 1;
  
  lagvec(index) = k;
  
  if any([m, n] == 1)
    corrfunc(index) = crossmat(k+1);
  else  
    corrfunc(index) = sum(diag(crossmat, k));
  end
  
end


if exist('diagnostics') & strcmpi(diagnostics, 'diagnostics')
  h = figure; 
  set(h, 'name', '''specx_time'' results');
  
  a(1) = subplot(3,2,1);
  imagesc(spec1); axis xy; colorbar; title('spec 1');
  a(2) = subplot(3,2,3);
  imagesc(spec2); axis xy; colorbar; title('spec 2');
  xlim1 = get(a(1), 'xlim');
  xlim2 = get(a(2), 'xlim');
  xlim = minmax([xlim1, xlim2]);
  set(a, 'xlim', xlim);
  
  subplot(3,2,[2,4]);
  imagesc(crossmat); colorbar; title('spec cross');
  
  subplot(3,2,[5,6]);
  plot(lagvec, corrfunc); title('corr func');
  
end



return;