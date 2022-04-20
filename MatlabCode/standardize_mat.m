
function matrix1 = standardize_mat(matrix1, normtype);

%%
%% 'standardize_mat'
%%
%% Kathryn A. Cortopassi, August 2004- KJP modified
%%
%% Standardize the matrix as specified
%%
%% syntax:
%% -------
%%   matrix1 = standardize_mat(matrix1, normtype);
%%
%% input:
%% ------
%%   matrix1 == a matrix of values (typically a spectrogram)
%%   normtyp == standardization method to use, either
%%              'colnorm'; 'column' for column-wise normalization
%%              'matnorm'; 'matrix' for matrix-based normalization
%%
%% output:
%% -------
%%   matrix1 == standardized matrix of values
%%



if sum(strcmpi(normtype, {'colnorm'; 'column'}))
  %% prepare the matrix by processing each column (in a spectrogram, each time slice, or FFT)
  %% by subtracting it's mean and dividing by it's standard deviation
  
  N = size(matrix1, 1);
  onesCol = ones(N, 1);
  matrix1 = ( matrix1 - (onesCol * mean(matrix1)) ) ./ ( sqrt(N-1) * (onesCol * std(matrix1)) );

  

  
  
elseif sum(strcmpi(normtype, {'matnorm'; 'matrix'}))
  %% prepare the matrix by subtracting its means and dividing by its
  %% standard deviation
  
  matrix1 = ( matrix1 - mean(matrix1(:)) ) ./ ( sqrt( length(matrix1(:))-1 ) * std(matrix1(:)) );
elseif sum(strcmp(normtype, {'Rownorm'; 'row'}))
    
      N = size(matrix1, 2);
  onesCol = ones(N,1);
  matrix1 = ( matrix1' - (onesCol * mean(matrix1,2)') ) ./ (  sqrt(N-1) * (onesCol * std(matrix1,[],2)') );
  matrix1=matrix1';
  

  
    
else
  
  error(sprintf('Normalization flag ''%s'' not recognized', normtype));
  return;
  
end


%% subtracting the mean biases the matrix or column around zero and provides
%% for both positive and negative correlation values
%%
%% dividing by the standard deviation times the square-root of NP-1 is the
%% same as dividing by the norm of the mean-adjusted matrix, 
%% this provides for correlation values between -1 and +1
%%
%% dividing by the standard deviation times the square-root of N-1 is the
%% same as dividing by the norm of the mean-adjusted column (time slice), 


return;