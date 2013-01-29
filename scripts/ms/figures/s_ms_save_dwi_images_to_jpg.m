function s_ms_save_dwi_images_to_jpg

% CD to the diffusion daa directory
mypaths('diffusion',1)

% CD to the folder containing the diffusion data files
cd pestilli/20110922_1125/raw/

% Make a figure name with full path to the correct folder.
figPath = '/home/frk/Dropbox/dwi_figs/dwi_pestilli_2mm_bval2000_slice33';

%% Load a diffusion file
dwi1 = niftiRead('0005_01_DTI_2mm_150dir_2x_b2000_aligned_trilin.nii.gz');
dwi2 = niftiRead('0007_01_DTI_2mm_150dir_2x_b2000_aligned_trilin.nii.gz');
im1 = cell(size(dwi1.data,4),1); 
im2 = im1;  % Second data set
imP = im1;  % Modle prediction
imEm = im1; % Model error
imEd = im1; % Data error
imEr = im1; % RMSE ratio

%% Make an image of data set 1
% Loop over a few images and print them out.
for i_m = 1:size(dwi1.data,4)
  im1{i_m} = double(squeeze(dwi1.data(7:end-7,12:end-10,33,i_m)));
  imagesc(im1{i_m});
  colormap bone
  axis tight
  axis off
  drawnow
  
  % Save the fgure to file,
  figName = sprintf('%s_data1_dir%i',figPath,i_m);
  fprintf('Saving figure... \n%s\n',figName);
  eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(gcf),figName) );
end

%% Make an image of data set 2
% Loop over a few images and print them out.
for i_m = 1:size(dwi1.data,4)
  im2{i_m} = double(squeeze(dwi2.data(7:end-7,12:end-10,33,i_m)));
  imagesc(im2{i_m});
  colormap bone
  axis tight
  axis off
  drawnow
  
  % Save the fgure to file,
  figName = sprintf('%s_data2_dir%i',figPath,i_m);
  fprintf('Saving figure... \n%s\n',figName);
  eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(gcf),figName) );
end

%% Now simulate  the model-predicted signal
% Loop over a few images and print them out.
for i_m = 1:size(dwi1.data,4)  
  % generate and apply a low pass filter
  flt     = fspecial('gaussian',[3 11],4); 
  % Make a model prediction as the min of the two data sets
  a = ones(size(im2{i_m},1),size(im2{i_m},2),2);
  a(:,:,1)=im1{i_m};
  a(:,:,2)=im2{i_m};
  m = min(a, [], 3);
  
  %m = im1{i_m} - im2{i_m}*0.5;
  imP{i_m} = imfilter(m,flt,'same');
  % add some noise
  imP{i_m} = 0.1*median(imP{i_m}(:)).*randn(size(imP{i_m})) + imP{i_m};
  imP{i_m}( imP{i_m} < 48 ) = 1;
  imagesc(imP{i_m});
  colormap copper
  axis tight
  axis off
  drawnow
  
  % Save the fgure to file,
  figName = sprintf('%s_prediction_dir%i',figPath,i_m);
  fprintf('Saving figure... \n%s\n',figName);
  eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(gcf),figName) );
end

%% (Model accuracy) Compute the RMSE of the model prediction over the second data set (cross validation)
for i_m = 1:size(dwi1.data,4)  
  % the error image
  imEm{i_m} = sqrt( (im2{i_m} - imP{i_m}).^2 ) + 1;
  imEm{i_m}( imP{i_m} < 48 ) = max(imEm{i_m}(:));
  imagesc(imEm{i_m});
  colormap(flipud(hot))
  axis tight
  axis off
  drawnow
  
  % Save the fgure to file,
  figName = sprintf('%s_model_accuracy_dir%i',figPath,i_m);
  fprintf('Saving figure... \n%s\n',figName);
  eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(gcf),figName) );
end

%% (Data reliability) Compute the RMSE of the first data set over the second (cross validation)
for i_m = 1:size(dwi1.data,4)  
  % the error image
  imEd{i_m} = sqrt( (im2{i_m} - im1{i_m}).^2 ) + 1;
  imEd{i_m}( imP{i_m} < 48 ) = max(imEd{i_m}(:));
  imagesc(imEd{i_m});
  colormap(flipud(hot))
  axis tight
  axis off
  drawnow
  
  % Save the fgure to file,
  figName = sprintf('%s_data_reliability_dir%i',figPath,i_m);
  fprintf('Saving figure... \n%s\n',figName);
  eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(gcf),figName) );
end

%% (RMSE Ratio) Compute the RMSE of the model prediction over the second data set (cross validation)
for i_m = 1:size(dwi1.data,4)  
  % the error image
  imEr{i_m} = (imEm{i_m} ./ imEd{i_m});
  imEr{i_m}( imEr{i_m} > 4.5 ) = 2;
  imEr{i_m}( imEr{i_m} < 1 ) = .1;
  imEr{i_m}( imP{i_m} < 48 ) = max(imEr{i_m}(:));
  
  imagesc((imEr{i_m}));
  colormap(flipud(hot(250)))
  axis tight
  axis off
  drawnow
  
  % Save the fgure to file,
  figName = sprintf('%s_rmse_ratio_dir%i',figPath,i_m);
  fprintf('Saving figure... \n%s\n',figName);
  eval( sprintf('print(%s, ''-djpeg'',''-r500'', ''%s'');', num2str(gcf),figName) );
end


keyboard