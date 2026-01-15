function [ress_ts] = evecs2ress(data, evecs, compnr)
% reconstruct RESS component time series
 ress_ts = zeros(size(data, 2),size(data, 3));
 for ti=1:size(data,3)
     if length(size(data)) < 3
         ress_ts = zeros(size(data,2),1);
         ress_ts(:) = evecs(:,compnr)'*data; % 
     else  
        ress_ts(:,ti) = evecs(:,compnr)'*squeeze(data(:,:,ti)); % 
     end     
 end  
end