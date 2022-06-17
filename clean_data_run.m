% Load data (cleaned by Pablo)

data = readtable("data_to_run.csv");

sub = table2cell(data);
subset2 = cell2mat(sub);

subset = subset2;

subset(:,1) = [];
%%
for i=1:600 
    FEATURES{i,1}=subset(:,i)';
end

F_names = data.Properties.VariableNames;

%%
data_to_paint = readtable("data_U.csv");

data_to_paint(:,602) = [];
sub0 = table2cell(data_to_paint);
subset20 = cell2mat(sub0);

%subset = subset20(:,[2 5 8 19 16 17 19 25 28 32 43 44 73 81 99 130 170 178 179 190 199 212 220 303 309 334 335 353 355 356 360 363 369 372 373 374 378 380 396 400 401 405 408 414 416 486 495 502 524 525]);

%names = data_to_paint(:,[2 5 8 19 16 17 19 25 28 32 43 44 73 81 99 130 170 178 179 190 199 212 220 303 309 334 335 353 355 356 360 363 369 372 373 374 378 380 396 400 401 405 408 414 416 486 495 502 524 525]);
%% Remove preclamsia data

subset20(104:129,:) = [];
%% SELECT 
for i=1:50
    FEATURES{i,1}=subset(:,i)';
end

F_names = names.Properties.VariableNames;

%%
subset_idclust = [subset idx_cluster];
clin = subset20(:,602:605);
final_m = [subset_idclust clin];
writematrix(final_m,'sub94_clust.csv') 
%%
writecell(F_names,'F_names.csv')
