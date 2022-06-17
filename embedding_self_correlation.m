function [self_corr,n_dims_min_out] = embedding_self_correlation(embedding1)
% This code is used to assess the number of relevant dimensions to
% consider from an MKL output space. To do so, it ranks the neighbors for
% each of the data entry (patient) and compares it to the same ranking but
% computed with one more dimension. This way, when the ranking of nearest
% neighbors does not change, we can argue that the space has stabilized.

N = length(embedding1(:,1));
n_dims_analyze = length(embedding1(1,:));
n_dims_min_out = n_dims_analyze;

self_corr = zeros(1,n_dims_analyze);
stop_aux = 0;
for count_dims = 2:n_dims_analyze
    disp(['Analyzing ',num2str(count_dims),' dimensions'])
    
    rank_comp_vec = zeros(1,N);
    k_1 = squareform(pdist(embedding1(:,1:count_dims-1))); % Distance matrix with n-1 dimensions
    k_2 = squareform(pdist(embedding1(:,1:count_dims))); % Distance matrix with n dimensions

            % For all data entries
            parfor aux_knn = 1:length(k_1)
                % Compute the ranking with n-1 dimensions
                [~,idx_1] = sort(k_1(aux_knn,:),'ascend');
                [~,r_1]=sort(idx_1);        
                % Compute the ranking with n dimensions
                [~,idx_2] = sort(k_2(aux_knn,:),'ascend');
                [~,r_2]=sort(idx_2);        
                
                A = [r_1', r_2'];           % Create a 2 column array containing the two rankings 
                B = sortrows(A,1);          % Order the rankings based on the first column
                
                % Order based Spearman/Kendall rank correlation
                rank_comp_vec(aux_knn) = corr(B(:,1),B(:,2),'type','Spearman'); 
            end
            sim_metric = mean(rank_comp_vec); % Average the correlation of all patient rankings

            %% Minimum number of dimensions
            % Set as minimum when the correlation is over 0.9
            self_corr(count_dims) = sim_metric;
            if self_corr(count_dims) > .9 && stop_aux ~= 1
                n_dims_min_out = count_dims;
                stop_aux = 1;
            end
    end
end

