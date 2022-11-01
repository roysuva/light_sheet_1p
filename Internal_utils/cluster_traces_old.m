function clust = cluster_traces(resp, neuronObj, varargin)
 
 
    p = inputParser;
    addParameter(p,'cluster_method','Hierarchical', @ischar);
    addParameter(p,'numclust',10, @isnumeric); 
    parse(p,varargin{:});
    params = p.Results; 
    
    
    
    % Normalize traces 
    %ind_all = 1:numel(resp.cell); 
    %outlier_idx = extract_delta(neuronObj, ind_all); 
    %ind_trim = setxor(ind_all, outlier_idx);  
    [traces1_y,traces1_y_norm] = deal(zeros(numel(resp.tracefilt_psth(:,1)),numel(resp.cell))); 
    for k=1:numel(resp.cell)
        %traces1_y(:,k) = resp.tracefilt_psth(:,ind_trim(k));
        traces1_y(:,k) = resp.trace_psth(:,k);
        traces1_y_norm(:,k) = traces1_y(:,k)./norm(traces1_y(:,k)); 
    end
    traces1_t = resp.t_axis_trial(1:size(traces1_y,1)); 
    
    
    switch params.cluster_method
        case 'Hierarchical'
            
            % Cluster data 
            mxcl = params.numclust; % max number of clusters
%             clusT = clusterdata(traces1_y_norm','linkage','complete',...
%                 'Distance','correlation','maxclust',mxcl);
            D = pdist(traces1_y_norm','seuclidean');
            Z = linkage(D,'complete');
            idx = cluster(Z,'Maxclust',mxcl);
            leafOrder = optimalleaforder(Z,D);
            
            % Generate figure 
            hf = figure; set(hf,'position',[8  140  1391  657]);
            nr = ceil(mxcl/5); nc = ceil(mxcl/nr);  
            for k=1:mxcl
                matrixmat = traces1_y_norm(:,find(idx==k))';
                subplot(nr,nc,k);
                imagesc(matrixmat); axis square; 
                title(sprintf('Cluster %d',k));
            end
            sgtitle(sprintf('Method: %s','Hierarchical Agglomerative'));  
            
            
            % Return values 
            clust.type = params.cluster_method; 
            clust.numclust = params.numclust;
            clust.cell_ind = 1:numel(resp.cell); 
            clust.clust_id = idx; 
            clust.D = D; 
            clust.Z = Z; 
            clust.leafOrder = leafOrder; 
            
            
            
        case 'GMM'
            
%           % PCA             
%             opt = optimset('MaxIter',1000,'MaxFunEvals',1000,'Display','iter');
%             [coeff, score, latent] = pca(traces1_y_norm,'Algorithm','svd','Centered',true,...
%                 'NumComponents',3,'Options',opt);
            
            % Use SVD instead of PCA 
            [U,S,V] = svd(traces1_y_norm);
            
            % Get projections
            nprincomp = 5;
            ncells = size(traces1_y_norm,2);
            proj = zeros(ncells,nprincomp);
            for j=1:nprincomp
                for k=1:ncells
                    proj(k,j) = dot(U(:,j),traces1_y_norm(:,k));
                end
            end
            
            % Fit GMM
            options = statset('Display','final','UseParallel',true,'MaxIter',3000);
            ngmmodel = params.numclust;
            gmmodel = fitgmdist(proj,ngmmodel,'CovarianceType','diagonal','RegularizationValue',0.001,'Options',options);
            gmmclust = cluster(gmmodel, proj); % cluster index
            clustname = cell(1,ngmmodel);
            for ng=1:ngmmodel
                clustname{ng} = sprintf('gmm clust %d',ng);
            end
            [sortedClust,sortedI] = sort(gmmclust,'ascend');
            
            % Generate figure 
            mxcl = params.numclust; 
            hf = figure; set(hf,'position',[8  140  1391  657]);
            for k=1:mxcl
                matrixmat = traces1_y_norm(:,find(gmmclust==k))';
                subplot(2,5,k);
                imagesc(matrixmat); axis square
                title(sprintf('Cluster %d',k));
            end
            sgtitle(sprintf('Method: %s','GMM'));
            
            % Return values 
            clust.type = params.cluster_method; 
            clust.numclust = params.numclust;
            clust.cell_ind = 1:numel(resp.cell); 
            clust.clust_id = gmmclust; 
            clust.gmmodel = gmmodel; 
            clust.svd_numcomp = nprincomp; 
            clust.svd_U = U; 
            clust.svd_S = S; 
            clust.svd_V = V; 
            clust.proj = proj; 
 
 
        case 'Spectral'
            
            mxcl = params.numclust;
            [idx,V,D] = spectralcluster(traces1_y_norm',mxcl,'Distance','correlation','SimilarityGraph','knn',...
                'ClusterMethod','kmeans','LaplacianNormalization','symmetric');
            
            hf = figure; set(hf,'position',[8  140  1391  657]);
            for k=1:mxcl
                matrixmat = traces1_y_norm(:,find(idx==k))';
                subplot(2,5,k);
                imagesc(matrixmat); axis square
                title(sprintf('Cluster %d',k));
            end
            sgtitle(sprintf('Method: %s','Spectral'));  
            
            clust.type = params.cluster_method; 
            clust.numclust = params.numclust;
            clust.cell_ind = 1:numel(resp.cell); 
            clust.clust_id = idx; 
 
    end
    
end
