function clust = cluster_traces(resp, neuronObj, varargin)
 
 
    p = inputParser;
    addParameter(p,'cluster_method','Hierarchical', @ischar);
    addParameter(p,'numclust',10, @isnumeric); 
    addParameter(p,'clustrng_for_BIC',[], @isnumeric); 
    addParameter(p,'which_trials',1, @(x)(isnumeric(x) || ischar(x))); 
    addParameter(p,'show_outlier_plot',true, @islogical); 
    parse(p,varargin{:});
    params = p.Results; 
    
    if isempty(params.clustrng_for_BIC)
        params.clustrng_for_BIC = [1 params.numclust+5]; 
    end
   
   
   
    % Normalize traces 
    if ischar(params.which_trials) && strcmpi(params.which_trials,'all')
        [traces1_y,traces1_y_norm] = deal(zeros(numel(resp.tracefilt_psth(:,1)),numel(resp.cell)));
        for k=1:numel(resp.cell)
            %traces1_y(:,k) = resp.tracefilt_psth(:,ind_trim(k));
            traces1_y(:,k) = resp.trace_psth(:,k);
            %traces1_y(:,k) = resp.cell(k).raster(1,:)';
            traces1_y_norm(:,k) = traces1_y(:,k)./norm(traces1_y(:,k));
            %traces1_y_norm(:,k) = traces1_y(:,k)./max(traces1_y(:,k));
        end
        traces1_t = resp.t_axis_trial(1:size(traces1_y,1));
    elseif params.which_trials==1
        [traces1_y,traces1_y_norm] = deal(zeros(numel(resp.tracefilt_psth(:,1)),numel(resp.cell)));
        for k=1:numel(resp.cell)
            %traces1_y(:,k) = resp.tracefilt_psth(:,ind_trim(k));
            %traces1_y(:,k) = resp.trace_psth(:,k);
            traces1_y(:,k) = resp.cell(k).raster(1,:)';
            traces1_y_norm(:,k) = traces1_y(:,k)./norm(traces1_y(:,k));
            %traces1_y_norm(:,k) = traces1_y(:,k)./max(traces1_y(:,k));
        end
        traces1_t = resp.t_axis_trial(1:size(traces1_y,1));
    end

    % Remove outliers 
    outliers = extract_outliers('traces',traces1_y_norm,'show_plot',params.show_outlier_plot);
    non_outliers = setxor(1:size(traces1_y_norm,2), outliers); 
    traces1_y_norm_trunc = traces1_y_norm(:,non_outliers); 
    
   
    
    switch params.cluster_method
        case 'Hierarchical'
            
            % Cluster data 
            mxcl = params.numclust; % max number of clusters
            Z = linkage(traces1_y_norm_trunc','ward','minkowski');     
            idx = cluster(Z,'Maxclust',mxcl);


%             idx = clusterdata(traces1_y_norm_trunc','linkage','complete',...
%                 'Distance','correlation','maxclust',mxcl);
%             D = pdist(traces1_y_norm_trunc,'correlation'); %
%             Z = linkage(D,'complete','correlation');
%             idx = cluster(Z,'Maxclust',mxcl);
%             leafOrder = optimalleaforder(Z,D);
            
            
            % Get silhouette values for cluster optimality 
            figure; 
            [s,h] = silhouette(traces1_y_norm_trunc', idx,'correlation');
            silhouette_stat.single_run.svals = s; 
            silhouette_stat.single_run.clust_range = 1:params.numclust; 
            test_clust_range = 2:15; 
            
            silhouette_stat.niter = 50; 
            evals = cell(1,silhouette_stat.niter);
            myfun = @(X,K) (cluster(linkage(X,'ward','minkowski'),'Maxclust',k)); 
            parfor_progress(silhouette_stat.niter); 
            parfor niter=1:silhouette_stat.niter    
                evaluation = evalclusters(traces1_y_norm_trunc',myfun,'Silhouette',...
                    'Distance','correlation','klist',test_clust_range); 
                evals{niter} = evaluation; 
                parfor_progress; 
            end
            parfor_progress(0); 
            silhouette_stat.test_clust_range = test_clust_range; 
            silhouette_stat.evals = evals; 




            % Generate figure 
            hf = figure; set(hf,'position',[8  140  1391  657]);
            nr = ceil(mxcl/5); nc = ceil(mxcl/nr);  
            for k=1:mxcl
                matrixmat = traces1_y_norm_trunc(:,find(idx==k))';
                subplot(nr,nc,k);
                imagesc(matrixmat); axis square;  
                title(sprintf('Cluster %d',k));
            end
            if ~isempty(which('sgtitle')) % accomodating older Matlab versions 
                sgtitle(sprintf('Method: %s','Hierarchical Agglomerative'));  
            end


            % Reassign cluster id, leafOrder 
            [reassign_clustid, reassign_lO] = deal(NaN(size(traces1_y_norm,2),1)); 
            reassign_clustid(non_outliers) = idx;  
            %reassign_lO(non_outliers) = leafOrder; 
            
            % Return values 
            clust.type = params.cluster_method; 
            clust.numclust = params.numclust;
            clust.cell_ind = 1:numel(resp.cell); 
            %clust.leafOrder = reassign_lO; 
            
            clust.type = params.cluster_method; 
            clust.numclust = params.numclust;
            clust.cell_ind = 1:numel(resp.cell); 
            clust.clust_id = reassign_clustid; 
            clust.non_outliers = non_outliers; 
            clust.outliers = outliers; 
            if exist('silhouette_stat','var')
                clust.silhouette_stat = silhouette_stat; 
            end

            
            
        case 'GMM'
            
%           % PCA             
%             opt = optimset('MaxIter',1000,'MaxFunEvals',1000,'Display','iter');
%             [coeff, score, latent] = pca(traces1_y_norm_trunc,'Algorithm','svd','Centered',true,...
%                 'NumComponents',3,'Options',opt);
            
            % Use SVD instead of PCA 
            [U,S,V] = svd(traces1_y_norm_trunc);
            nl = 1; eigs = diag(S);
            while sum(eigs(1:nl))/sum(eigs)<0.4 
                nl = nl+1; 
            end
            nprincomp = nl-1; 
            
            % Get projections
            ncells = size(traces1_y_norm_trunc,2);
            proj = zeros(ncells,nprincomp);
            for j=1:nprincomp
                for k=1:ncells
                    proj(k,j) = dot(U(:,j),traces1_y_norm_trunc(:,k));
                end
            end
            
            % Fit GMM
            options = statset('Display','final','UseParallel',true,'MaxIter',10000);
            ngmmodel = params.numclust;
            gmmodel = fitgmdist(proj,ngmmodel,'CovarianceType','diagonal','RegularizationValue',0.01,'Options',options);
            gmmclust = cluster(gmmodel, proj); % cluster index
            clustname = cell(1,ngmmodel);
            for ng=1:ngmmodel
                clustname{ng} = sprintf('gmm clust %d',ng);
            end
            [sortedClust,sortedI] = sort(gmmclust,'ascend');

            % Get silhouette values for cluster optimality 
            figure; 
            [s,h] = silhouette(proj, gmmclust);
            silhouette_stat.single_run.svals = s; 
            silhouette_stat.single_run.clust_range = 1:params.numclust; 
            test_clust_range = 2:15; 
            
            silhouette_stat.niter = 50; 
            evals = cell(1,silhouette_stat.niter); 
            parfor_progress(silhouette_stat.niter); 
            for niter=1:silhouette_stat.niter
                myfun = @(X,K) (cluster(fitgmdist(X,K,'CovarianceType','diagonal','RegularizationValue',0.01,'Options',options), X));
                evaluation = evalclusters(proj,myfun,'Silhouette','Distance','correlation','klist',test_clust_range); 
                evals{niter} = evaluation; 
                parfor_progress; 
            end
            parfor_progress(0); 
            silhouette_stat.test_clust_range = test_clust_range; 
            silhouette_stat.evals = evals; 
                

            % Get BIC (Bayesian Information Criterion for model fit) 
            bicrepeats = 50; 
            bic = zeros(numel(params.clustrng_for_BIC(1):params.clustrng_for_BIC(2)),bicrepeats); 
            for ki=1:bicrepeats % repeats 
                for bc=params.clustrng_for_BIC(1):params.clustrng_for_BIC(2) 
                    gmmodel_ = fitgmdist(proj,bc,'CovarianceType','full','RegularizationValue',0.01,...
                        'SharedCovariance',true,'Options',options);
                    bic(bc,ki) = gmmodel_.BIC; 
                end
            end
            bic_cluster_rng = params.clustrng_for_BIC(1):params.clustrng_for_BIC(2); 




            % Generate figure 
            mxcl = params.numclust; 
            hf = figure; set(hf,'position',[8  140  1391  657]);
            for k=1:mxcl
                matrixmat = traces1_y_norm_trunc(:,find(gmmclust==k))';
                subplot(2,5,k);
                imagesc(matrixmat); axis square;  
                title(sprintf('Cluster %d',k));
            end
            if ~isempty(which('sgtitle')) % accomodating older Matlab versions 
                sgtitle(sprintf('Method: %s','GMM'));
            end
           
            
            % Reassign cluster id, projection values 
            reassign_clustid = NaN(size(traces1_y_norm,2),1); 
            reassign_clustid(non_outliers) = gmmclust;  
            proj_full = NaN(size(traces1_y_norm,2),nprincomp);  
            proj_full(non_outliers,:) = proj; 

            % Return values 
            clust.type = params.cluster_method; 
            clust.numclust = params.numclust;
            clust.cell_ind = 1:numel(resp.cell); 
            clust.clust_id = reassign_clustid; 
            clust.gmmodel = gmmodel; 
            clust.svd_numcomp = nprincomp; 
            clust.svd_U = U; 
            clust.svd_S = S; 
            clust.svd_V = V; 
            clust.proj = proj; 
            clust.proj_full = proj_full; 
            clust.non_outliers = non_outliers; 
            clust.outliers = outliers; 
            clust.silhouette_stat = silhouette_stat; 

            if exist('bicrepeats','var')
                clust.bic.nrepeats = bicrepeats; 
            end
            if exist('bic','var')
                clust.bic.vals = bic; 
                clust.bic.cluster_rng_test = bic_cluster_rng; 
            end
           
 
 
        case 'Spectral'
            
            mxcl = params.numclust;
            [idx,V,D] = spectralcluster(traces1_y_norm_trunc',mxcl,'Distance','euclidean','SimilarityGraph','knn',...
                'ClusterMethod','kmedoids','LaplacianNormalization','symmetric','KNNGraphType','complete');
            

%             % Get silhouette values for cluster optimality 
%             figure; 
%             [s,h] = silhouette(traces1_y_norm_trunc', idx);
%             silhouette_stat.single_run.svals = s; 
%             silhouette_stat.single_run.clust_range = 1:params.numclust; 
%             test_clust_range = 2:15; 
%             
%             silhouette_stat.niter = 50; 
%             evals = cell(1,silhouette_stat.niter);
%             myfun = @(X,K) (spectralcluster(X,K,'Distance','euclidean','SimilarityGraph','knn',...
%                 'ClusterMethod','kmedoids','LaplacianNormalization','symmetric','KNNGraphType','complete'));
%             parfor_progress(silhouette_stat.niter); 
%             parfor niter=1:silhouette_stat.niter    
%                 evaluation = evalclusters(traces1_y_norm_trunc,myfun,'Silhouette',...
%                     'Distance','euclidean','klist',test_clust_range); 
%                 evals{niter} = evaluation; 
%                 parfor_progress; 
%             end
%             parfor_progress(0); 
%             silhouette_stat.test_clust_range = test_clust_range; 
%             silhouette_stat.evals = evals; 



            % Generate figure
            hf = figure; set(hf,'position',[8  140  1391  657]);
            for k=1:mxcl
                matrixmat = traces1_y_norm_trunc(:,find(idx==k))';
                rng_rec = [prctile(matrixmat(:),1) prctile(matrixmat(:),0.99)];  
                subplot(2,5,k);
                imagesc(matrixmat); axis square; 
                title(sprintf('Cluster %d',k));
            end
            if ~isempty(which('sgtitle')) % accomodating older Matlab versions 
                sgtitle(sprintf('Method: %s','Spectral')); 
            end
            
            % Reassign cluster id 
            reassign_clustid = NaN(size(traces1_y_norm,2),1); 
            reassign_clustid(non_outliers) = idx; 
            
            
            clust.type = params.cluster_method; 
            clust.numclust = params.numclust;
            clust.cell_ind = 1:numel(resp.cell); 
            clust.clust_id = reassign_clustid; 
            clust.non_outliers = non_outliers; 
            clust.outliers = outliers; 
            if exist('silhouette_stat','var')
                clust.silhouette_stat = silhouette_stat; 
            end


        case 'Wavelet'

            mxcl = params.numclust;
            S = mdwtcluster(traces1_y_norm_trunc','maxclust',mxcl,'pdist','Minkowski','linkage','ward');
            IdxCLU = S.IdxCLU;
            idx = IdxCLU(:,1); 
            
           
            % Get silhouette values for cluster optimality 
            figure; 
            [s,h] = silhouette(traces1_y_norm_trunc', idx);
            silhouette_stat.single_run.svals = s; 
            silhouette_stat.single_run.clust_range = 1:params.numclust; 
            test_clust_range = 2:15; 
            
            myfun = @(X,K) (spectralcluster(X,K,'Distance','correlation','SimilarityGraph','knn',...
                'ClusterMethod','kmeans','LaplacianNormalization','symmetric','KNNGraphType','complete'));

            silhouette_stat.niter = 50; 
            evals = cell(1,silhouette_stat.niter); 
%             myfun = @(X,K) (spectralcluster(X,K,'Distance','correlation','SimilarityGraph','knn',...
%                 'ClusterMethod','kmeans','LaplacianNormalization','symmetric','KNNGraphType','complete'));
            parfor_progress(silhouette_stat.niter); 
            for niter=1:silhouette_stat.niter 
                evaluation = evalclusters(traces1_y_norm_trunc',myfun_,'Silhouette',...
                    'Distance','correlation','klist',test_clust_range); 
                evals{niter} = evaluation; 
                parfor_progress; 
            end
            parfor_progress(0); 
            silhouette_stat.test_clust_range = test_clust_range; 
            silhouette_stat.evals = evals; 
        
% 
%             function [idx_] = myfun_(X,K)
%                 S_ = mdwtcluster(X,'maxclust',K,'pdist','Minkowski','linkage','ward');
%                 IdxCLU_ = S_.IdxCLU;
%                 idx_ = IdxCLU_(:,1);
%             end

    
            % Generate figure
            hf = figure; set(hf,'position',[8  140  1391  657]);
            for k=1:mxcl
                matrixmat = traces1_y_norm_trunc(:,IdxCLU(:,1)==k)'; 
                subplot(2,5,k);
                imagesc(matrixmat); axis square; 
                title(sprintf('Cluster %d',k));
            end
            if ~isempty(which('sgtitle')) % accomodating older Matlab versions 
                sgtitle(sprintf('Method: %s','Wavelet')); 
            end
            

            % Reassign cluster id 
            reassign_clustid = NaN(size(traces1_y_norm,2),1); 
            reassign_clustid(non_outliers) = idx; 
            
            
            clust.type = params.cluster_method; 
            clust.numclust = params.numclust;
            clust.cell_ind = 1:numel(resp.cell); 
            clust.clust_id = reassign_clustid; 
            clust.non_outliers = non_outliers; 
            clust.outliers = outliers; 
            clust.silhouette_stat = silhouette_stat;
 

        case 'DTW' % dynamic time warping 


            mxcl = params.numclust;
            [idx,c,sumd,d] = kmedoids(traces1_y_norm_trunc',mxcl,'Distance',@dtwf); 
            
            % Get silhouette values for cluster optimality 
            figure; 
            [s,h] = silhouette(traces1_y_norm_trunc', idx);
            silhouette_stat.single_run.svals = s; 
            silhouette_stat.single_run.clust_range = 1:params.numclust; 
            test_clust_range = 2:15; 

            myfun = @(X,K) (kmedoids(X,K,'Distance',@dtwf));

            silhouette_stat.niter = 50; 
            evals = cell(1,silhouette_stat.niter); 
            parfor_progress(silhouette_stat.niter); 
            for niter=1:silhouette_stat.niter 
                evaluation = evalclusters(traces1_y_norm_trunc',myfun,'Silhouette',...
                    'Distance','correlation','klist',test_clust_range); 
                evals{niter} = evaluation; 
                parfor_progress; 
            end
            parfor_progress(0); 
            silhouette_stat.test_clust_range = test_clust_range; 
            silhouette_stat.evals = evals; 

            
            % Generate figure
            hf = figure; set(hf,'position',[8  140  1391  657]);
            for k=1:mxcl
                matrixmat = traces1_y_norm_trunc(:,idx==k)'; 
                subplot(2,5,k);
                imagesc(matrixmat); axis square; 
                title(sprintf('Cluster %d',k));
            end
            if ~isempty(which('sgtitle')) % accomodating older Matlab versions 
                sgtitle(sprintf('Method: %s','Wavelet')); 
            end


            % Reassign cluster id 
            reassign_clustid = NaN(size(traces1_y_norm,2),1); 
            reassign_clustid(non_outliers) = idx; 
            
            clust.type = params.cluster_method; 
            clust.numclust = params.numclust;
            clust.cell_ind = 1:numel(resp.cell); 
            clust.clust_id = reassign_clustid; 
            clust.non_outliers = non_outliers; 
            clust.outliers = outliers; 
            clust.silhouette_stat = silhouette_stat;



    end


    function dist_ = dtwf(x_,y_)
        m2 = size(y_,1);
        dist_ = zeros(m2,1);
        for im=1:m2
            dist_(im) = dtw(x_,y_(im,:));
        end
    end
 
end


