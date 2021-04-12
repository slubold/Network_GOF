function [error] = clustering_error(Z, label)
% Computing the clustering error. Modified based on the 
% Matlab code from Universal Latent Space Model Fitting 
% for Large Networks with Edge Covariates by Zhuang Ma
% , Zongming Ma, and Hongsong Yuan.
%
% INPUT
%   Z               latent vector matrix, shape = n * p
%	label    		true label, length = n

% OUTPUT
%   error           average clustering error

n = length(label);
labelset = unique(label);
k = length(labelset);
errortemp = 0;
for index = 1:100
    [id, C] = kmeans(Z, k, 'Replicates', 200);
    P = perms(labelset);
    N = size(P, 1);
    confusion = zeros(1, N);
    for i = 1 : N
        for node = 1 : n
            if (label(node) ~= P(i, id(node)))
                   confusion(i) = confusion(i) + 1;
            end
        end
    end
    errortemp = errortemp + min(confusion) / n;
end

error = errortemp / 100;



