function [ clustering ] = k_means_rows( Q )
    clustering_pre = kmeans(Q, 2, 'MaxIter', 1000);
    clustering = 2 * (clustering_pre - 1) - 1;
end
