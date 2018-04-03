function [centers, membershipmat] = fcm1d(D, number_of_clusters, exponent, iteration_cap,  improvement_thresh)
%  Compute 1D fuzzy c-mean clustering centers.
%
%  Version info:
%  $Rev: 450 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2015-11-20 11:34:20 -0500 (Fri, 20 Nov 2015) $:    Date of last commit
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

% exponent can't be 1 in this implementation but papers almost always use 2 anyways
objective_function = zeros(iteration_cap, 1);	% Array for objective function
centers=zeros(number_of_clusters,1);
for i=1:number_of_clusters
    centers(i)=min(D)+(max(D)-min(D))/(number_of_clusters-1)*(i-1);
end
dist = distance(centers, D);       % compute distance from cluster centers      
membershipmat =(dist.^(-2/(exponent-1)))./(ones(number_of_clusters, 1)*sum(dist.^(-2/(exponent-1))));
membershipmat(isnan(membershipmat))=1; %to handle rare cases of infinity/infinity

% now iteratively update
for i = 1:iteration_cap,
	centers = ((membershipmat.^exponent)*D)./(ones(size(D, 2), 1)*sum((membershipmat.^exponent),2));
    dist = distance(centers, D); 
    objective_function(i) = sum(sum((dist.^2).*(membershipmat.^exponent)));
    membershipmat =(dist.^(-2/(exponent-1)))./(ones(number_of_clusters, 1)*sum(dist.^(-2/(exponent-1))));
	% check termination condition
	if i > 1,
		if abs(objective_function(i) - objective_function(i-1)) < improvement_thresh, break; end,
	end
end


function out = distance(centers, D)
out = zeros(size(centers, 1), size(D, 1));
for k = 1:size(centers, 1),
    out(k, :) = abs(centers(k)-D)';
end
