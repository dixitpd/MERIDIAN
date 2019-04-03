function [thx] = sample_with_phi(phi,mu,sg,nS)

%
% This script samples parameters kdeg and kdegs given the weights phi, the 
% locations of the Gaussian basis function mu and their widths sg
%

thx = [];
for i=1:nS
id = discretesample(phi,1);
    fx  = mvnrnd(mu(id,:),diag(sg(id,:).*sg(id,:)),5);
    thx = [thx;fx];
end

end


