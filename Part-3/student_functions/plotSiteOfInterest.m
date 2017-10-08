%% function plot Site of Interest
%
% Author Kuan Chun Hwang

function [map2, siteOfinterestPos] = plotSiteOfInterest(sitesofInterest,map2,cellDim,mapDim)

% Initialize site numbers
site_n = 1;
for sites = 1:length(sitesofInterest(:,1))
    
    % Convert the positions of site of interests into cell positions
    sites_cell = pos2cell(sitesofInterest(sites,:),cellDim,mapDim);
    sites_cellY =round(sites_cell(1));
    sites_cellX =round(sites_cell(2));
    
    % The site of Interest cannot be travelled to if it is an obstacle
    if map2(sites_cellX,sites_cellY) ~= 4
        
        % Show sites of interest on map and store the site cell positions
        map2(sites_cellX,sites_cellY) = 3;
        siteOfinterestPos(site_n,:) = sitesofInterest(sites,:);
        site_n = site_n +1;
    end
    
end

end