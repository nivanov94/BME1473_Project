%
% new_data = laplacian_filter(data, locations, nbrs)
%
% input:
%   data      - (CxD) matrix of eeg data to downsample where
%                C is the number of channels
%                D is the number of samples per epoch
%   locations - EEG channel locations
%   nbrs      - number of nearest neighbors to consider
% output:
%   new_data - (CxD) matrix of filtered input data
%

function [new_data] = laplacian_filter(data, locations, nbrs) 

    C = size(data, 1); % Number of channels
    D = size(data, 2); % Number of samples per epoch

    new_data = data;

    % Get channel locations in 3D space
    x = [ locations.X ];
    y = [ locations.Y ];
    z = [ locations.Z ];

    % Create the distance matrix for each channel
    dist_mat = zeros(C);
    for i=1:C
        for j=1:C
            dist_mat(i,j) = sqrt( (x(i)-x(j))^2 + (y(i)-y(j))^2 + (z(i)-z(j))^2 );
        end
    end

    % Get the n closest channels
    [~, sorted_index] = sort(dist_mat, 2);
    nearest_nbr = num2cell(sorted_index(:,2:(nbrs+1)), 2);  

    % Get the normalized nearest neighbor inverse distance matrix
    g = zeros(C, nbrs);
    for i=1:C
        g(i,:) = 1 ./ dist_mat(i,nearest_nbr{i});
        g(i,:) = g(i,:) ./ sum(g(i,:));
    end
    
    for d=1:D
        % Take the sum of n-closest channel values weighted by normalized inverse distance
        sum_volt = zeros(C, 1);
        for c=1:C
            sum_volt(c) = g(c, :) * new_data(nearest_nbr{c}, d);
        end
    end
    
    % Subtract this from the sample
    new_data(:, d) = new_data(:, d) - sum_volt;
end
