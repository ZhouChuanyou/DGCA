% ===================================================================================
% An Effective Method for Digital Rock Reconstruction with Enhanced Pore Connectivity
% Author: CY Zhou
% Date: 2025.05.13
% Description: Processes raw digital rock data, removes small isolated pores,
%              and connects remaining isolated pore clusters to the main pore.
% ===================================================================================
clear all;clc;close all;dbstop if error;
tic;
id = fopen('DigitalRock.raw','r','b'); % Basic digital rock file
data = fread(id,'uint8');
fclose(id);
matrix = reshape(data,[400 400 400]);

% Display pores in the initial matrix (0 represents a pore)
% Connected-component analysis: identify and label all connected pore regions
CC = bwconncomp(matrix == 0, 18);
% Calculate the number of voxels (size) in each connected component
voxelSizes = cellfun(@numel, CC.PixelIdxList);
% Set a threshold to remove small isolated pores, such as pores with fewer than 50 voxels
minVoxelSize = 26;  % Adjust according to the actual situation
smallPores = find(voxelSizes < minVoxelSize);  % Find all isolated small pores

% Remove small pores
for i = 1:length(smallPores)
    matrix(CC.PixelIdxList{smallPores(i)}) = 1;  % Set these pores as matrix
end

% Perform connected-component analysis again to process the remaining pores
CC = bwconncomp(matrix == 0, 18);
Sizes = size(matrix);
NewCC = true; % CC is genuinely updated each time
% Use a while loop to process all isolated pores
while length(CC.PixelIdxList) > 1
    voxelSizes = cellfun(@numel, CC.PixelIdxList);
    [maxPoreSize, maxPoreIndex] = max(voxelSizes);
    maxPore = CC.PixelIdxList{maxPoreIndex};    

    if NewCC
        isolatedPoreIndex = find((1:length(CC.PixelIdxList)) ~= maxPoreIndex, 1);
    else
        isolatedPoreIndex = find((1:length(CC.PixelIdxList)) ~= maxPoreIndex, 2);
        isolatedPoreIndex = isolatedPoreIndex(2);
        NewCC = true;
    end
    
    if isempty(isolatedPoreIndex)
        break;
    end

    % Get the voxel coordinates of the isolated pore
    poreVoxels = CC.PixelIdxList{isolatedPoreIndex};
    [x, y, z] = ind2sub(size(matrix), poreVoxels);

    % Calculate and round the equivalent center point of the isolated pore
    centerX = ceil(mean(x));
    centerY = ceil(mean(y));
    centerZ = ceil(mean(z));
    % Compare (centerX, centerY, centerZ) with all (x, y, z) combinations
    % Determine whether the center point is actually inside the isolated pore
    isCenterInPoreVoxels = any(x == centerX & y == centerY & z == centerZ);
    if ~isCenterInPoreVoxels
        centerX = x(1); 
        centerY = y(1); 
        centerZ = z(1);
    end

    % Calculate and round the equivalent radius of the isolated pore
    numVoxels = numel(poreVoxels);  % Pore volume
    equivalentRadius = ceil((3 * numVoxels / (4 * pi))^(1/3));  % Calculate equivalent radius

    % Find the nearest voxel in the largest connected domain
    nearestPore = find_nearest_pore(Sizes, maxPore, centerX, centerY, centerZ);

    % Connect with a random curve of tortuosity 1.1; line width is half the equivalent sphere radius
    lineWidth = ceil(equivalentRadius / 2);
    if lineWidth < 8
        lineWidth = 8;
    end

    % Perform the connection
    matrix = connect_pores_with_curvature(matrix, [centerX, centerY, centerZ], nearestPore, lineWidth);

    PreCC=CC;
    % Update the largest pore system after each connection
    CC = bwconncomp(matrix == 0, 18);
    if isequal(PreCC,CC)
        NewCC = false;
    end    
end
toc;

% Save the processed three-dimensional matrix
saveRawData('DigitalRock2.raw', matrix);
save('DigitalRock2.mat', 'matrix');
% Print the porosity before and after processing
PoreVolume = sum(data == 0);  % Count the number of pores before connection
totalVolume = numel(data);
Porosity = PoreVolume / totalVolume;  % Calculate porosity
fprintf('Porosity before connection: %.4f\n', Porosity);

finalPoreVolume = sum(matrix(:) == 0);  % Count the number of pores after connection
finalPorosity = finalPoreVolume / totalVolume;
fprintf('Porosity after connection: %.4f\n', finalPorosity);

outInteger = matlab2geoeas(matrix);
% Create and open a file for writing data
fileID = fopen('DigitalRock2.out', 'w');

% Write the file header
fprintf(fileID, 'DigitalRock2\n');
fprintf(fileID, '1\n');               
fprintf(fileID, 'Realization\n'); 

% Write matrix data (one row at a time)
fprintf(fileID, '%d\n', outInteger);

% Close the file
fclose(fileID);

%% Function to find the nearest pore
function nearestPore = find_nearest_pore(matrixSize, maxPore, centerX, centerY, centerZ)
    % Get the voxel coordinates of the largest connected domain
    [conn_x, conn_y, conn_z] = ind2sub(matrixSize, maxPore);

    % Calculate and round distances from the equivalent center point to the largest pore voxels
    distances = ceil(sqrt((conn_x - centerX).^2 + (conn_y - centerY).^2 + (conn_z - centerZ).^2));

    % Find the largest-pore voxel coordinates corresponding to the minimum distance
    [~, idx] = min(distances);
    nearestPore = [conn_x(idx), conn_y(idx), conn_z(idx)];
    nearestPoreW = nearestPore;
    % Get points in the nearby region (centered at the nearest point, search within a specified range)
    searchRadius = 15;
    nearbyPores = [];

    for i = 1:length(conn_x)
        if ceil(sqrt((conn_x(i) - nearestPore(1))^2 + (conn_y(i) - nearestPore(2))^2 + (conn_z(i) - nearestPore(3))^2)) <= searchRadius
            nearbyPores = [nearbyPores; conn_x(i), conn_y(i), conn_z(i)];
        end
    end

    % Calculate and round the center point of the nearby points
    if ~isempty(nearbyPores)
        nearestPore = ceil(mean(nearbyPores, 1));
    end
    % Check whether nearestPore is in the conn_x, conn_y, conn_z lists
    inList = any(conn_x == nearestPore(1) & conn_y == nearestPore(2) & conn_z == nearestPore(3));
    if ~inList
        nearestPore = nearestPoreW;
    end
end



%% Function to connect isolated pores (connection line + line-width control)
function matrix = connect_pores_with_curvature(matrix, startPore, endPore, lineWidth)
    [x1, y1, z1] = deal(startPore(1), startPore(2), startPore(3));
    [x2, y2, z2] = deal(endPore(1), endPore(2), endPore(3));

    % Calculate and round the straight-line distance between the two points
    distance = ceil(sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2));

    % % Generate the connection line
    numPoints = 2*distance;
    ConnPoints = generate_random_curve([x1, y1, z1], [x2, y2, z2], numPoints, 1.1);

    % Generate a connection of width lineWidth along the line; use lineWidth/2 as the radius
    for i = 1:size(ConnPoints, 1)
        currentPoint = ConnPoints(i, :);
        matrix = create_sphere(matrix, currentPoint(1), currentPoint(2), currentPoint(3), ceil(lineWidth / 2));
    end
end

%% Function to generate a random curve with a target tortuosity
function GenConnPoints = generate_random_curve(startPoint, endPoint, numPoints, tortuosity)
    t = linspace(0, 1, max(2, numPoints))';
    direction = endPoint - startPoint;
    straightDistance = norm(direction);

    if straightDistance == 0
        GenConnPoints = repmat(startPoint, numel(t), 1);
        return;
    end

    direction = direction / straightDistance;

    % Construct two random unit vectors perpendicular to the endpoint direction.
    randomVector = randn(1, 3);
    randomVector = randomVector - dot(randomVector, direction) * direction;
    if norm(randomVector) < eps
        randomVector = null(direction)';
        randomVector = randomVector(1, :);
    end
    perpendicular1 = randomVector / norm(randomVector);
    perpendicular2 = cross(direction, perpendicular1);
    perpendicular2 = perpendicular2 / norm(perpendicular2);

    % Use a smooth, endpoint-preserving random offset. The amplitude is
    % calibrated below so that the polyline length is 1.1 times the
    % straight-line distance.
    randomAngle = 2 * pi * rand;
    offsetDirection = cos(randomAngle) * perpendicular1 + ...
        sin(randomAngle) * perpendicular2;

    lowerAmplitude = 0;
    upperAmplitude = straightDistance;
    for iteration = 1:40
        amplitude = (lowerAmplitude + upperAmplitude) / 2;
        candidatePoints = startPoint + t * direction * straightDistance + ...
            amplitude * sin(pi * t) .* offsetDirection;
        segmentLengths = sqrt(sum(diff(candidatePoints, 1, 1).^2, 2));
        candidateTortuosity = sum(segmentLengths) / straightDistance;

        if candidateTortuosity < tortuosity
            lowerAmplitude = amplitude;
        else
            upperAmplitude = amplitude;
        end
    end

    amplitude = (lowerAmplitude + upperAmplitude) / 2;
    GenConnPoints = startPoint + t * direction * straightDistance + ...
        amplitude * sin(pi * t) .* offsetDirection;
    GenConnPoints(1, :) = startPoint;
    GenConnPoints(end, :) = endPoint;
end

%% Function to generate spherical voxels with a specified width along the path
function matrix = create_sphere(matrix, xc, yc, zc, radius)
    [rows, cols, depths] = size(matrix);    
    for x = max(1, ceil(xc - radius)):min(rows, ceil(xc + radius))
        for y = max(1, ceil(yc - radius)):min(cols, ceil(yc + radius))
            for z = max(1, ceil(zc - radius)):min(depths, ceil(zc + radius))
                if ceil(sqrt((x - xc)^2 + (y - yc)^2 + (z - zc)^2)) <= radius
                    matrix(x, y, z) = 0;
                end
            end
        end
    end
end

%% Custom function: save data in .raw format
function saveRawData(filename, data)
    fileID = fopen(filename, 'w');
    fwrite(fileID, data, 'uint8');
    fclose(fileID);
end
