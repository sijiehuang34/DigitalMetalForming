%% Method1_ThreeCases.m
clc; clear all; 

% Case 1: Perfect top view (Acquisition500)
% Case 2: Showing underside (Acquisition660)
% Case 3: Showing tool head (Acquisition1270)
img = imread("Acquisition500.jpg");

% Turn into grayscale for pixel histogram analysis
gray_img = im2gray(img);

% Binarize (retains top side of the punch only)
bw = imbinarize(gray_img); 

% Denoise (removes small objects)
minSize = 175; % pxl threshold through trial & error
bw = bwareaopen(bw,minSize);

% Fill holes (regionprops is used to estimate the area enclosed later)
bw = imfill(bw,"holes");

% Make known, fixed-position, irrelevant regions = 0 pxl
bw(1:60, :) = 0;
bw(1:200, 500:end) = 0;

% Now denoise again (removes any remaining small objects)
minSize = 175; % pxl threshold through trial & error
bw = bwareaopen(bw,minSize);

% Look for exterior boundaries
    % "noholes" = don't search for inner contours
    % B = a cell storing the boundary coordinates for every bounded shape
    % L = a label matrix of each shape
[B,L] = bwboundaries(bw,"noholes");

% Display each shape in a different color
    % 0.5 sets the background gray
    % @cool = a type of color map
imshow(label2rgb(L,@cool,[.5 .5 .5]))

hold on

% Make the boundaries white for visualization

for i = 1:length(B) % for i-many shapes
  boundary = B{i}; % the boundary coordinates for the i-th shape
  plot(boundary(:,2),boundary(:,1),"w",LineWidth=2) % plot coord. in white
end

% Obtain shape properties
stats = regionprops(L,"Extrema","Centroid","Area");

% Single out the punch, store & display properties
for i = 1:length(stats)

    % Area (pxl)
    area = stats(i).Area;
    if area ~= max([stats.Area]) && area ~= min([stats.Area]) 
    % should work even when punch isn't visible

        % Extremas
        extrema = stats(i).Extrema; % get the extrema for the i-th region
        plot(extrema(:,1), extrema(:,2), 'rx', 'LineWidth', 1, ...
            "MarkerSize", 8); % plot the extrema points in red

        % Outline the leftmost edge
        x_Lbottom = stats(i).Extrema(7,1); % left-bottom x
        x_Ltop = stats(i).Extrema(8,1); % left-top x
        y_Lbottom = stats(i).Extrema(7,2); % left-bottom y
        y_Ltop = stats(i).Extrema(8,2); % left-top y
        line([x_Lbottom, x_Ltop], [y_Lbottom, y_Ltop], 'Color', 'k', ...
            'LineWidth', 2);

        % Calculate the angle deviation
            % Positive angle: tilted towards the right (anticlockwise)
            % Negative angle: tilted towards the left (clockwise)
        height = y_Lbottom - y_Ltop;
        base = x_Lbottom - x_Ltop; 
        angle_dev = round(rad2deg(atan(base / height)), 2);
        angleString = "Angle deviation (Punch): " + angle_dev + "°";
        text(400, 50, angleString, 'FontSize', 10, 'Color', 'k', ...
            'FontWeight', 'bold');

        % Centroid
        centroid = stats(i).Centroid;
        plot(centroid(1,1), centroid(1,2), 'r*', 'LineWidth', 1, ...
            "MarkerSize", 8)
        centroidString = sprintf('Centroid (Punch): (%.4f, %.4f)', ...
            centroid(1,1), centroid(1,2));
        text(400, 25, centroidString, 'FontSize', 10, 'Color', 'k', ...
            'FontWeight', 'bold');

    end
end

% Match & append the boundary data to the stats structure
[stats.Boundary] = deal([]);

for i = 1:length(stats)
    area = stats(i).Area;

    % BY THE SCANNING ORDER OF THE LABEL MATRIX: LEFT TO RIGHT

    % 1 = PART
    if area == min([stats.Area])
        stats(i).Boundary = B{1};
    % 3 = DIE
    elseif area == max([stats.Area])
        stats(i).Boundary = B{3};
    % 2 = PUNCH
    else 
        stats(i).Boundary = B{2};
    end
    
end

% Find midpoints of the part
part_midpoints = []; % initialize a matrix to store the midpoint coords

totalCol = length(bw);

for currentCol = 1:totalCol % loop through every column of the image
    % Find if the current column is in the x-coordinates of the part
    if any(B{1}(:,2) == currentCol)
        % Find all corresponding y-coords for this x-coord (column)
        index = find(B{1}(:,2) == currentCol);
        y_coords = B{1}(index, 1);
        % Define & store midpoints
        midY = ( max(y_coords) + min(y_coords) ) / 2;
        midPt = [currentCol, midY];
        part_midpoints = [part_midpoints; midPt];
    end
end

% Draw the midpoints
for i = 1:10:length(part_midpoints)
    plot(part_midpoints(i,1), part_midpoints(i,2), 'r.', 'MarkerSize', 10)
end

hold off;


