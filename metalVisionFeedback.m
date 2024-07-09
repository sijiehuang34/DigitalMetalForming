function metalVisionFeedback(bw, B, L) %#codegen

% Input bw: the binarized image, 
% Input B: the boundary coordinates for every bounded shape
% Input L: the label matrix of each shape
% No output (final image will be displayed on a fig window automatically)

% Display each shape in a different color
% 0.5 sets the background gray
% DEBUG NOTES: 2nd argument 'colormap'; must use standard syntax nx3 double
n = max(L(:)); % Number of unique labels in L (THIS ENSURES EVERY LABEL GETS A UNIQUE COLOR)
cmap = jet(n); % Generate a jet colormap with n colors
imshow(label2rgb(L,cmap,[.5 .5 .5]));

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
    % DEBUG NOTES: Must use standard syntax to ensure output is scalar
    maxArea = max([stats.Area], [], "all");
    minArea = min([stats.Area], [], "all");
    if area ~= maxArea && area ~= minArea
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
        % DEBUG NOTES: Only supports round(x) syntax
        angle_dev = round(rad2deg(atan(base / height)));
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
% DEBUG NOTES: Cannot directly access a field
% [stats.Boundary] = deal([]);

boundaries = cell(length(stats), 1);

for i = 1:length(stats)
    area = stats(i).Area;

    % BY THE SCANNING ORDER OF THE LABEL MATRIX: LEFT TO RIGHT

    % 1 = PART
    if area == min([stats.Area], [], "all")
        boundaries{i} = B{1};
    % 3 = DIE
    elseif area == max([stats.Area], [], "all")
        boundaries{i} = B{3};
    % 2 = PUNCH
    else 
        boundaries{i} = B{2};
    end
    
end

% for i = 1:length(stats)
%     area = stats(i).Area;
% 
%     % BY THE SCANNING ORDER OF THE LABEL MATRIX: LEFT TO RIGHT
% 
%     % 1 = PART
%     if area == min([stats.Area], [], "all")
%         stats(i).Boundary = B{1};
%     % 3 = DIE
%     elseif area == max([stats.Area], [], "all")
%         stats(i).Boundary = B{3};
%     % 2 = PUNCH
%     else 
%         stats(i).Boundary = B{2};
%     end
%     
% end

% Find midpoints of the part
part_midpoints = []; % initialize a matrix to store the midpoint coords

totalCol = length(bw);

for currentCol = 1:totalCol % loop through every column of the image
    % Find if the current column is in the x-coordinates of the part
    if any(B{1}(:,2) == currentCol)
        % Find all corresponding y-coords for this x-coord (column)
        index = B{1}(:,2) == currentCol;
        y_coords = B{1}(index, 1);
        % Define & store midpoints
        midY = ( max(y_coords) + min(y_coords) ) / 2;
        midPt = [currentCol, midY];
        part_midpoints = [part_midpoints; midPt];
    end
end

% Draw the midpoints
for i = 5:10:length(part_midpoints)
    plot(part_midpoints(i,1), part_midpoints(i,2), 'r.', 'MarkerSize', 10)
end

hold off;

% NOT SUPPORTED BY C GENERATION
% % Capture the figure window as an image
% frame = getframe(gcf); % gcf gets the current figure handle///fig gets the entire figure
% % Convert the frame to an image matrix
% processedImgMat = frame2im(frame);