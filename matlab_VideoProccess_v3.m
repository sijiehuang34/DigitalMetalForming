%% Video processing

clc; clear; close all;

test_name = getTestName();

%path = '/home/dmf/Documents/GitHub/rt-dmf/SourceCode/DatasetExp/Exp3_ILCTest/Test1';
path = test_name;
speedVideo = 1;
%test_name = path;
imageDir = [test_name  '/images'];
%imageDir = 'images';
imageFiles = dir(fullfile(imageDir, 'Acquisition*.jpeg'));
[~, reindex] = sort(str2double(regexp({imageFiles.name}, '\d+', 'match', 'once')));
imageFiles = imageFiles(reindex);
numFiles = length(imageFiles);
processedImages = cell(1, numFiles);
angleVec = zeros(1, length(imageFiles));


cpp_file = [test_name  '/CPPlog.csv'];
data = readtable(cpp_file);
rows_to_replace = data.EndTimeProcess_sec_ == 0;
data(rows_to_replace, :) = array2table(NaN(sum(rows_to_replace), width(data))); % Replace entire rows with NaNs
timeVector = mergeTimeData(data.StartTime_sec_,data.StartTime_nanosec_);
timeVector = timeVector - timeVector(1);

data.BendAngle_deg_(data.BendAngle_deg_ == 0) = NaN;
validIndices = ~isnan(data.BendAngle_deg_);  % Indices of non-NaN values
time_stamp = timeVector(validIndices); % Corresponding time values
angleVec = data.BendAngle_deg_(validIndices); % Corresponding data values


%parfor i = 1:1:floor(numFiles/10)
for i = 1:1:floor(numFiles/speedVideo)
%for i =200:1:201  
    new_idx = i*speedVideo;
    filePath = fullfile(imageDir, imageFiles(new_idx).name);
    img = imread(filePath);
    local_dataset = data(new_idx,:);
   
   [processedImg,angle] = applyOpenCVOperations(img,angleVec(new_idx),time_stamp(new_idx),speedVideo); % You need to define this function
   % processedImages{i} = processedImg;
    processedImages{i} = processedImg;
   % angleVec(i) = angle;
   if mod(i, 10) == 0
            disp(['Image Convertion Progress is ', num2str((i*speedVideo/numFiles*100)*0.5), '%']);
   end
    
end

cpp_file = [test_name  '/CPPlog.csv'];
data = readtable(cpp_file);
videoSaveDir = [test_name  '/video.avi'];
v = VideoWriter(videoSaveDir);
v.FrameRate = 1/data.SimulatedDt(1);
open(v);

for i = 1:1:floor(numFiles/speedVideo)
    %% cindyfx: writeVideo does not take 'logical'
    grayscaleImage = uint8(processedImages{i}) * 255;
    writeVideo(v, grayscaleImage);
   if mod(i, 10) == 0
            disp(['Save Progress is ', num2str(50+(i*speedVideo/numFiles*100)*0.5), '%']);
   end
end

close(v);
disp(['Conversion complete. File saved at:',videoSaveDir]);

function [image,angle] = applyOpenCVOperations(image,angle,timestamp,speedVideo)
    %% cindyfx () starts;

    %% Pre-processing
    % Turn into grayscale for pixel histogram analysis
    image = im2gray(image);

    % Binarize (retains top side of the punch only)
    image = imbinarize(image); 

    % Denoise (removes small objects)
    minSize = 175; % pxl threshold through trial & error
    image = bwareaopen(image,minSize);

    % Fill holes (regionprops is used to estimate the area enclosed later)
    image = imfill(image,"holes");

    % Make known, fixed-position, irrelevant regions = 0 pxl
    image(1:60, :) = 0;
    image(1:200, 500:end) = 0;

    % Now denoise again (removes any remaining small objects)
    minSize = 175; % pxl threshold through trial & error
    image = bwareaopen(image,minSize);

    %% Look for exterior boundaries
    % "noholes" = don't search for inner contours
    % B = a cell storing the boundary coordinates for every bounded shape
    % L = a label matrix of each shape
    [B,L] = bwboundaries(image,"noholes");

    % Display each shape in a different color
    % 0.5 sets the background gray
    % @jet = a type of color map
    imshow(label2rgb(L,@jet,[.5 .5 .5]))

    hold on

    % Make the boundaries black for visualization

    for i = 1:length(B) % for i-many shapes
      boundary = B{i}; % the boundary coordinates for the i-th shape
      plot(boundary(:,2), boundary(:,1), "k", LineWidth = 3) % plot coord. in white
    end

    %% Obtain shape properties
    stats = regionprops(L,"Extrema","Centroid","Area");

    % Single out the punch, store & display properties
    for i = 1:length(stats)
    
        % Area (pxl)
        area = stats(i).Area;
        if area ~= max([stats.Area]) && area ~= min([stats.Area]) 
        % should work even when punch isn't visible
    
            % Extremas
            extrema = stats(i).Extrema; % get the extrema for the i-th region
            plot(extrema(:,1), extrema(:,2), 'rx', 'LineWidth', 2, ...
                "MarkerSize", 10); % plot the extrema points in red
    
            % Outline the leftmost edge
            x_Lbottom = stats(i).Extrema(7,1); % left-bottom x
            x_Ltop = stats(i).Extrema(8,1); % left-top x
            y_Lbottom = stats(i).Extrema(7,2); % left-bottom y
            y_Ltop = stats(i).Extrema(8,2); % left-top y
            line([x_Lbottom, x_Ltop], [y_Lbottom, y_Ltop], 'Color', 'r', ...
                'LineWidth', 2);
    
            % Calculate the angle deviation
                % Positive angle: tilted towards the right (anticlockwise)
                % Negative angle: tilted towards the left (clockwise)
            height = y_Lbottom - y_Ltop;
            base = x_Lbottom - x_Ltop; 
            angle_dev = round(rad2deg(atan(base / height)), 2);
            angleString = "Angle (Punch): " + angle_dev + "°";
            text(400, 50, angleString, 'FontSize', 14, 'Color', 'k', ...
                'FontWeight', 'bold');
    
            % Centroid
            centroid = stats(i).Centroid;
            plot(centroid(1,1), centroid(1,2), 'r*', 'LineWidth', 2, ...
                "MarkerSize", 10)
            centroidString = sprintf('Centroid (Punch): (%.2f, %.2f)', ...
                centroid(1,1), centroid(1,2));
            text(400, 25, centroidString, 'FontSize', 14, 'Color', 'k', ...
                'FontWeight', 'bold');
    
        end
    end

    %% Match & append the boundary data to the stats structure
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

    %% Find midpoints of the part
    part_midpoints = []; % initialize a matrix to store the midpoint coords
    
    totalCol = 720;
    
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

    %% Draw the midpoints
    for i = 1:10:length(part_midpoints)
        plot(part_midpoints(i,1), part_midpoints(i,2), 'r.', 'MarkerSize', 10)
    end
    
    hold off;

    %% Display encoder data for the part
    partString = "Angle (Part): " + angle + "°";
    text(20, 400, partString, 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
    
    timeString = "Time: " + timestamp + " sec";
    text(20, 500, timeString, 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');

    speedString = "Video Speed Up: " + speedVideo + "X";
    text(20, 525, speedString, 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');

    %% Let whatever is shown in the figure window become the new image
    % Capture the figure as an image
    frame = getframe(gcf); % gcf gets the current figure handle///fig gets the entire figure
    % Convert the frame to an image matrix
    imageMatrix = frame2im(frame);
    % Store the image matrix in a variable
    image = imageMatrix;

    %% cindyfx () ends;

    % image = insertText(image, [20, 60], ['Angle (Part): ', num2str(angle), '°'], 'FontSize', 18, 'TextColor', 'white', 'BoxOpacity', 0);
    % image = insertText(image, [20, 90], ['Time: ', num2str(timestamp), ' sec'], 'FontSize', 18, 'TextColor', 'white', 'BoxOpacity', 0);
    % image = insertText(image, [20, 120], ['Video Speed Up ', num2str(speedVideo), 'X'], 'FontSize', 18, 'TextColor', 'white', 'BoxOpacity', 0);

end

function showImage(image)
    figure, imshow(image); title('Canny Edge Detection on Masked Region');
    ax = gca;
    ax.Visible = 'on';
    xlabel('X-axis (px)');
    ylabel('Y-axis (px)');
end

function [coefficients, angle, isVertical] = fitLineNormalDistance(x, y)
    isVertical = var(x) < var(y); % Check if the line is nearly vertical
    if isVertical
        % The line is closer to vertical. Switch x and y for fitting.
        A = [y ones(length(y), 1)];
        b = x;
        coefficients = (A' * A) \ A' * b;
        angle =- 90-atan2(-coefficients(1), 1) * (180 / pi);
   else
        % The line is closer to horizontal. Proceed as normal.
        A = [x ones(length(x), 1)];
        b = y;
        coefficients = (A' * A) \ A' * b;
        angle = atan2(-coefficients(1), 1) * (180 / pi);
    end
end


function totalSeconds = mergeTimeData(secondsVec, nanosecondsVec)
    secondsVec(end) = [];
    nanosecondsVec(end) = [];

    if any(secondsVec < 0) || any(nanosecondsVec < 0) || any(nanosecondsVec > 999999999)
        error('Input data contains values outside the valid range.');
    end

    nanosecondsConverted = double(nanosecondsVec) * 1e-9;
    totalSeconds = double(secondsVec) + nanosecondsConverted;
end


function test_name = getTestName()
    % getTestName Reads the test name from a specified configuration file.
    %
    % INPUT:
    % config_file_path - Path to the configuration file.
    %
    % OUTPUT:
    % test_name - The test name read from the configuration file.
    config_file_path = 'datasetloc.txt';
    % Open the file
    fileID = fopen(config_file_path, 'r');

    % Check if the file was opened successfully
    if fileID == -1
        error('Cannot open config file: %s', config_file_path);
    end

    % Read the test name from the file
    test_name = fgetl(fileID);

    % Close the file
    fclose(fileID);

    % Display the test name (optional)
    disp(['Test name: ', test_name]);
end

