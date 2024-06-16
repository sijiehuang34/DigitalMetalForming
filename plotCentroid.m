%% plotCentroid.m
clc; clear all; close all;

load('centroidData.mat');
load('timestampData.mat');

centroid = cell2mat(centroidData);
time = cell2mat(timestampData);

plot(time, centroid(:, 1), 'b-', 'LineWidth', 2, 'DisplayName', 'Centroid X-Coordinate');
hold on
plot(time, centroid(:, 2), 'g-', 'LineWidth', 2, 'DisplayName', 'Centroid Y-Coordinate');

xlabel('Time (second)', 'FontSize', 14);
ylabel('Centroid Coordinate (pixel)', 'FontSize', 14);
title('Centroid Coordinate Over Time', 'FontSize', 14);

set(gca, 'FontSize', 14); % Set the font size for the axis numbers

legend('show', 'FontSize', 14, 'Location', 'best'); % This will display the legend with the specified labels
hold off

