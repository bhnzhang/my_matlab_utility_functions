% example code on using the video writer

clear; clc; close all;

% prepare a new AVI file
v = VideoWriter('peaks_slow.avi');
v.FrameRate = 10;
open(v);

% Generate initial data
Z = peaks;
surf(Z);
axis tight manual
set(gca, 'nextplot', 'replacechildren');    % clears figure

% Create set of frames and write each frame to the file
for k = 1:20
    
   surf( sin(2*pi*k/20)*Z, Z );
   frame = getframe( gcf );
   writeVideo(v, frame);
    
end

% close the video file
close(v);