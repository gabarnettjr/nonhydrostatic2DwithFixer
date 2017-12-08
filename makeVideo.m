F = plottingFromSavedFiles;

v = VideoWriter('./movies/bubbleTopoSmooth25.avi');
v.FrameRate = 10;
v.Quality = 100;

open(v);
writeVideo(v,F);
close(v);
