F = plottingFromSavedFiles;

v = VideoWriter('./movies/doubleStraka12p5.avi');
v.FrameRate = 10;
v.Quality = 100;

open(v);
writeVideo(v,F);
close(v);
