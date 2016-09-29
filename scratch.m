% scratch space

%% Example of how to add NaN separators to the triangulation points array
src = [0,1,2, ...
    4,5, ...
    7,8,9, ...
    10, 11, 12, 13];
len = [3, 2, 3, 4];
dst = nan(length(src)+length(len), 1);

i0 = 1; 
j0 = 1;
for kk = 1:length(len)
    i1 = i0+len(kk)-1;
    j1 = j0+len(kk)-1;
    
    fprintf('src(%d): %d, %d, dst(%d): %d, %d\n', ...
        length(src), i0, i1, length(dst), j0, j1);
    dst(j0:j1) = src(i0:i1);
    
    i0 = i1+1;
    j0 = j1+2;
end
    