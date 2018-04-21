%------------------
% Images to stitch.
%------------------
path1 = 'images/005/005.JPG';
path2 = 'images/005/006.JPG';

%-------------
% Read images.
%-------------
fprintf('Read images and SIFT matching\n');tic;
fprintf('> Reading images...');tic;
img1 = imresize(imread(sprintf('%s',path1)),1);
img2 = imresize(imread(sprintf('%s',path2)),1);
fprintf('done (%fs)\n',toc);

[ kp1,ds1 ] = sift(single(rgb2gray(img1)));
[ kp2,ds2 ] = sift(single(rgb2gray(img2)));
showkeys(img1, kp1);
matches_idx=match(kp1,ds1,kp2,ds2);
