close all
% Generate a RICE logo for the image

% Assuming a 26 x 128 grid

im = zeros(26, 128);

% R 
im(9:17, 8) = 1;
im(9, 8:20) = 1;
im(9:11, 20) = 1;
im(11, 8:20) = 1;

%%% R Diagonal
im(12, 9:10) = 1;
im(13, 11:12) = 1;
im(14, 13:14) = 1;
im(15, 15:16) = 1;
im(16, 17:18)= 1;
im(17, 19:20) = 1;
%%% 

% I
im(9:17, 40) = 1;

% C 
im(9:17, 72) = 1;
im(9, 72:84) = 1;
im(17, 72:84) = 1;


% E 
im(9:17, 108) = 1;
im(9, 108:120) = 1;
im(17, 108:120) = 1;
im(13, 108:120) = 1;



imagesc(im);