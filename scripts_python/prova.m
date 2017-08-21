
pkg load image

I = imread('../Ground_truth/occlusions_sintel/market_2/frame_0032.png');
I;
%imwrite(double(cat(3, I, I, I)), 'imatge_prova.png')
imwrite(uint8(I), 'imatge_prova.png')
I2 = imread('imatge_prova.png'); 