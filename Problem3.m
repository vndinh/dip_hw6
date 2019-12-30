clear; clc;

Problem_3();

function Problem_3()
% Template for EE535 Digial Image Processing
% Insert the code in the designated area below
%% Loading directory for image files
imgdir = uigetdir('D:/KAIST/Courses/dip/hw/hw6/Test_images');
file = fopen(fullfile(imgdir,'\cameraman_gray_256x256.raw'),'rb');
gray_cam = fread(file,fliplr([256,256]),'*uint8')';
fclose(file);
%%
%%---------------------Insert code below ----------------------%%
gray_cam = double(gray_cam);
% Problem 3.1. Perform 8x8-block DCT
dct8_cam = DCT2D(gray_cam,8);
dct8_mag = 20*log(abs(dct8_cam));

% Problem 3.2. Perform quantization
Q = cell(31,1);
for i = 1:31
  Q{i,1} = Quantization(dct8_cam,i);
end

Q1 = 20*log(abs(Q{1,1}));
Q16 = 20*log(abs(Q{16,1}));
Q31 = 20*log(abs(Q{31,1}));

% Problem 3.1. Perform dequantization and 8x8-block inverse DCT
deQ = cell(31,1);
invDCT8 = cell(31,1);
psnr_val = zeros(31,1);
for i = 1:31
  deQ{i,1} = Dequantization(Q{i,1},i);
  invDCT8{i,1} = IDCT2D(deQ{i,1},8);
  psnr_val(i) = PSNR(gray_cam,invDCT8{i,1});
end

deQ1 = 20*log(abs(deQ{1,1}));
deQ16 = 20*log(abs(deQ{16,1}));
deQ31 = 20*log(abs(deQ{31,1}));

%% Displaying figures
figure('Name','Problem 3.1 - Perform the 8x8 block DCT');
subplot(1,2,1); imshow(uint8(gray_cam)); title('Gray Cameraman');
subplot(1,2,2); imshow(uint8(dct8_mag)); title('8x8-block DCT');

figure('Name','Problem 3.2 - Perform quantization');
subplot(1,3,1); imshow(uint8(Q1)); title('Quantization with parameter = 1');
subplot(1,3,2); imshow(uint8(Q16)); title('Quantization with parameter = 16');
subplot(1,3,3); imshow(uint8(Q31)); title('Quantization with parameter = 31');

figure('Name','Problem 3.3 - Dequantization and inverse 8x8-block DCT');
subplot(2,3,1); imshow(uint8(deQ1)); title('Dequantization with parameter = 1');
subplot(2,3,2); imshow(uint8(deQ16)); title('Dequantization with parameter = 15');
subplot(2,3,3); imshow(uint8(deQ31)); title('Dequantization with parameter = 31');
subplot(2,3,4); imshow(uint8(invDCT8{1,1})); title('8x8-IDCT with quantization parameter = 1');
subplot(2,3,5); imshow(uint8(invDCT8{16,1})); title('8x8-IDCT with quantization parameter = 16');
subplot(2,3,6); imshow(uint8(invDCT8{31,1})); title('8x8-IDCT with quantization parameter = 31');

figure('Name','Problem 3.3 - PSNR Graph'); plot(1:31,psnr_val); grid on;
xlabel('Quantization parameters');
ylabel('PSNR [dB]');
title('PSNR graph');

%%---------------------------------------------------------------%%
end

%% Inner Function
% 2-Dimension DCT Function
function A = DCT2D(X,blkSize) % 2D DCT transform
  [m,n] = size(X);
  DCT = cos((pi/(2*blkSize))*(0:blkSize-1)'*(1:2:2*blkSize-1));
  A = zeros(m);
  for i = 1:blkSize:m-blkSize+1   
    for j = 1:blkSize:n-blkSize+1  
      A(i:i+blkSize-1,j:j+blkSize-1) = DCT*X(i:i+blkSize-1,j:j+blkSize-1)*DCT';
      A(i:i+blkSize-1,j:j+blkSize-1) = (2/blkSize)*A(i:i+blkSize-1,j:j+blkSize-1);
      A(i,j:j+blkSize-1) = A(i,j:j+blkSize-1)/sqrt(2);
      A(i:i+blkSize-1,j) = A(i:i+blkSize-1,j)/sqrt(2);       
    end
  end
end

function A = IDCT2D(X,blkSize)
  [m,n] = size(X);
  DCT = cos((pi/(2*blkSize))*(0:blkSize-1)'*(1:2:2*blkSize-1));
  IDCT = DCT';
  A = zeros(256);
  for i = 1:blkSize:m-blkSize+1   
    for j = 1:blkSize:n-blkSize+1  
      B = (2/blkSize)*X(i:i+blkSize-1,j:j+blkSize-1);
      B(1,:) = B(1,:)/sqrt(2);
      B(:,1) = B(:,1)/sqrt(2);
      A(i:i+blkSize-1,j:j+blkSize-1) = IDCT*B*IDCT';
    end
  end
end

function A = Quantization(X,p)
  A = round(X/(2*p));
end

function A = Dequantization(X,p)
  A = X * 2 * p;
end

% PSNR function
function psnr = PSNR(Orig,Dist)
	[m, n, p] = size(Orig);
	Orig = double(Orig);
	Dist = double(Dist);
	error = Orig - Dist;
	MSE = sum(sum(sum(error.^2)))/(m*n*p);
	if MSE > 0
    psnr = 20*log10(max(max(max(Orig)))) - 10*log10(MSE);
	else
    psnr = 99;
	end
end
