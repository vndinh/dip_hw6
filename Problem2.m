clear; clc;

Problem_2();

function Problem_2()
% Template for EE535 Digial Image Processing
% Insert the code in the designated area below
%% Loading directory for image files
imgdir = uigetdir('D:/KAIST/Courses/dip/hw/hw6/Test_images');
file = fopen(fullfile(imgdir,'\bridge_gray_256x256.raw'),'rb');
gray_bridge = fread(file,fliplr([256,256]),'*uint8')';
fclose(file);
%%
%%---------------------Insert code below ----------------------%%
hist_bridge = histPDF(gray_bridge);
m = length(hist_bridge);
C = Huffman(hist_bridge);
code = char(GetCode(C,m));
entropy_huff = sum(code' ~= ' ')*hist_bridge; % Huffman coding entropy
H = -(log(hist_bridge)/log(2))'*hist_bridge;  % Source entropy
compression_efficiency = 100*H/entropy_huff;  % Huffman coding compression efficiency

%% Displaying figures
figure('Name','Problem 2.1 - Huffman Coding');
subplot(1,2,1); imshow(uint8(gray_bridge)); title('The original image');
subplot(1,2,2); stem((0:1:255),hist_bridge); title('Histogram');
xlim([-20 270]); xlabel('Symbol'); ylabel('Probability');

fprintf('Original entropy = %f\n', H);
fprintf('Huffman entropy = %f\n', entropy_huff);
fprintf('Compression efficiency = %f %%\n', compression_efficiency);

fileID = fopen('huffman_code.txt','w');
fprintf(fileID,' %12s %12s \r \n','Grayscale','Codeword');
for i = 1:m
  fprintf(fileID,'\t %3d \t\t\t\t %15s \n',i-1,code(i,:));
end
%%---------------------------------------------------------------%%
end

%% Inner Function
% Histogram of image
function A = histPDF(X)
  [m,n] = size(X);
  A = zeros(256,1);
  for i = 1:256
    A(i) = sum(sum(X==i-1))/(m*n);
  end
end

% Huffman table (tree)
function A = Huffman(X)
  len = length(X);
  
  % Create 1-256 symbols corresponding with 0-255 grayscales
  A = cell(len,1);
  for i = 1:len
    A{i} = i;
  end
  
  while (length(A)>2)     % Repeat until only two branches
    [X,Idx] = sort(X);    % Sort to ascending probability
    A = A(Idx);           % Reorder Huffman table
    A{2} = {A{1}, A{2}};  % Merge 2 smallest branches
    A(1) = [];            % Remove the first branch
    X(2) = X(1) + X(2);   % Sum two smallest probabilities
    X(1) = [];            
  end
end

% Get Huffman code
function y = GetCode(a,n)
  global y
  y = cell(n,1);
  trace(a,[])
end

% Assign 0 and 1 to each pair of branches
function trace(a,codeword)
  global y
  if isa(a,'cell')
    trace(a{1},[codeword 0]);
    trace(a{2},[codeword 1]);
  else   
    y{a} = char(48 + codeword);   
  end
end
