clear; clc;

Problem_1();

function Problem_1()
% Template for EE535 Digial Image Processing
% Insert the code in the designated area below
%% Loading directory for image files
imgdir = uigetdir('D:/KAIST/Courses/dip/hw/hw6/Test_images');
file = fopen(fullfile(imgdir,'\shepp_proj_180x256.raw'),'rb');
gray_image = fread(file,fliplr([180,256]),'*float')';
fclose(file);
%%
%%---------------------Insert code below ----------------------%%
% 1.1. Perform the back projection
bp_img = BackProjection(gray_image,256,180);

% 1.2. Perform the filtered back projection
ram_lak_sinogram_180 = RamLakFilter(gray_image,180);
ram_lak_bp_180 = BackProjection(ram_lak_sinogram_180,256,180);

shepp_logan_sinogram_180 = SheppLoganFilter(gray_image,180);
shepp_logan_bp_180 = BackProjection(shepp_logan_sinogram_180,256,180);

% 1.3. Compare results
bp_cv = CutView(bp_img);
ram_lak_cv = CutView(ram_lak_bp_180);
shepp_logan_cv = CutView(shepp_logan_bp_180);

% 1.4. Repeat 1.2 with number of views equal to 60
ram_lak_sinogram_60 = RamLakFilter(gray_image,60);
ram_lak_bp_60 = BackProjection(ram_lak_sinogram_60,256,60);

shepp_logan_sinogram_60 = SheppLoganFilter(gray_image,60);
shepp_logan_bp_60 = BackProjection(shepp_logan_sinogram_60,256,60);

%% Displaying figures
% Problem 1.1
figure('Name','Problem 1.1 - Perform the back projection');
subplot(1,2,1); imshow(gray_image,[]); title('The original image');
subplot(1,2,2); imshow(abs(bp_img)/255,[]); title('Back projection image');

% Problem 1.2
figure('Name','Problem 1.2 - Perform the filtered back projection');
subplot(1,2,1); imshow(ram_lak_bp_180,[0.44 0.48]); title('Ram-Lak back projection');
subplot(1,2,2); imshow(shepp_logan_bp_180,[0.44 0.48]); title('Shepp-Logan back projection');

% Problem 1.3
figure('Name','Problem 1.3 - Cut View');
subplot(1,3,1); imshow(bp_cv,[]); title('Cutview of back projection image');
subplot(1,3,2); imshow(ram_lak_cv,[]); title('Cutview of Ram-Lak back projection image');
subplot(1,3,3); imshow(shepp_logan_cv,[]); title('Cutview of Shepp-Logan back projection image');

% Problem 1.4
figure('Name','Problem 1.4 - Perform the filtered back projection ');
subplot(2,2,1); imshow(ram_lak_bp_180,[0.44 0.48]); title('Ram-Lak with 180 views');
subplot(2,2,2); imshow(shepp_logan_bp_180,[0.44 0.48]); title('Shepp-Logan with 180 views');
subplot(2,2,3); imshow(ram_lak_bp_60+0.855,[]); title('Ram-Lak with 60 views');
subplot(2,2,4); imshow(shepp_logan_bp_60+0.855,[]); title('Shepp-Logan with 60 views');

%%---------------------------------------------------------------%%
end

%% Inner Function
% Back Projection
function A = BackProjection(X,Ns,Nv)
	[~,n] = size(X);
	A = zeros(Ns,Ns);
  t = round(180/Nv);
  theta_d = (1:t:180)';
  len_theta = length(theta_d);
  theta_r = theta_d*pi/180;
	for i = -n/2:n/2-1
    for j = -n/2:n/2-1
      s = j*cos(theta_r) + i*sin(theta_r);
      s = round(s+n/2+1);
      for k = 1:len_theta
        if s(k)>0 && s(k)<=n
          A(n/2-i,j+n/2+1) = A(n/2-i,j+n/2+1) + X(k,s(k));
        end
      end
    end
  end
end

% Ram-Lak Filter
function A = RamLakFilter(X,Nv)
  [~,n] = size(X);
  len = 2 * n + 1;
  H = zeros(1,len);
  H(n+2:2:len) = -1./((1:2:n).^2*pi^2);
  H(n+1) = 1/4;
  H(n:-2:1) = H(n+2:2:len);
  
  A = zeros(Nv,n);
  for i = 1:Nv
    for j = 1:n
      K = X(i,:).*H(n+2-j:2*n+1-j);
      A(i,j) = sum(K);
    end
  end
end

% Shepp-Logan Filter
function A = SheppLoganFilter(X,Nv)
  [~,n] = size(X);
	H = -2./(pi^2*(4*(-n:n).^2-1));

	A = zeros(Nv,n);
	for i = 1:Nv
    for j = 1:n
      K = X(i,:).*H(n+2-j:2*n+1-j);
      A(i,j) = sum(K);
    end
	end
end

% Cut view
function A = CutView(X)
  [~,n] = size(X);
	A = zeros(n,n);
	for i = 1:n
    max_val = max(max(X));
    j = n-round(X(n/2,i)/max_val*n/2);
		A(j,i) = 1;
	end
end
