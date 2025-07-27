clc;
clear all;
close all;

%% Parameters
M = 512;
N = 512;

x1 = linspace(-32,32,M);
y1 = linspace(-32,32,N);
[x,y] = meshgrid(x1,y1);  
%% Cell Simulation
N_cells = 5;
pixel_to_um = 3; %3;%1.5;
sigma_nx = zeros(1,N_cells);
sigma_ny = zeros(1,N_cells);
sigma_cx = zeros(1,N_cells);
sigma_cy = zeros(1,N_cells);
sigma_tx = zeros(1,N_cells);
sigma_ty = zeros(1,N_cells);
theta = zeros(1,N_cells);
for i = 1:N_cells
sigma_nx(i) = ((2.5 + (2.9-2.5) .* abs(rand(1)))/(2*3))*pixel_to_um;
sigma_ny(i) = ((1.7 + (1.9-1.7) .* abs(rand(1)))/(2*3))*pixel_to_um;
sigma_cx(i) = ((2.7 + (4.5-2.7) .* abs(rand(1)))/(2*3))*pixel_to_um;
sigma_cy(i) = ((1.9 + (2.1-1.9) .* abs(rand(1)))/(2*3))*pixel_to_um;
sigma_tx(i) = ((0.9 + (1.5-0.9) .* abs(rand(1)))/(2*3))*pixel_to_um;
sigma_ty(i) = ((0.95+ (1.05-0.95) .* abs(rand(1)))/(2*3))*pixel_to_um;
theta(i) = 360*rand(1);
end
xyRange = [-30 30] ;
minimumDistance = 60 ;
attempt_counter = 1 ; 
Distances = 0 ;
while any(Distances < minimumDistance) && attempt_counter < 100000
    attempt_counter = attempt_counter + 1 ;
    Pxy = randi(xyRange, N_cells, 2) ; 
    Distances = pdist(Pxy) ;
end
if attempt_counter < 100000
    plot(Pxy(:,1), Pxy(:,2),'bo') ;
else
    disp('No positions found.') ;
end



% a = ((cosd(Theta)^2) / (2*sigma1^2)) + ((sind(Theta)^2) / (2*sigma2^2));
% b = -((sind(2*Theta)) / (4*sigma1^2)) + ((sind(2*Theta)) / (4*sigma2^2));
% c = ((sind(Theta)^2) / (2*sigma1^2)) + ((cosd(Theta)^2) / (2*sigma2^2));

Ph = zeros(M,N,N_cells);
for i = 1: size(Pxy,1)
Ph(:,:,i) =   2*exp((-1/(sigma_nx(i)).^2)*(((x-Pxy(i))/sigma_nx(i)).^2+((y-Pxy(i,2))/sigma_ny(i)).^2)) ...
              + 2*exp((-1/(sigma_cx(i)).^2)*(((x-Pxy(i))/sigma_cx(i)).^2+((y-Pxy(i,2))/sigma_cy(i)).^2)) ...
              + 1*exp((-1/(sigma_cx(i)).^2)*(((x-Pxy(i)-(2*sigma_cx(i)))/sigma_tx(i)).^2+((y-Pxy(i,2))/sigma_ty(i)).^2)) ...
              +  1*exp((-1/(sigma_cx(i)).^2)*(((x-Pxy(i)+(2*sigma_cx(i)))/sigma_tx(i)).^2+((y-Pxy(i,2))/sigma_ty(i)).^2));
P_rotate = imrotate(Ph(:,:,i), theta(i));
P_rotate = imresize(P_rotate, [M,N]);
Ph(:,:,i) = P_rotate;
% Ph(:,:,i) =   6*exp((-1/(sigma_nx(i)).^2)*(((x-Pxy(i))/sigma_nx(i)).^2+((y-Pxy(i,2))/sigma_ny(i)).^2)) ...
%               + 4*exp((-1/(sigma_cx(i)).^2)*(((x-Pxy(i))/sigma_cx(i)).^2+((y-Pxy(i,2))/sigma_cy(i)).^2)) ...
%               + 3*exp((-1/(sigma_cx(i)).^2)*(((x-Pxy(i)-(2.5*sigma_cx(i)))/sigma_tx(i)).^2+((y-Pxy(i,2))/sigma_ty(i)).^2)) ...
%               +  3*exp((-1/(sigma_cx(i)).^2)*(((x-Pxy(i)+(2.5*sigma_cx(i)))/sigma_tx(i)).^2+((y-Pxy(i,2))/sigma_ty(i)).^2));
end
phase = sum(Ph,3);
figure, imagesc(phase);
colormap jet;axis tight; axis on;
phase_max = max(max(phase));
% figure, imagesc(wrap(phase));
% figure; plot(phase(230,:));
% line([30 94], [164 158]);
% improfile
%%
% x = linspace(-10,10,N);
% y = linspace(-10,10,M);
% [X,Y] = meshgrid(x,y);
% a2 = 6; b2 = 6; a3 = 2.4; b3 = 2.4; c1 = 8;
% aberr = exp(1i*((a2*X^2+b2*Y^2) + (a3*X^3 + b3*Y^3) + c1*X*Y));
% figure, imagesc(angle(aberr)); colormap jet;


%%
x = linspace(-0.8,0.8,1024);
[X,Y] = meshgrid(x,x);
[M,N] = size(X);
[theta,r] = cart2pol(X,Y);
idx = r<=1;
pOrder = 35;
p =0:pOrder-1;
z = zeros(size(X));
y = zernfun2(p,r(idx),theta(idx));
zMat = zeros(M,N,pOrder);
for k = 1:length(p)
  z(idx) = y(:,k);  
  zMat(:,:,k) =  z;
end
Is = zeros(M,N,3);
% phase_aberr =  10*(zMat(:,:,8)+zMat(:,:,5)+zMat(:,:,6));
% phase_aberr =  10*zMat(:,:,8)+ 5*zMat(:,:,1)+10*zMat(:,:,1) + 10*zMat(:,:,9);
% phase_aberr = sqrt(zMat(:,:,4).^2 + zMat(:,:,6).^2);
z = 25 * (zMat(:,:,13)); %10=low, 25 = high
% ph_aberr = (pi/2) - 0.5*atan(zMat(:,:,6)./zMat(:,:,4));
z = z(256:end-257,256:end-257);
figure;imagesc(z); colormap jet
%%
aberrated_image = phase + z;
figure, imagesc(aberrated_image); colormap jet;
%%
conj_phase = exp( 1i*( aberrated_image ) );
figure, imagesc(angle(conj_phase)); colormap jet;
%%
tic
[filtered, sigma] = funcDRT(conj_phase,5,2);
toc
phase_1 = real(filtered);
figure, imagesc(phase_1); colormap jet;

%%
True_phase_simu = conj_phase ./ filtered;
figure, imagesc(angle(True_phase_simu)); colormap jet;

%%
phase_unwrap = Unwrap_TIE_DCT_Iter(angle(True_phase_simu));
figure, imagesc(phase_unwrap);colormap jet;
%%
figure, 
plot(phase(:,81),'LineWidth',1.5,'color',[0 0 1]);
axis tight;
hold on
plot(aberrated_image(:,81),'LineWidth',1.5,'color',[1 0 0]);
plot(phase_unwrap(:,81),'LineWidth',1.5,'color',[0 1 0]);
legend({'Phase without aberration';'Phase with aberration'; 'DRT'},'FontSize',11,'FontWeight','bold');
hold off

%%
% Generate phase aberration using Zernike polynomials
% % zernike_order = 4; % Quadratic aberration
% % aberration_mag = 2;
% % F = Zernike_func_Cartesian_Cordinate (x);
% % 
% % z = aberration_mag * ( F(:,:,4) );
% % z = z - min(z(:));
% % % Apply aberration to the image
% % aberrated_image = phase + z; %exp(1i * 2*pi * z);
% % 
% % % Display results
% % subplot(1,2,1);
% % imagesc(phase);
% % title('Original image');
% % axis image;
% % subplot(1,2,2);
% % imagesc(abs(aberrated_image).^2);
% % title('Aberrated image');
% % axis image;
% % 
% % figure, imagesc(aberrated_image); colormap jet;