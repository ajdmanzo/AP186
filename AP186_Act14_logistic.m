clc; clear all; close all;

%data

X = [];
ru = matfile('rr_um.mat'); ru =  ru.rr_um;
gu = matfile('gg_um.mat'); gu = gu.gg_um;
bu = matfile('bb_um.mat'); bu = bu.bb_um;
%ripe
rr = matfile('rr_m.mat'); rr =  rr.rr_m;
gr = matfile('gg_m.mat'); gr = gr.gg_m;
br = matfile('bb_m.mat'); br = br.bb_m;
x1 = [ru;rr];
x2 = [gu;gr];
x3 = [bu;br];
for i=1:length(x1)
    Xx = [1 x1(i) x2(i) x3(i)];
    X = [X;Xx];
end
%% Perceptron 
% 
% %vector x
% X = [];
% %ripe
% x11  = rr;
% x22 = gr;
% x33 = br; 
% % unripe
% xx1 = ru;
% xx2 = gu;
% xx3 = bu;
% 
% x1 = [x11;xx1];
% x2 = [x22; xx2];
% x3 = [x33, xx3];
% 
% for i=1:length(x1)
%     Xi = [1 x1(i) x2(i)];
%     X = [X;Xi];
%     
% end
% 
nu = 0.0001
%initialize weights
% w = rand(length(x1),3); w0 = w(:,1); w1=(:,2);w2=w(:,3);
% w = rand(1,3);
w = [0.3 0.7 0.6 0.0001]; 
w0 = w(1); w1=w(2);w2=w(3); w3 =  w(4);
% w2 = rand(20,3);
a = []; z = [];
d1 = zeros(length(ru),1); %unripe
d2 = ones(length(rr),1); %ripe
d = [d1;d2];
%weight updates for each ii

% figure()
% plot(rr,gr,'y.')
% hold on
% plot(ru,gu,'r*')
% plot(xx,(yy),'g')
% xlabel('r');
% ylabel('eccentricity');
% hold off

% w = [1 9 5]; 
% w0 = w(1); w1=w(2);w2=w(3); %w3=w(4);
% a = []; z = [];  
% d1 = zeros(length(x11)*2,1);
% d2 = ones(length(xx1)*2,1);
% d = [d1;d2];
%weight updates for each ii
for ii=1:length(x1)
   
%     a = (X(ii,:))*transpose(w);
    a = w*transpose(X(ii,:));
    z = sigmoid(a);
    delta_w0 = nu*(d(ii)-z)*(X(ii,1));
    delta_w1 = nu*(d(ii)-z)*(X(ii,2));
    delta_w2 = nu*(d(ii)-z)*(X(ii,3));
    delta_w3 = nu*(d(ii)-z)*(X(ii,4));
    w = [w0+delta_w0 w1+delta_w1 w2+delta_w2 w3+delta_w3];
    disp(w)
end

% %===============================TEST=======================================
test_img = im2double(imread('um6.jpg')); 
[r, c, num] = size(test_img);
im = imcrop(test_img, [(c/2)-90 (r/2)-90 100 100]);
figure(); imshow(im);

R = im(:,:,1);
G = im(:,:,2); 
B = im(:,:,3); 
I = R+G+B;
I(find(I==0))=100000;
%normalized chromaticity coordinates
r = R./I;
g = G./I;
b = B./I;
% 
%initial weights - weights obtained from training 
x_test = [1 mean2(r) mean2(g) mean2(b)]; % mean2(b)]
%a_test = x_test*transpose(w); 
a_test = w*transpose(x_test); 
z_test = sigmoid(a_test); 

