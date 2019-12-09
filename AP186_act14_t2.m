clc; clear all; close all;

%data
r_ripe = table2array(readtable('C:\Users\Asus\Documents\Applied Physics 186\act 14\results\r1_ripe.csv'));
g_ripe = table2array(readtable('C:\Users\Asus\Documents\Applied Physics 186\act 14\results\g1_ripe.csv'));
b_ripe = table2array(readtable('C:\Users\Asus\Documents\Applied Physics 186\act 14\results\b1_ripe.csv'));

rr = r_ripe(:,1);
gr = g_ripe(:,1);
br = b_ripe(:,1); 

r_unripe = table2array(readtable('C:\Users\Asus\Documents\Applied Physics 186\act 14\results\r1_unripe.csv'));
g_unripe = table2array(readtable('C:\Users\Asus\Documents\Applied Physics 186\act 14\results\g1_unripe.csv'));
b_unripe = table2array(readtable('C:\Users\Asus\Documents\Applied Physics 186\act 14\results\b1_unripe.csv'));

ru = r_unripe(:,1);
gu = g_unripe(:,1);
bu = b_unripe(:,1); 

X = [];
ru = 
gu = 
bu = 
%ripe
rr = 
gr = 
br = 

x1 = [ru;rr];
x2 = [gu;gr];
x3 = [bu;br];
for i=1:length(x1)
    Xx = [1 x1(i) x2(i) x3(i)];
    X = [X;Xm];
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
% nu = 0.12
% %initialize weights
% % w = rand(length(x1),3); w0 = w(:,1); w1=(:,2);w2=w(:,3);
% % w = rand(1,3);
% w = [0.3 0.7 0.6]; 
% w0 = w(1); w1=w(2);w2=w(3);
% % w2 = rand(20,3);
% a = []; z = [];
% d1 = ones(length(x11),1);
% d2 = -1.*(ones(length(xx1),1));
% d = [d1;d2];
% %weight updates for each ii
% for ii=1:length(x1)
% %     d1 = ones(20,1); 
% %     d2 = -1.*(ones(22,1));
%     d = [d1;d2];
%     a = (X(ii,:))*transpose(w);
%     z = g(a);
%     delta_w0 = nu*(d(ii)-z)*(X(ii,1));
%     delta_w1 = nu*(d(ii)-z)*(X(ii,2));
%     delta_w2 = nu*(d(ii)-z)*(X(ii,3));
%     w = [w0+delta_w0 w1+delta_w1 w2+delta_w2];
% %     z = [z;z1];
%     
% end
% 
% m = w(2)/w(3);
% b = -w(1)/w(3);
% xx = [0:0.01:1];
% yy = m*xx +b;
% 
% figure()
% plot(rr,gr,'y.')
% hold on
% plot(ru,gu,'r*')
% plot(xx,(yy),'g')
% xlabel('r');
% ylabel('eccentricity');
% hold off

%% Logistic Regression

%===============================TRAINING=======================================
%vector x
X = [];
%unripe
x11  = ru;
x22 = gu;
x33 = bu; 
%ripe
xx1 = rr;
xx2 = gr;
xx3 = br;

x1 = [x11; x11; xx1; xx1];
x2 = [x22; x22; xx2; xx2];
x3 = [x33; x33; xx3; xx3];

for i=1:length(x1)
    Xi = [1 x1(i) x2(i)]; % x3(i)];
    X = [X;Xi];
    
end

nu = 0.1
%initialize weights
w = [1 9 5]; 
w0 = w(1); w1=w(2);w2=w(3); %w3=w(4);
a = []; z = [];  
d1 = zeros(length(x11)*2,1);
d2 = ones(length(xx1)*2,1);
d = [d1;d2];
%weight updates for each ii
for ii=1:length(x1)
    d1 = ones(20,1); 
    d2 = -1.*(ones(22,1));
    d = [d1;d2];
    a = (X(ii,:))*transpose(w);
    z = sigmoid(a);
    delta_w0 = nu*(d(ii)-z)*(X(ii,1));
    delta_w1 = nu*(d(ii)-z)*(X(ii,2));
    delta_w2 = nu*(d(ii)-z)*(X(ii,3));
    delta_w3 = nu*(d(ii)-z)*(X(ii,4));
    w = [w0+delta_w0 w1+delta_w1 w2+delta_w2]; %w3+delta_w3];
    z = [z;z1];
    
end

%===============================TEST=======================================
test_img = im2double(imread('C:\Users\Asus\Documents\Applied Physics 186\act 14\test\mango4.png')); 
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

%initial weights - weights obtained from training 
x_test = [1 mean2(r) mean2(g)]; % mean2(b)]
%a_test = x_test*transpose(w); 
a_test = w*transpose(x_test); 
z_test = sigmoid(a_test); 

