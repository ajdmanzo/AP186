clc; clear all; close all;


%orange
eo = matfile('eccen_o.mat'); eo = eo.eccen_m; %load matrices
ro = matfile('rr_o.mat'); ro = ro.rr_m;
go = matfile('gg_o.mat'); go = go.gg_m;
%mango
em = matfile('eccen_m.mat'); em = em.eccen_m;
rm = matfile('rr_m.mat'); rm = rm.rr_m;
gm = matfile('gg_m.mat'); gm = gm.gg_m;
%banana
eb = matfile('eccen_b.mat'); eb = eb.eccen_m;
rb = matfile('rr_b.mat'); rb = rb.rr_m;
gb = matfile('gg_b.mat'); gb = gb.gg_m;
eb1 = matfile('eccen_b1.mat'); eb1 = eb1.eccen_m;
rb1 = matfile('rr_b1.mat'); rb1 = rb1.rr_m;
gb1 = matfile('gg_b1.mat'); gb1 = gb1.gg_m;

gb = [gb;gb1];
rb = [rb;rb1];
eb = [eb;eb1];

%Number of fruit classes
M = 3;
%mean of each class
mean1 = mean(em);
mean2 = mean(eb); 
mean3 = mean(eo);
%covariance of each class
cov1 = cov(em); %mango
cov2 = cov(eb); %banana
cov3 = cov(eo); %orange

%initialize probabilities
P1 = 1/M; P2 = P1; P3 = P1;
d = 2;
[X,Y] = meshgrid(em,rm);
%pdf 
% function [pdf] = pdf(x,mean,cov)
%     pdf = (1/((2*pi)^(d/2))*sqrt(cov))*exp((-1/2)*(transpose(x-mean))*(x-mean)/cov);
% end

pdf = @(x,mean,cov) (1/((2*pi)^(d/2))*sqrt(cov))*exp((-1/2)*(transpose(x-mean))*(x-mean)/cov);
for i = 1:22
    pdf1 =  arrayfun(@(x) pdf(X,mean1,cov1));
    
end
figure()
plot(em,pdf1,'b.')


figure()
contour(X,Y,pdf1);


