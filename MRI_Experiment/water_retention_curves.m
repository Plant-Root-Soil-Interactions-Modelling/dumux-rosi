clear all; 
close all; 
fontsize = 16; 

%1= soil, 2= sand
ths = [0.472; 0.39]; 
thr = [0.015; 0.02]; 
alpha = [0.0271; 0.013]; 
n = [1.77; 12.7]; 
m = 1-1./n; 

names = {'soil', 'sand'}; 
h = logspace(0, 7, 1000).'; 
theta = zeros(size(h)); 

for hh = 1:size(thr,1)
    for ii = 1:size(h, 1) 
            theta(ii,hh) = thr(hh)+(ths(hh)-thr(hh))/(1+(alpha(hh)*h(ii))^n(hh))^m(hh); 
    end 
end

figure; 
for hh = 1:size(n,1) 
    tt(hh) = semilogy(theta(:,hh),h,'LineWidth',1.5); hold all; 
    ylim([10^0, 1.5*10^4]);  
    set(gca, 'ytick',logspace(0,4,5)); 
    set(gca, 'xlim', [0 0.5]);
    xlabel('water content (cm^3/cm^3)')
    ylabel('-\psi (cm)')
    set(gca, 'FontSize', fontsize); 
end
legend(tt, names)

