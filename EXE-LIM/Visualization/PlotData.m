
% Visualize the data in the EXE-LIM network
% Guoshi Li (5/24/2020)


clc;
close all;
clear all;

load EC;  % Estimated (optimized) EC for 98 NC and 96 MDD subjects
load SC;  % Structural connectivity with reduced models

NR = 9;   % Number of ROIs
NP = 71;  % Number of estimated parameters
NM = 2;   % Number of networks

Nb = [4 5]; % Number of ROIs in each network

n1 = Count_NC;  % Number of NC subjects
n2 = Count_MDD; % Number of MDD subjects

%==========================================================================
%              Individual Connection Analysis
%==========================================================================

% Average weight
Wm_NC = mean(W_NC);
Wm_MD = mean(W_MDD);

STD_NC = std(W_NC);
STD_MD = std(W_MDD);

SE_NC = STD_NC/sqrt(n1);
SE_MD = STD_MD/sqrt(n2);

% Put NC and MDD together
Wm_NC_MD  = [Wm_NC;   Wm_MD];
STD_NC_MD = [STD_NC; STD_MD];
SE_NC_MD  = [SE_NC; SE_MD];

% Divide the parameters into four sets
Wm_EE = Wm_NC_MD(:, 1:NR);
Wm_IE = Wm_NC_MD(:, NR+1:2*NR);
Wm_GC = Wm_NC_MD(:, 2*NR+1:NP-1);
Wm_PI = Wm_NC_MD(:, NP);

Wm_GC1 = Wm_GC(:, 1:N2/2);
Wm_GC2 = Wm_GC(:, N2/2+1:end);

% Standard Error
EE_ER = SE_NC_MD(:, 1:NR);
IE_ER = SE_NC_MD(:, NR+1:2*NR);
GC_ER = SE_NC_MD(:, 2*NR+1:NP-1);

GC_ER_NC = GC_ER(1,:);
GC_ER_MD = GC_ER(2,:);


% Reconstruct the inter-regional NR*NR matrix

m = 1;

GC1 = Wm_GC(1,:);
GC2 = Wm_GC(2,:);

MGC_NC = zeros(NR, NR);
MGC_MD = zeros(NR, NR);

MSE_NC = zeros(NR, NR);
MSE_MD = zeros(NR, NR);

for i=1:NR
    for j=1:NR
       
       if (MAP2(i,j)==1)
           MGC_NC(i,j) = GC1(m);
           MGC_MD(i,j) = GC2(m);
           
           MSE_NC(i,j) = GC_ER_NC(m);
           MSE_MD(i,j) = GC_ER_MD(m);
                      
           m = m+1;
       end      
        
    end  
end

% Weights
VGC_NC = MGC_NC(:);
VGC_MD = MGC_MD(:);

VGC_NC = VGC_NC(VGC_NC~=0);
VGC_MD = VGC_MD(VGC_MD~=0);

VGC_NC_MD = [VGC_NC VGC_MD];

VGC1 = VGC_NC_MD(1:N2/2, :);
VGC2 = VGC_NC_MD(N2/2+1:end, :);

% Standard errors
VSE_NC = MSE_NC(:);
VSE_MD = MSE_MD(:);

VSE_NC = VSE_NC(VSE_NC~=0);
VSE_MD = VSE_MD(VSE_MD~=0);

VSE_NC_MD = [VSE_NC VSE_MD];

VSE1 = VSE_NC_MD(1:N2/2, :);
VSE2 = VSE_NC_MD(N2/2+1:end, :);

%==========================================================================
%              Network Analysis
%==========================================================================

% Intra-regional weights
WNee = [mean(Wm_EE(:,1:4), 2) mean(Wm_EE(:,5:end), 2)  ];
WNie = [mean(Wm_IE(:,1:4), 2) mean(Wm_IE(:,5:end), 2)  ];

Wee_NC = W_NC(:,1:NR);
Wie_NC = W_NC(:,NR+1:2*NR);
Wgc_NC = W_NC(:,2*NR+1:NP-1);
SPI_NC = W_NC(:,NP);

Wee_MDD = W_MDD(:,1:NR);
Wie_MDD = W_MDD(:,NR+1:2*NR);
Wgc_MDD = W_MDD(:,2*NR+1:NP-1);
SPI_MDD = W_MDD(:,NP);

WeeN_NC = [mean( Wee_NC(:,1:4), 2)  mean(Wee_NC(:, 5:end), 2)  ];
WieN_NC = [mean( Wie_NC(:,1:4), 2)  mean(Wie_NC(:, 5:end), 2)  ];

WeeN_MDD = [mean(Wee_MDD(:, 1:4), 2)   mean(Wee_MDD(:, 5:end), 2)  ];
WieN_MDD = [mean(Wie_MDD(:, 1:4), 2)   mean(Wie_MDD(:, 5:end), 2)  ];

WN_NC  = [WeeN_NC  WieN_NC  SPI_NC];
WN_MDD = [WeeN_MDD WieN_MDD SPI_MDD];

% Standard errors
STD_NET_NC = std(WN_NC);
STD_NET_MD = std(WN_MDD);

SE_NET_NC = STD_NET_NC/sqrt(n1);
SE_NET_MD = STD_NET_MD/sqrt(n2);

SE_NET_NC_MD = [SE_NET_NC; SE_NET_MD];

SE_EE_NET = SE_NET_NC_MD(:,1:2);
SE_IE_NET = SE_NET_NC_MD(:,3:4);
SE_PI = SE_NET_NC_MD(:,end);


% Inter-regional weights
GCN_NC = zeros(NM, NM);
GCN_EXC_NC = zeros(NM, NM);
GCN_INH_NC = zeros(NM, NM);

GCN_MD = zeros(NM, NM);
GCN_EXC_MD = zeros(NM, NM);
GCN_INH_MD = zeros(NM, NM);

k1 = 1;
l1 = 1;

for i = 1:NM
    
    l1 = 1;
    
    for j = 1:NM
   
      k2 = k1+Nb(i)-1;  
      l2 = l1+Nb(j)-1;
      
      BLOCK = MGC_NC(k1:k2, l1:l2);
      [a, b] = find(BLOCK~=0);
      
      GCN_NC(i, j) = sum(sum(BLOCK))/length(a);
            
      PBLOCK = BLOCK(:);
      
      ID1 = find(PBLOCK>0);
      ID2 = find(PBLOCK<0);
      
      if ( length(ID1)~=0 ) 
         TEMP1 = PBLOCK(ID1);
         GCN_EXC_NC(i, j) = sum(TEMP1)/length(ID1);
      end
      
      if ( length(ID2)~=0 )  
         TEMP2 = PBLOCK(ID2);
         GCN_INH_NC(i, j) = sum(TEMP2)/length(ID2);
      end      
        
      l1 = l1+Nb(j);
      
    end
    
    k1 = k1+Nb(i);
end


k1 = 1;
l1 = 1;

for i = 1:NM
    
    l1 = 1;
    
    for j = 1:NM
   
      k2 = k1+Nb(i)-1;  
      l2 = l1+Nb(j)-1;
      
      BLOCK = MGC_MD(k1:k2, l1:l2);
      [a, b] = find(BLOCK~=0);
      
      GCN_MD(i, j) = sum(sum(BLOCK))/length(a);
            
      PBLOCK = BLOCK(:);
      
      ID1 = find(PBLOCK>0);
      ID2 = find(PBLOCK<0);
      
      if ( length(ID1)~=0 ) 
         TEMP1 = PBLOCK(ID1);
         GCN_EXC_MD(i, j) = sum(TEMP1)/length(ID1);
      end
      
      if ( length(ID2)~=0 )  
         TEMP2 = PBLOCK(ID2);
         GCN_INH_MD(i, j) = sum(TEMP2)/length(ID2);
      end      
        
      l1 = l1+Nb(j);
      
    end
    
    k1 = k1+Nb(i);
end



%==========================================================================

FVm_NC  = mean(FVAL_NC);
FVm_MD = mean(FVAL_MDD);
FVm   =  [FVm_NC FVm_MD];

SE_FV_NC = std(FVAL_NC)/sqrt(n1);
SE_FV_MD = std(FVAL_MDD)/sqrt(n2);
SE_FV = [SE_FV_NC SE_FV_MD];

SD_FV = [std(FVAL_NC) std(FVAL_MDD)];

%==========================================================================
%              Ploting
%==========================================================================

Wm_EE  = Wm_EE';
Wm_IE  = Wm_IE';
Wm_GC  = Wm_GC';
Wm_GC1 = Wm_GC1';
Wm_GC2 = Wm_GC2';

Wm_PI = Wm_PI';

WNee = WNee';
WNie = WNie';

FLIP_VGC_NC_MD = flip(VGC_NC_MD);


% STD
EE_ER = EE_ER';
IE_ER = IE_ER';
GC_ER = GC_ER';

GC1_ER = GC_ER(1:N2/2,:);
GC2_ER = GC_ER(N2/2+1:end,:);

SE_EE_NET = SE_EE_NET';
SE_IE_NET = SE_IE_NET';
SE_PI = SE_PI';

FLIP_VSE_NC_MD = flip(VSE_NC_MD);



ROI_LABEL = {'L.dlPFC', 'R.dlPFC', 'L.SPC', 'R.SPC'...
               'Thal', 'L.Amyg', 'R.Amyg', 'L.HPC', 'R.HPC'};
           
ROI_LABEL_NET = {'EXE', 'LIM'};             
           
           
INTER_LABEL1 = {'L.dlPFC-R.SPC', 'L.dlPFC-Thal', 'L.dlPFC-L.HPC', 'R.dlPFC-L.SPC', 'R.dlPFC-R.SPC',...
                 'R.dlPFC-Thal',  'R.dlPFC-R.Amyg', 'R.dlPFC-R.HPC', 'L.SPC-R.dlPFC', 'L.SPC-R.SPC',...
                 'L.SPC-Thal',   'L.SPC-R.Amyg', 'L.SPC-L.HPC', 'L.SPC-R.HPC', 'R.SPC-L.dlPFC',...
                 'R.SPC-R.dlPFC','R.SPC-L.SPC', 'R.SPC-Thal', 'R.SPC-R.Amyg', 'R.SPC-L.HPC',...
                 'R.SPC-R.HPC',  'Thal-L.dlPFC',  'Thal-R.dlPFC', 'Thal-L.SPC', 'Thal-R.SPC',...
                 'Thal-L.Amyg'};    
             
INTER_LABEL2 = { 'Thal-R.Amyg', 'Thal-L.HPC', 'Thal-R.HPC', 'L.Amyg-Thal', 'L.Amyg-L.HPC',...
                 'L.Amyg-R.HPC', 'R.Amyg-R.dlPFC', 'R.Amyg-L.SPC', 'R.Amyg-R.SPC', 'R.Amyg-Thal',...
                 'R.Amyg-L.HPC', 'R.Amyg-R.HPC', 'L.HPC-L.dlPFC', 'L.HPC-L.SPC', 'L.HPC-R.SPC',...
                 'L.HPC-Thal', 'L.HPC-L.Amyg', 'L.HPC-R.Amyg', 'L.HPC-R.HPC', 'R.HPC-R.dlPFC',...
                 'R.HPC-L.SPC', 'R.HPC-R.SPC', 'R.HPC-Thal', 'R.HPC-L.Amyg','R.HPC-R.Amyg','R.HPC-L.HPC'};        
             
INTER_LABEL = [INTER_LABEL1 INTER_LABEL2]; 

FLIP_INTER_LABEL = flip(INTER_LABEL);
           

ROI_LABEL_NET = {'EXE', 'LIM'};          
                    
          
Xtick = 1:NR;        

figure;
bar(Xtick, Wm_EE);
hold on;
ngroups = size(Wm_EE, 1);
nbars = size(Wm_EE, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, Wm_EE(:,i), EE_ER(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 12);
title('Recurrent Excitation', 'FontSize', 16);
set(gca, 'XTick', Xtick, 'XTickLabel', ROI_LABEL);
ylabel('Weight', 'FontSize', 14);
xtickangle(45);
axis([0 NR+1 2.5 3.6]);
set(gca, 'YTick', [2.5:0.5:4]); 
box('off');
legend('NC', 'MDD');



figure;
bar(Xtick, Wm_IE);
hold on;
ngroups = size(Wm_IE, 1);
nbars = size(Wm_IE, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, Wm_IE(:,i), IE_ER(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 12);
title('Recurrent Inhibition', 'FontSize', 16);
set(gca, 'XTick', Xtick, 'XTickLabel', ROI_LABEL);
ylabel('Weight', 'FontSize', 14);
xtickangle(45);
axis([0 NR+1 2.0 3.6]);
set(gca, 'YTick', [2:0.5:4]); 
box('off');
legend('NC', 'MDD');



%==========================================================================
Xgc = 1:N2/2;


%==========================================================================
% Spontaneous input
Xn = 1:2;

figure;
bar(Xn, Wm_PI);
hold on;
er = errorbar(Xn, Wm_PI, SE_PI, '.');
er.Color = [0 0 0];
er.LineWidth = 1;
er.LineStyle = 'none'; 
hold off
set(gca, 'FontSize', 12);
title('Spont. Input', 'FontSize', 16);
set(gca, 'XTick', [1, 2], 'XTickLabel', {'NC', 'MDD'});
axis([0 3 0.2 0.4]);
box('off');
% legend('NC', 'MDD');


%==========================================================================
% Objective function

figure;
bar(Xn, FVm);
hold on;
er = errorbar(Xn, FVm, SD_FV, '.');
er.Color = [0 0 0];
er.LineWidth = 1;
er.LineStyle = 'none'; 
hold off
set(gca, 'FontSize', 14);
title('EXE-LIM Network');
ylabel('Fitness Value'); 
set(gca, 'XTick', [1, 2], 'XTickLabel', {'NC', 'MDD'});
axis([0 3 0.0 1.0]);
box('off');



%==========================================================================
figure;
bar(Xn, WNee);
hold on;
ngroups = size(WNee, 1);
nbars = size(WNee, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, WNee(:,i), SE_EE_NET(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 12);
title('Recurrent Excitation', 'FontSize', 16);
set(gca, 'XTick', Xn, 'XTickLabel', ROI_LABEL_NET);
ylabel('Weight', 'FontSize', 14);
axis([0.5 2.5 2.5 3.5]);
set(gca, 'YTick', [2:0.5:4]); 
box('off');
legend('NC', 'MDD');


figure;
bar(Xn, WNie);
hold on;
ngroups = size(WNie, 1);
nbars = size(WNie, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, WNie(:,i), SE_IE_NET(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 12);
title('Recurrent Inhibition', 'FontSize', 16);
set(gca, 'XTick', Xn, 'XTickLabel', ROI_LABEL_NET);
ylabel('Weight', 'FontSize', 14);
axis([0.5 2.5 2.5 3.2]);
set(gca, 'YTick', [2:0.5:4]); 
box('off');
legend('NC', 'MDD');



Xgch = 1:N2;

figure;
barh(Xgch, FLIP_VGC_NC_MD);
hold on;
ngroups = size(FLIP_VGC_NC_MD, 1);
nbars = size(FLIP_VGC_NC_MD, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(FLIP_VGC_NC_MD(:,i), x, FLIP_VSE_NC_MD(:,i), '.','horizontal' );
    er.Color = [0 0 0];
    er.LineWidth = 0.5;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 10);
title('Inter-regional Coupling', 'FontSize', 14);
set(gca, 'YTick', Xgch, 'YTickLabel', FLIP_INTER_LABEL );
xlabel('Weight', 'FontSize', 14);
% xtickangle(90);
% axis([0 ngroups+1 -1.0 1.5]);
% set(gca, 'YTick', [-1:0.5:1.5]); 
box('off');
legend('NC', 'MDD');


%==========================================================================

CLIM = [-0.5 1.5];

figure;
h1 = heatmap(ROI_LABEL, ROI_LABEL, round(MGC_NC, 3), 'colormap', jet, 'ColorLimits', CLIM );
% set(gca, 'Fontsize', 12);
h1.FontSize = 12;
h1.CellLabelFormat = '%.2g';
title('NC');

figure;
h2 = heatmap(ROI_LABEL, ROI_LABEL, round(MGC_MD, 3), 'colormap', jet, 'ColorLimits', CLIM );
% set(gca, 'Fontsize', 12);
h2.FontSize = 12;
h2.CellLabelFormat = '%0.2g';
title('MDD');


%==========================================================================









