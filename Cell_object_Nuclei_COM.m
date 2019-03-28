clear all
close all
clc

% Name: Matthew Watson
% Date: March 28th, 2019
%
% Description: Replaces GFP COM with Nuclei COM_Data

%% Load COM Data Files

load('Nuclei_35kPa_COM.mat')
load('GFP_35kPa_COM.mat')
unique_GFP = unique(waterL_GFP);
unique_Nuclei = unique(waterL_Nuclei);

counter = 1;
%COM_GFP_New = COM_GFP;



for n = 2:max(unique(waterL_GFP))
    
    if ismember(n, unique(waterL_GFP)) == 1
        for i = 1:max(unique_Nuclei)-1
            
            Nuc_COM_x = COM_Nuclei(i,1);
            Nuc_COM_y = COM_Nuclei(i,2);
            waterL_GFP(Nuc_COM_y, Nuc_COM_x)
            
            if waterL_GFP(Nuc_COM_y, Nuc_COM_x) == n
                COM_GFP_New(counter,1) = Nuc_COM_x;
                COM_GFP_New(counter,2) = Nuc_COM_y;
                cellNumber(counter) = n;
            end
        end
        counter = counter + 1;
    end
    
end
cellNumber = cellNumber';
dataTable = table(cellNumber,COM_GFP_New,'VariableNames',{'Cell', 'COM'});
save('GFP_with_Nuclei_35kPa_COM.mat','dataTable','waterL_GFP')

%%
nameGFP = '35kPa_rMSC_well1_w2GFP Widefield_s2.TIF';
I = imread(nameGFP);
I = imadjust(I);

I_marker_old = insertMarker(I,COM_GFP,'+','color','black','size',10);
I_marker_new = insertMarker(I_marker_old,COM_GFP_New,'+','color','red','size',10);

figure
imshowpair(I_marker_old,I_marker_new,'montage')