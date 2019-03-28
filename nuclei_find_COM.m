clear all
close all
clc

% Name: Matthew Watson
% Date: March 28th, 2019
%
% Description: Finds the center of masses of viable nuclei.

%% Load Nuclei Images
nameDAPI = '35kPa_rMSC_well1_w1DAPI Widefield_s2.TIF';
I_nuc = imread(nameDAPI);
bw = imbinarize(I_nuc, .0025);
bw2 = imfill(bw, 'holes');
bw3 = imopen(bw2, ones(10,10));
bw4 = bwareaopen(bw3, 10);
bw4_perim = bwperim(bw4);


imshow(bw4)
figure
imshow(I_nuc,[])
%% Marker-based watershed transformation
mask_em = imextendedmax(I_nuc,20);

% Clean up and overlay to check
mask_em = imclose(mask_em, ones(5,5));
mask_em = imfill(mask_em, 'holes');
mask_em = bwareaopen(mask_em, 30);


% Implement marker and then take watershed transform
I_complement = imcomplement(I_nuc);
I_mod = imimposemin(I_complement, ~bw4 | mask_em);
imshow(I_mod);
%%
L = watershed(I_mod);
figure
imshow(label2rgb(L));

waterL_Nuclei = double(L);
[height, width] = size(waterL_Nuclei);
unique(waterL_Nuclei)
%% Calculate sizes of elements
for n = 2:max(unique(waterL_Nuclei))
    i = n - 1;
    count = 0;
    for y = 1:height
        for x = 1:width
            if waterL_Nuclei(y, x) == n
                count = count + 1;
            end
        end
    end
    cell_size(i) = count;
end

%% Removes elements under certain size:

sizeConstraint = 100; % change 100 for size filter
for n = 2:max(unique(waterL_Nuclei))
    i = n - 1;
    count = 0;
    for y = 1:height
        for x = 1:width
            if waterL_Nuclei(y, x) == n && cell_size(i) < sizeConstraint 
                waterL_Nuclei(y,x) = 0;
            end
        end
    end
end
unique(waterL_Nuclei)
%% Get rid of cells on border
waterL_Nuclei = imclearborder(waterL_Nuclei);
unique(waterL_Nuclei)
%% Find center of mass nuclei
for n = 2:max(unique(waterL_Nuclei))
    i = 1;
    xCoord = 1;
    yCoord = 1;


    for y = 1:height
        for x = 1:width
            if waterL_Nuclei(y, x) == n
                    xCoord(i) = x;
                    yCoord(i) = y;
                    i = i + 1;
            end
        end
    end
    x_sum = sum(xCoord);
    y_sum = sum(yCoord);
    massTotal = length(xCoord);

    xCom = int16(x_sum/massTotal);
    yCom = int16(y_sum/massTotal);
    x_y(n - 1, 1) = xCom;
    x_y(n - 1, 2) = yCom;

end

COM_Nuclei = double(x_y);

% Saves .mat file
save('Nuclei_35kPa_COM.mat','COM_Nuclei','waterL_Nuclei');