% Name: Matthew Watson
% Date: March 28th, 2019
% Description: Find Cell Objects with Nuclei Mask for MSCs on 35 kPa PAAm
% gels

clear all
close all
clc

%% Load image
nameGFP = '35kPa_rMSC_well1_w2GFP Widefield_s2.TIF';   % Change this for your own images
I = imread(nameGFP);
I = imadjust(I);

I_gauss = imgaussfilt(I,[3 3]);                        % Slight blur for thresholding


%% Binarize the image
bw = imbinarize(I_gauss, .35);                         % 0.35 is a manual threshold
bw2 = imfill(bw, 'holes');
bw3 = imopen(bw2, ones(10,10));
bw4 = bwareaopen(bw3, 10);
bw4_perim = bwperim(bw4);

%% Uncomment figures to check that you have a good threshold
overlay1 = imoverlay(I_gauss, bw4_perim, [1 0 0]);
%figure
%imshow(overlay1)

%% Load in Nuclei image
nameDAPI = '35kPa_rMSC_well1_w1DAPI Widefield_s2.TIF';     % Name of image with nuclei stain
I_nuc = imread(nameDAPI);


%% Binarize Nuclei
nuc_bw = imbinarize(I_nuc, .0025);
nuc_bw2 = imfill(nuc_bw, 'holes');
nuc_bw3 = imopen(nuc_bw2, ones(10,10));
nuc_bw4 = bwareaopen(nuc_bw3, 10);
mask_em = nuc_bw4;
%figure
%imshow(mask_em)
%% Mask nuclei and cell objects
overlay2 = imoverlay(I_gauss, bw4_perim | mask_em, [0.3 1 0.3]);
figure
imshow(overlay2)
I_gauss_c = imcomplement(I_gauss);
I_mod = imimposemin(I_gauss_c, ~bw4 | mask_em);
%% Segment cell bodies

waterL_GFP = watershed(I_mod);
L = label2rgb(waterL_GFP);
figure;
imshowpair(I_gauss, L,'montage')
%% Set some constrictions on what are actual cells


waterL_GFP = double(waterL_GFP);
[height, width] = size(waterL_GFP);

% Find the area of cell objects

%Calculate sizes of elements
for n = 2:max(unique(waterL_GFP))
    i = n - 1;
    count = 0;
    for y = 1:height
        for x = 1:width
            if waterL_GFP(y, x) == n
                count = count + 1;
            end
        end
    end
    cell_size(i) = count;
end

unique(waterL_GFP)
%% Removes elements under certain area size:

maxCellSize = 1000;

for n = 2:max(unique(waterL_GFP))
    i = n - 1;
    count = 0;
    for y = 1:height
        for x = 1:width
            if waterL_GFP(y, x) == n && cell_size(i) < maxCellSize
                waterL_GFP(y,x) = 0;
            end
        end
    end
end

unique(waterL_GFP)
waterL_GFP = imclearborder(waterL_GFP);
unique(waterL_GFP)
% figure
% imagesc(waterL_GFP)

% Find center of mass GFP stain
for n = 2:max(unique(waterL_GFP))
    i = 1;
    xCoord = 0;
    yCoord = 0;


    for y = 1:height
        for x = 1:width
            if waterL_GFP(y, x) == n
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

COM_GFP = double(x_y);

waterL_Marker = insertMarker(waterL_GFP,x_y,'+','color','black','size',10);
waterL_Marker_BW = im2bw(waterL_Marker);
se = strel('disk',2,0);
waterL_Marker_BW = imcomplement(waterL_Marker_BW);
waterL_Marker_BW = imdilate(waterL_Marker_BW,se);
%figure
%imshowpair(waterL_GFP,waterL_Marker_BW,'montage')
I_gray = mat2gray(I);
nameGFP = '35kPa_rMSC_well1_w2GFP Widefield_s2.TIF';
IGFP = imread(nameGFP);
IGFP_gray = mat2gray(IGFP);

centerOverlay = imoverlay(I,waterL_Marker_BW,'red');

I_with_Marker = insertMarker(I,x_y,'+','color','black','size',10);
%figure
%imshow(I_with_Marker)


%figure
%imshow(centerOverlay)

save('GFP_35kPA_COM.mat','COM_GFP','waterL_GFP')



