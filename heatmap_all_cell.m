clear all
close all
clc

% Name: Matthew Watson
% Date: March 28th, 2019
% 
% Description: Plots a surface of the cell

arrayFactor = 1;            % Do not change this

stiffness = 35;             % Change this for different gel stiffness example

% does a file exist for n
for n = 2:8
    filename = ['all_cell_',num2str(n),'_',num2str(stiffness),'kPa_protein.mat']; % change file name for different stains
    contained = exist(filename,'file');

    if contained == 2
        cellData = load(filename);

        % make a data table of x, y, and intensity.
        x = cos(cellData.phase*pi/180).*cellData.comPtDist;
        y = sin(cellData.phase*pi/180).*cellData.comPtDist;
        z = cellData.norm_intensity;
        
        indexArrayStart = arrayFactor*2000-1999;
        indexArrayEnd = arrayFactor*2000;
        index = indexArrayStart:indexArrayEnd;

        data(index,1) = x;
        data(index,2) = y;
        data(index,3) = z; 

        data_norm(index,1) = cellData.normal_r;
        data_norm(index,2) = cellData.phase;
        data_norm(index,3) = cellData.norm_intensity/(max(cellData.norm_intensity));

        % Increase arrayFactor
        arrayFactor = arrayFactor+1;
    end
end


    data = double(data);
    data_norm = double(data_norm);
    % Make phase data go from 0 - 360;

    for i = 1:length(data_norm(:,2))
        if data_norm(i,2) < 0
            data_norm(i,2) = 360 + data_norm(i,2);
        end
    end

    radial_bins = 10;
    phase_bins = 12;

    radial_bin_size = 1/radial_bins;
    phase_bin_size = 360/phase_bins;

    bin_number = 1;

    counts = 0;
    for i = 1:radial_bins
        for j = 1:phase_bins
            element_number = 1;
            element = [];
            for k = 1:length(data_norm)
                % binning
                if i == 1 && j == 1
                      if (data_norm(k,1) >= (i-1)*radial_bin_size && data_norm(k,1) <= i*radial_bin_size &&...
                        data_norm(k,2) >= (j-1)*phase_bin_size && data_norm(k,2)<= j*phase_bin_size)

                        % r bin min value
                        element(element_number,1) = (i-1)*radial_bin_size;
                        % r bin max value
                        element(element_number,2) = (i)*radial_bin_size;
                        % phase bin min value
                        element(element_number,3) = (j-1)*phase_bin_size;
                        % phase bin max value
                        element(element_number,4) = (j)*phase_bin_size;
                        % bin intensity value
                        element(element_number,5) = data_norm(k,3);

                        % increase element array index
                        element_number = element_number + 1;

                        % increase number of points collected
                        counts = counts + 1;
                      end
                elseif i == 1 && j ~= 1
                    if (data_norm(k,1) >= (i-1)*radial_bin_size && data_norm(k,1) <= i*radial_bin_size &&...
                        data_norm(k,2) > (j-1)*phase_bin_size && data_norm(k,2)<= j*phase_bin_size)

                        % r bin min value
                        element(element_number,1) = (i-1)*radial_bin_size;
                        % r bin max value
                        element(element_number,2) = (i)*radial_bin_size;
                        % phase bin min value
                        element(element_number,3) = (j-1)*phase_bin_size;
                        % phase bin max value
                        element(element_number,4) = (j)*phase_bin_size;
                        % bin intensity value
                        element(element_number,5) = data_norm(k,3);

                        % increase element array index
                        element_number = element_number + 1;

                        % increase number of points collected
                        counts = counts + 1;
                    end
                elseif j == 1 && i ~= 1
                    if (data_norm(k,1) > (i-1)*radial_bin_size && data_norm(k,1) <= i*radial_bin_size &&...
                        data_norm(k,2) >= (j-1)*phase_bin_size && data_norm(k,2)<= j*phase_bin_size)

                        % r bin min value
                        element(element_number,1) = (i-1)*radial_bin_size;
                        % r bin max value
                        element(element_number,2) = (i)*radial_bin_size;
                        % phase bin min value
                        element(element_number,3) = (j-1)*phase_bin_size;
                        % phase bin max value
                        element(element_number,4) = (j)*phase_bin_size;
                        % bin intensity value
                        element(element_number,5) = data_norm(k,3);

                        % increase element array index
                        element_number = element_number + 1;

                        % increase number of points collected
                        counts = counts + 1;
                    end
                else
                    if (data_norm(k,1) > (i-1)*radial_bin_size && data_norm(k,1) <= i*radial_bin_size &&...
                        data_norm(k,2) > (j-1)*phase_bin_size && data_norm(k,2)<= j*phase_bin_size)

                        % r bin min value
                        element(element_number,1) = (i-1)*radial_bin_size;
                        % r bin max value
                        element(element_number,2) = (i)*radial_bin_size;
                        % phase bin min value
                        element(element_number,3) = (j-1)*phase_bin_size;
                        % phase bin max value
                        element(element_number,4) = (j)*phase_bin_size;
                        % bin intensity value
                        element(element_number,5) = data_norm(k,3);

                        % increase element array index
                        element_number = element_number + 1;

                        % increase number of points collected
                        counts = counts + 1;
                    end
                end
            end

            bin_data{bin_number} = element;
            bin_data_mean(bin_number,:) = mean(bin_data{bin_number});
            bin_data_std(bin_number,:) = std(bin_data{bin_number});
            bin_number = bin_number + 1;
        end

    end

    theta = bin_data_mean(:,4) - (bin_data_mean(:,4) - bin_data_mean(:,3))/2;
    r = bin_data_mean(:,2);
    x = r.*cos(theta);
    y = r.*sin(theta);

    x_y = [x,y];
    R = 0:1/1000:1;
    THETA = 0:360/1000:360;
    THETA = THETA*pi/180;

    X = -2:2/1000:2;
    Y = -2:2/1000:2;
    Z = NaN(length(X),length(Y));



%% Generate Heat Map of Mean values
for k = 1:length(bin_data_mean)
    for i = 1:length(X)
        for j = 1:length(Y)
            current_angle = atan2(Y(j),X(i))*180/pi;
            current_r = sqrt((X(i)^2 + Y(j)^2));
            if current_angle < 0
                current_angle = current_angle + 360;
            end           
            
            if current_r == 0
                
                Z(j,i) = nanmean(Z(j,i) + bin_data_mean(k,5));
                
            elseif current_r > bin_data_mean(k,1) && current_r <= bin_data_mean(k,2) && ...
                    current_angle > bin_data_mean(k,3) && current_angle <= bin_data_mean(k,4)

                Z(j,i) = bin_data_mean(k,5);

            else
                
                Z(j,i) = Z(j,i);
                
            end
                        
        end
    end
end
%%
figure(1)
h = polar(X,Y);
hold on
set(gca,'fontsize',24)
b = imagesc(X,Y,Z);
set(b,'AlphaData',~isnan(Z))
set(h,'Visible','off')
axis([-2 2 -2 2.5])
grid on
colormap(jet)
cb = colorbar;
cb.Location = 'SouthOutside';
cb.Label.String = 'Normalized Intensity';
t = title([num2str(stiffness),' kPa YAP1 Localization']);
t.Position = [0 2.5 0];



%% Generate Heat Map of STD Values

Z2 = NaN(length(X),length(Y));
for k = 1:length(bin_data_mean)
    for i = 1:length(X)
        for j = 1:length(Y)
            current_angle = atan2(Y(j),X(i))*180/pi;
            current_r = sqrt((X(i)^2 + Y(j)^2));
            if current_angle < 0
                current_angle = current_angle + 360;
            end           
            
            if current_r == 0
                
                Z2(j,i) = nanmean(Z2(j,i) + bin_data_std(k,5));
                
            elseif current_r > bin_data_mean(k,1) && current_r <= bin_data_mean(k,2) && ...
                    current_angle > bin_data_mean(k,3) && current_angle <= bin_data_mean(k,4)

                Z2(j,i) = bin_data_std(k,5);

            else
                
                Z2(j,i) = Z2(j,i);
                
            end
                        
        end
    end
end

%%
figure(2)
h = polar(X,Y);
set(gca,'fontsize',24)
hold on
b = imagesc(X,Y,Z2);
set(b,'AlphaData',~isnan(Z2))
set(h,'Visible','off')
axis([-2 2 -2 2.5])
grid on
colormap(jet)
cb = colorbar;
cb.Location = 'SouthOutside';
cb.Label.String = 'Normalized Intensity';
t = title(['Cell ',num2str(n),' Protein Localization STD']);
t.Position = [0 2.5 0];

