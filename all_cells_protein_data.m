clear all
close all
clc

% Author: Matthew Watson
% Date: March 15th, 2018
% 
% Description: Secondary antibody intensity qualification for all cells.

%

%% Load Data

load('GFP_with_Nuclei_35kPa_COM.mat')                                                     % Loads the Cell Object COM Data
Stain_image = imread('35kPa_rMSC_well1_w3TRITC Widefield_s2.TIF');                        % Reads in the TRITC channel image
Stain_image = double(Stain_image);                                                        % MATLAB likes to manipulate doubles so convert it to a double
[height, width] = size(Stain_image);                                                      % Determines the height and the widht of the image

%% Get Rid of zeros in COM data
counter = 1;                                                                          % Start indexing at 1

for i = 1:length(dataTable.Cell)                                                      % Go through each center of mass data
    
    if dataTable.COM(i,1) ~= 0 && dataTable.COM(i,2) ~= 0                             % If both x and y center of masses are not equal to zero, then do stuff
        
        COM_GFP_new(counter,1) = dataTable.COM(i,1);                                  % new indexed COM x is non zero COM
        COM_GFP_new(counter,2) = dataTable.COM(i,2);                                  % new indexed COM y is non zero COM
        counter = counter + 1;                                                        % Increase counter
        
    end
    
end
%%


for n = 2:max(dataTable.Cell)
    contained = 0;
    location = 0;
    [contained, location] = ismember(n, dataTable.Cell);
    
    if contained == 1    
        % Pick some boundaries
        COM = COM_GFP_new(location,:);                                                    % The center of mass of the cell is one index behind the index of the cell, due to the unique 0 cell index, that is in fact, not a cell.
        top_bound = COM(1,2) - 200;                                                         % top bound -100 pixels of COM                    
        bottom_bound = COM(1,2) + 200;                                                      % bottom bound +100
        left_bound = COM(1,1) - 200;
        right_bound = COM(1,1) + 200;

        %% Make sure boundaries are valid indices
        if left_bound <= 0
            left_bound = 1;
        end

        if right_bound > 1344
            right_bound = 1344;
        end

        if top_bound <= 0
            top_bound = 1;
        end

        if bottom_bound > 1024
            bottom_bound = 1024;
        end

        % Start testing points in the cell
         rand_point_x = 1;
         rand_point_y = 1;
         counter = 1;

         while counter <= 2000
             while waterL_GFP(rand_point_y, rand_point_x) ~= n
                 rand_point_x = randi([left_bound, right_bound],1);
                 rand_point_y = randi([top_bound,bottom_bound],1);
                 if rand_point_x == COM(1,1) && rand_point_y == COM(1,2)
                     rand_point_x = randi([left_bound, right_bound],1);
                     rand_point_y = randi([top_bound, bottom_bound],1);
                 end              
             end

             if waterL_GFP(rand_point_y, rand_point_x) == n
                 
                  x1 = rand_point_x;                                                        % x position of random point inside cell.
                  x2 = COM(1,1);                                                            % x position of COM of cell.
                  y1 = rand_point_y;                                                        % y position of random point inside cell.
                  y2 = COM(1,2);                                                            % y position of COM of cell.
                  points(counter,:) = [x1 y1];                                              % Collect and save the position of the random point inside cell.
                  comPtDist(counter) = sqrt((x1 - x2)^2 + (y1 - y2)^2);                     % Euclidean distance between COM and random point inside cell.

                  % Figure out the vector 
                  if ~(x1 == x2 && y1 == y2)
                      vector(counter,:) = [(x1 - x2) (y1 - y2)];                                % Random point inside cell - COM of cell
                      phase(counter) = atan2(y1-y2,x1-x2)*180/pi;
                      % Normalize vector
                      vector_mag(counter) = sqrt(vector(counter,1)^2 + vector(counter,2)^2);     
                      vector_norm(counter,:) = vector(counter,:)./vector_mag(counter);           % Normalize the vector by dividing it by the vector magnitude. 
                      
                      % Use the vector to find nearest edge
                      new_point_x = rand_point_x;                                               % declare x coordinate of dummy point so that we can get into a loop to start moving forward one unit x component, searching for an edge.
                      new_point_y = rand_point_y;                                               % declary y coordinate of dummy point so that we can get into a loop to start moving forward one unit y component, search for an edge.

                      while waterL_GFP(int16(new_point_y), int16(new_point_x)) == n             % While the new point is still contained in the cell

                          new_point_x = new_point_x + vector_norm(counter,1);                   % Move forward one unit x component
                          new_point_y = new_point_y + vector_norm(counter,2);                   % move forward one unit y component

                          % If new point is now at an edge take intensity
                          % move to next point
                          if waterL_GFP(int16(new_point_y), int16(new_point_x)) ~= n            % If new point is no longer in the cell

                              % calculate the distance from COM to edge point
                              x_edge = double(int16(new_point_x));                              % make the new x edge point an integer, and then make it a double so that matlab can manipulate it.
                              y_edge = double(int16(new_point_y));                              % make the new y edge point an integer, and then make it a double so that matlab can manipulate it.
                              edge_points(counter,:) = [x_edge y_edge];                         % make a new data structure for the edge points and keep it.          
                              edge_distance(counter) = sqrt((x_edge - x2)^2 + (y_edge - y2)^2); % Distance from the edge point and the COM of the cell in pixels.

                              % Calculate a normalized radius from center of mass
                              normal_r(counter) = comPtDist(counter)/edge_distance(counter);    % The magnitude of the vector between the COM and point where data is being collected, normalized by the distance between COM and nearest edge along that vector.
                          end
                      end
                  end
                  
                   % If x1 and y1 = x2 and y2 r* = 0
                  if x1 == x2 && y1 == y2
                      normal_r(counter) = 0;
                  end

                  % Take an intensity measurement
                  intensity(counter) = Stain_image(rand_point_y, rand_point_x);             % Take a reading of "intensity" (image value) from the random point inside the cell and save it.

                  % Increase counter
                  counter = counter + 1;


                  % Reset random point
                  rand_point_x = 1;
                  rand_point_y = 1;
             end
         end

         median_intensity = median(median(Stain_image));                                                 % mean intensity of collection of points from inside the cell
         norm_intensity = (intensity-median_intensity);                                       % each individual intensity normalized by the mean intensity
         sub_intensity = intensity - median_intensity;                                       % each individual intensity - normalized intensity


         %%
         save(['all_cell_',num2str(n),'_35kPa_protein.mat'],'normal_r','norm_intensity','comPtDist','phase');
    end
end