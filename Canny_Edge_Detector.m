% Canny Edge Detector
% Author: Rohit Ainapure

x1 = imread('/home/csgrad/rohitain/Desktop/Lena.jpg');
x2 = imread('/home/csgrad/rohitain/Desktop/House.jpg');

x_db = im2double(x1);

% sqr function to square the values of pixels for Gradient calculations
sqr = @(x) x.^2;

[nr nc] = size(x_db);

h_test = [2 4  5  4  2;
          4 9  12 9  4;
          5 12 15 12 5;
          4 9  12 9  4;
          2 4  5  4  2];
      
h_test = h_test.*(1/159);

h = fspecial('gaussian',5,2);
x_gaussian = conv2(x_db,h,'same');

% figure,
% imshow(x_gaussian);
x_db = padarray(x_db,[1,1],0,'pre');
first_derivative_x(1:1:nr,1:1:nc) = (sqr(((x_db(2:1:nr+1,2:1:nc+1) - (x_db(1:1:nr,2:1:nc+1))))) + (sqr((x_db(2:1:nr+1,2:1:nc+1)) - (x_db(2:1:nr+1,1:1:nc)))));
G1 = sqrt(first_derivative_x);

% Vertical Gradient
sobel_x = [-1 0 1;
           -2 0 2;% thresh_high = 0.1;
% thresh_low = 0.01;

           -1 0 1];
       
% Horizontal Gradient
sobel_y = [-1 -2 -1;% thresh_high = 0.1;
% thresh_low = 0.01;

           0  0  0;
           1  2  1];
       
G_X = conv2(x_gaussian,sobel_x,'same');
% figure,
% imshow(G_X);

G_Y = conv2(x_gaussian,sobel_y,'same');
% figure,
% imshow(G_Y);


% Computes the Gradient
G = sqrt((G_X.^2) + (G_Y.^2));
% figure,
% imshow(G);

G2 = G1-G;
Theta = zeros(nr,nc);

% Computes the Edge Direction in degrees
for i=1:nr
    for j=1:nc
        Theta(i,j) = atand(G_Y(i,j) ./ G_X(i,j));
        if Theta(i,j) < 0
            Theta(i,j) = 360 + Theta(i,j);
        end
    end
end% ideal = edge(x_db,'canny');
% figure,
% imshow(ideal);


    
% figure,
% imshow(Theta);

% Define ranges for Theta for approximation of edge direction as follows:
% if 0 <= Theta < 22.5 OR 157.5 <= Theta <= 180 then 
% Theta = 0

% if 22.5 <= Theta < 67.5 then 
% Theta = 45

% if 67.5 <= Theta < 112.5 then
% Theta = 90

% if 112.5 <= Theta < 157.5 then% thresh_high = 0.1;
% thresh_low = 0.01;

% Theta = 135

[nrt nct] = size(Theta);
X_EDGE_DIRECTION = zeros(nr,nc);

for i=1:nrt
    for j=1:nct
        if Theta(i,j) >337.5 || Theta(i,j) >=0 && Theta(i,j) < 22.5 || Theta(i,j) > 157.5 && Theta(i,j) <= 202.5
            X_EDGE_DIRECTION(i,j) = 0;
        end
        
        if (Theta(i,j) >= 22.5 && Theta(i,j) < 67.5) || (Theta(i,j) > 202.5 && Theta(i,j) <= 247.5)
            X_EDGE_DIRECTION(i,j) = 45;
        end
        
        if (Theta(i,j) >= 67.5 && Theta(i,j) < 112.5) || (Theta(i,j) > 247.5 && Theta(i,j) <= 292.5) 
            X_EDGE_DIRECTION(i,j) = 90;
        end 
        
        if (Theta(i,j) >= 112.5 && Theta(i,j) <= 157.5) || (Theta(i,j) > 292.5 && Theta(i,j) <= 337.5)
            X_EDGE_DIRECTION(i,j) = 135;
        end
    end
end

% figure,
% imshow(X_EDGE_DIRECTION);

[nrtx nctx] = size(X_EDGE_DIRECTION);
X_Mark = G;

G = padarray(G,[1,1],0,'both');

% Non Maximum Suppression, removal of weak edges
% Iterate through the X_EDGE_DIRECTION matrix. For each pixel,
% consider the direction of orientation. Compare the pixels values
% along that direction of orientation. For Ex. is orientation is 45deg
% then check the pixel vaules at North east & South west.
% If the current pixel value is largest, then preserve it else remove it.

for i=2:nr+1
    for j=2:nc+1
        
        if X_EDGE_DIRECTION(i-1,j-1) == 0
            if (G(i,j) < G(i,j+1) || (G(i,j) < G(i,j-1)))
                X_Mark(i-1,j-1) = 0;
            end
        end
        
        if X_EDGE_DIRECTION(i-1,j-1) == 45
            if (G(i,j) < G(i-1,j+1)) || (G(i,j) < G(i+1,j-1))
                X_Mark(i-1,j-1) = 0;
            end
        end
        
        if X_EDGE_DIRECTION(i-1,j-1) == 90
            if ((G(i,j) < G(i-1,j)) || (G(i,j) < G(i+1,j)))
                X_Mark(i-1,j-1) = 0;
            end
        end
        
         if X_EDGE_DIRECTION(i-1,j-1) == 135
            if (G(i,j) < G(i-1,j-1)) || (G(i,j) < G(i+1,j+1))
                X_Mark(i-1,j-1) = 0;
            end
         end
    end
end


X_NMS = X_Mark;       

% Result after Non Maximum Suppression
% figure,
% imshow(X_Mark);

% Thresholding
% Set two thresholds - Hysterisis

% For Local Thresholding
thresh_high = 0.1655;           
thresh_low = 0.13;              

X_Hyst = X_Mark;

% Hysterisis
X_Hyst_Pad = padarray(X_Hyst,[1 1],0,'both');

for i=2:nr+1
    for j=2:nc+1
        if X_Hyst_Pad(i,j) >= thresh_high
            X_Hyst(i-1,j-1) = 1;
        end
        if X_Hyst_Pad(i,j) < thresh_high && X_Hyst_Pad(i,j) >= thresh_low
            if X_Hyst_Pad(i,j+1) >= thresh_high || X_Hyst_Pad(i-1,j) >=thresh_high || X_Hyst_Pad(i,j-1) >= thresh_high || X_Hyst_Pad(i+1,j) >= thresh_high || X_Hyst_Pad(i+1,j+1) >=thresh_high || X_Hyst_Pad(i-1,j+1) >= thresh_high || X_Hyst_Pad(i-1,j-1) >= thresh_high || X_Hyst_Pad(i+1,j-1) >= thresh_high
                X_Hyst(i-1,j-1) = 1;
            else
                X_Hyst(i-1,j-1) = 0;
            end  
        end
        if X_Hyst_Pad(i,j) < thresh_low
            X_Hyst(i-1,j-1) = 0;
        end
    end
end

X_Final = X_Hyst;


%Ideal result
ideal = edge(x_db,'canny');
figure,
imshow(ideal);


%X_Diff = ideal - X_Final;

% figure,
% imshow(x_db);

% Result after Hysterisis
figure,
imshow(X_Final);