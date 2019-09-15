%% VIDEO READING (DO ONCE)
read = input('Do you need to read in the video file? (Type Y for yes, anything else for No): ','s');
videoname = input('Enter the filename, with extension, for the video to process: ','s');%'S1010005.MP4'; % MAKE SURE THE FILE IS IN YOUR WD
fileprefix = [videoname,'_'];
if read=='Y'
    video = VideoReader(videoname);% REMEMBER TO MODIFY THE FILENAME
    ii = 1;
    while hasFrame(video)
        img = readFrame(video);
        filename = [fileprefix,sprintf('%03d',ii), '.jpg'];% REMEMBER TO MODIFY THE FILENAME
        imwrite(img,filename);
        ii = ii+1;
    end
end

%% INITIALIZING AND CUTTING
%imageNames = dir(fullfile('images','*.jpg'));
imageNames = dir([fileprefix,'*.jpg']); % REMEMBER TO MODIFY THE FILENAME
imageNames = {imageNames.name};
maxframe = length(imageNames);
fontSize = 11;

xname = cell(1,maxframe); % all frames
yname = cell(1,maxframe); % all frames
brightThresh = 253;

%%
disp('Now the frames will be transformed to grayscale, and the LED markers detected in each frame.')
input('Press ENTER to continue: ')
for ii=1:maxframe
    
    readframe = [fileprefix,sprintf('%03d',ii),'.jpg'];
    [I,colormap] = imread(readframe);

    grayI = rgb2gray(I); % Convert to grayscale
    [rowsGray, columnsGray] = find(grayI == 255);

    binaryIgray = grayI(:,:) >= brightThresh;
    brightObjectsMask = binaryIgray;
    imshow(brightObjectsMask)
    [y,x] = find(brightObjectsMask);
    %x = transpose(x);
    %y = transpose(y);
    xname{1,ii} = x;
    yname{1,ii} = y;
    progress = uint16(ii/maxframe*100);
    disp(sprintf('Progress of binary masking: %02d%%',progress))
end


%% Manual Regioning
for ii = 1
    readframe = [fileprefix,sprintf('%03d',ii),'.jpg'];% REMEMBER TO MODIFY THE FILENAME
    [I,colormap] = imread(readframe);

    grayI = rgb2gray(I);
    [rowsGray, columnsGray] = find(grayI == 255);

    binaryIgray = grayI(:,:) >= brightThresh;
    brightObjectsMask = binaryIgray;
    
    input('Now you will begin demarcating regions containing each marker using the masked first frame.\nPress ENTER to continue: ');
    numberMarkers = input('How many markers are present in the image?: ');
    disp('You will use the data cursor tool to manually examine data values, and select the initial ')
    disp('locations for each marker region. ')
    disp('From the following image, select approximate marker locations for each marker, one by one. ')
    disp('Use alt-click to select each marker location, then right click and ')
    disp('save the location using Export to Workspace.')
    disp('Save as a variable "markers". ') 
    input('Press ENTER to begin: ');
    f = figure;
    imshow(brightObjectsMask, []);
    caption = sprintf('First Frame: Binary mask of bright regions');
    title(caption, 'FontSize', fontSize);
    datacursormode on
    waitfor(f);
    % If done correctly, markers(i).Position is the position (x,y) of
    % marker number i in the first frame. 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NB: Marker 1 will be the last marker clicked, so the first marker
    % clicked will be Marker # = numberMarkers.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



%% Point tracking using Regions (rather than splitting)
disp('Point tracking will now commence using binary masks and initial marker locations.')
input('Press ENTER to continue: ')
xs = cell(numberMarkers,maxframe);
ys = cell(numberMarkers,maxframe);
centers = cell(numberMarkers,maxframe);
counter = 0;
searchRadius = 100;
for ii = 1:maxframe
    progress = uint16(ii/maxframe*100);
    disp(sprintf('Progress of point tracking: %02d%%',progress))
    for jj = 1:length(xname{1,ii})
        
        if ii==1
            for kk = 1:numberMarkers
                if sqrt((xname{1,ii}(jj)-markers(kk).Position(1))^2+(yname{1,ii}(jj)-markers(kk).Position(2))^2)<=searchRadius
                    
                    xs{kk,ii} = [xs{kk,ii} xname{1,ii}(jj)];
                    ys{kk,ii} = [ys{kk,ii} yname{1,ii}(jj)];

                end 
            end 
        end
        if ii > 1
            for kk = 1:numberMarkers
                mm = 35;
                if sqrt((xname{1,ii}(jj)-mean(xs{kk,ii-1}))^2+(yname{1,ii}(jj)-mean(ys{kk,ii-1}))^2)<=searchRadius

                    xs{kk,ii} = [xs{kk,ii} xname{1,ii}(jj)];
                    ys{kk,ii} = [ys{kk,ii} yname{1,ii}(jj)];
                        
                end
            end
        end
    end
end
% Centers
centersxMatrix = zeros(numberMarkers,maxframe);
centersyMatrix = zeros(numberMarkers,maxframe);
for ii = 1:maxframe
    for kk = 1:numberMarkers
        centers{kk,ii} = [mean(xs{kk,ii}) mean(ys{kk,ii})];
        centersxMatrix(kk,ii) = mean(xs{kk,ii}); 
        centersyMatrix(kk,ii) = mean(ys{kk,ii});
    end
end
%% Spatial Calibrating 


% 10 'cm' should correspond to a true distance measured,
% corresponding to the distance between close corners of red markers
disp('Now spatial calibration will be applied to the frames to convert pixel distances ')
disp('to physical distances. ')
disp('You will be prompted for the two markers between which you measured the physical distance.')
input('Press ENTER to continue: ')

% The two white LEDs were chosen for the shake table test videos provided.
markerA = input('Enter the number of the first marker (using the order you selected the markers): ');
markerB = input('Enter the number of the second marker (using the order you selected the markers): ');

% The measured distance was ~20 cm for the shake table test videos provided.
units = input('Enter the physical units used: ','s'); 
distanceInUnits = input('Enter the measured physical distance between the two selected markers: '); 

distanceInPixels = sqrt((centers{markerB,1}(1)-centers{markerA,1}(1)).^2 + (centers{markerB,1}(2)-centers{markerA,1}(2)).^2);
distancePerPixel = distanceInUnits / distanceInPixels;

%xInUnits = cell(numberMarkers,maxframe);
dxInUnitsMatrix = zeros(numberMarkers,maxframe);
%yInUnits = cell(numberMarkers,maxframe);
dyInUnitsMatrix = zeros(numberMarkers,maxframe);
%deltaInUnits = cell(numberMarkers,maxframe);
%deltaxInUnits = cell(numberMarkers,maxframe);
%deltayInUnits = cell(numberMarkers,maxframe);

for i=2:maxframe
    for k=1:numberMarkers
        %deltaInUnits{k,i} = distancePerPixel*sqrt( (centers{k,i}(1)-centers{k,i-1}(1)).^2 + (centers{k,i}(2)-centers{k,i-1}(2)).^2);
        %deltaxInUnits{k,i} = distancePerPixel*(centers{k,i}(1)-centers{k,i-1}(1));
        %deltayInUnits{k,i} = distancePerPixel*(centers{k,i}(2)-centers{k,i-1}(2));
        
        %xInUnits{k,i} = distancePerPixel*(centers{k,i}(1));
        %yInUnits{k,i} = distancePerPixel*(centers{k,i}(2));
        
        % Below sets the first frame to 0 displacement and subsequent frames 
        % w.r.t that zero point.
        dxInUnitsMatrix(k,i) = distancePerPixel*(centers{k,2}(1))-distancePerPixel*(centers{k,i}(1));
        dyInUnitsMatrix(k,i) = distancePerPixel*(centers{k,2}(2))-distancePerPixel*(centers{k,i}(2));
    end
end

%% PLOTTING TIME SERIES
framerate = 1/60;
vidlength = framerate*maxframe;
t = (0:framerate:vidlength-framerate);

for k = 1:numberMarkers
    xmax = max(abs(dxInUnitsMatrix(k,:)));
    ymax = max(abs(dyInUnitsMatrix(k,:)));
    
    figure()
    plot(t,dxInUnitsMatrix(k,:))
    xlabel('Time (s)')
    ylabel(['x-displacement (',units,')'])
    title(sprintf('Horizontal motion time series for marker %01d',k))
    axis([0 vidlength -xmax xmax])
    
    figure()
    plot(t,dyInUnitsMatrix(k,:))
    xlabel('Time (s)')
    ylabel(['y-displacement (',units,')'])
    title(sprintf('Vertical motion time series for marker %01d',k))
    axis([0 vidlength -ymax ymax])
end
