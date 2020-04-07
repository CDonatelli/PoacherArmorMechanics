function sOut = FishBurrowing(videoName)

sOut = struct;
%vars = who;
%eval(['sOut=',cell2mat(vars(1))])

FileNamePrefix = VideoReader(videoName);
FrNum = FileNamePrefix.NumberOfFrames;

Lines(1).Frame = [];
Lines(1).MidLine = [];
X=[];
Y=[];

    ImStart = read(FileNamePrefix,1);
    
    % Get a rectangle to limit search area to area the fish is swimming in
    rect = CropVideo(ImStart);
    sOut.rect = rect;
    
    % Get the rough levels of the background and of the fish  
    % (assuming background is plain and uniform for now)
    % If you have a BW video instead of RGB swich the "GetLevels" version
    [BackLev, FishLev] = GetLevels(imcrop(ImStart,rect));
    % or
    %[BackLev, FishLev] = GetLevelsBW(imcrop(ImStart,rect));
    
    sOut.backLevel = BackLev; sOut.fishLevel = FishLev;
    
    %This bit assumes the fish is dark and the background is light
    ThreshLevel = median([BackLev(2),FishLev(1)])/255;

    % Use this if the fish is light and the back is dark
    %median(BackLev(1),FishLev(2))
    
    % Test if the threshhold is good and edit if it isnt
        choice = 1;
        RawImage = read(FileNamePrefix,1);%get the first image to allow user to click the fish    
        RawImage = imcrop(RawImage, rect);
        BinaryImage = ProcessImage(RawImage,ThreshLevel);
        imshow(BinaryImage)
        while choice ~= 2
            choice = input('Is this good? (1 = brighten fish, 0 = dim  fish, 2 = good): ');
            if choice == 0
                ThreshLevel = ThreshLevel - 0.025;
                BinaryImage = ProcessImage(RawImage,ThreshLevel);
                imshow(BinaryImage)
            elseif choice == 1
                ThreshLevel = ThreshLevel + 0.025;
                BinaryImage = ProcessImage(RawImage,ThreshLevel);
                imshow(BinaryImage)
            else
                ThreshLevel = ThreshLevel;
            end
        end
        
    answer = questdlg('Which tracking mode do you want to use?', ...
    'Tracking Modes', ...
    'FastTrack','SensiTrack','SensiTrack');
    switch answer
        case 'FastTrack'
            TrackOpt = 2;
        case 'SensiTrack'
            TrackOpt = 1;
    end

for Index = 1:FrNum
    RawImage = read(FileNamePrefix,Index);%get the first image to allow user to click the fish    
    RawImage = imcrop(RawImage, rect);
    BinaryImage = ProcessImage(RawImage,ThreshLevel);       %blur the image, threshold and invert
    LabelImage = bwlabeln(BinaryImage,4);       %label the image to use image props          
    imshow(BinaryImage);
    ImageStats = regionprops(LabelImage,'all'); %get stats on the labelled image

    %if this is the first frame then get teh nose, otherwise use teh front
    %point from the last image as teh temporary nose.
    if (size(X,1) == 0 || BinaryImage(round(Y),round(X)) == 0)
        [X Y] = ginput(1);  %get the location of the fish
    end
    
    FishRegion = LabelImage(round(Y),round(X)); %get the region number of the fish
    FishImage = BinaryImage;%.*(LabelImage==FishRegion);  %kill all the rest of the binary image
    imshow(FishImage)       %show just the fish to make sure all is well
    hold on;
    plot(X,Y,'or'); %show the dot the user clicked
    
    %figure out which way the rest of the fish lies.
    %going to assume it is in the direction of the centroid from the point
    %on the head.  Setting that general direction establishes a polarity
    %for the midline search to proceed down the animal rather than from the
    %head to the nose.
    XTemp=ImageStats(FishRegion).Centroid(1)-X;
    YTemp=ImageStats(FishRegion).Centroid(2)-Y;
    [AngleToNext,D] = cart2pol(XTemp, YTemp);
    
    %use teh general direction of teh rest of the body to find the 'nose'.
    %for ease this will be the point furthest from the user clicked point
    %in the opposite direction from the centroid.
    
    Nose = FindNose(FishImage, X, Y, AngleToNext+pi);
    X=Nose(1);
    Y=Nose(2);
    plot(X,Y,'og');
    
    % set the radius for the midline finding circle 
    Radius = RadFind(BinaryImage,X,Y);   
    
    %find a center for the drawn circle that is 2*Radius in the opposite
    %direction from the centroid
    [TempCenter(1), TempCenter(2)] = pol2cart(AngleToNext+pi,2*Radius);
    TempCenter(1) = TempCenter(1)+X;
    TempCenter(2) = TempCenter(2)+Y;
    plot(TempCenter(1),TempCenter(2),'og');
    
    %this finds a circle on the clicked point and plots it for debug
    %coordinates of a circle centered on the user point
    FullCircle = GetArc(TempCenter(1),TempCenter(2),3*Radius,0,2*pi);     
    plot(FullCircle(:,1),FullCircle(:,2),'.b');     %debug code
    hold on
    %180 degrees of arc centered on the user point
    FullArc = GetArc(TempCenter(1),TempCenter(2),3*Radius,...
        AngleToNext-pi/2,AngleToNext+pi/2); 
    %shows arc that crosses fish body posterior to current point
    plot(FullArc(:,1),FullArc(:,2),'.r');    
    
    FullArc = TrimArc(FishImage,FullArc);    % removes out of bounds values
    FishArc = FindWhite(FishImage,FullArc);  % narrow the list points to those on the fish

    %first point of midline is user point
        Centers=[X, Y];     
    %get the rest of the midline
        Centers = FindMidPoint(FishImage,FishArc,Centers);    
    %the auto tracking works better if the next starting point is a bit
    %back from the nose.
        X = Centers(TrackOpt,1); Y = Centers(TrackOpt,2);
    
    %%%% THIS IS MEANT TO FIX IT IF THE MIDLINE BOUNCES
%     if any(diff(Centers(:,1)) <= 0)
%         screwed = false; over = [];
%         for i = 2:length(Centers)
%             if Centers(i-1,1) > Centers(i,1)
%                 screwed = true;
%                 over = i-1;
%                 break
%             else
%                 screwed = false;
%             end
%         end
%         Centers(over:end,:) = [];
%     end
    %%%%%%%%%%%%%%%%%%%%%%
    
    
    Lines(Index).Frame=Index;       %save data in the output structure
    Lines(Index).MidLine=Centers;
    hold off    %allow the image to be redrawn
end

close   %close the image 
hold on %see multiple traces
for Index = 1:size(Lines,2)
    plot(Lines(Index).MidLine(:,1),Lines(Index).MidLine(:,2),'-k');
end

sOut.midLines = Lines;
eval([videoName(1:end-4), '= sOut'])
save(videoName(1:end-4), videoName(1:end-4));

%this finds a radius of a circle centered on a point that overlaps both
%sides of the fish. written by Cassandra Donatelli 2014

function R = RadFind(Image,X,Y)
        R = 5;
        changes = 0;
        [m,n] = size(Image);
        while changes < 2
            circ = GetArc(X,Y,R,0,2*pi);
            del = [];
            for i = 1:length(circ)
                if circ(i,1) > n || circ(i,1) <=0 || circ(i,2) > m ...
                        || circ(i,2) <=0
                    del = [del,i];
                end
            end
            circ(del,:) = []; 
            for i = 1:length(circ)
                vals(i) = Image(circ(i,2), circ(i,1));
            end
            for i = 1:length(vals)-1
                if vals(i) ~= vals(i+1)
                    changes = changes +1;
                else changes = changes;
                end
            end
            %plot(circ(:,1), circ(:,2)); hold on
            R = R+1;
        end
        
%returns the list of points in an arc centered at X,Y of Radius R. A 2014
%revision to this code removes the duplications that arise from the
%floor step.

function Arc = GetArc(X,Y,R,PhiStart,PhiFinish)
Arc=[];
for Theta = PhiStart:2*pi/720:PhiFinish    %make a reading each degree                                            
    [DX, DY] = pol2cart(Theta,R);   %get the cartesian coordinates of the polar expression
    Arc = [Arc;DX DY];      %save the coordinates in the list
end

Arc(:,1) = (Arc(:,1) + X);  %add the center X value
Arc(:,2) = (Arc(:,2) + Y);  %add the center Y value
Arc=floor(Arc); %round everything down
NewArc = [Arc(1,:)];
for i = 2:size(Arc,1)
    if Arc(i,1) ~= Arc(i-1,1) || Arc(i,2) ~= Arc(i-1,2)
        NewArc=[NewArc; Arc(i,1) Arc(i,2)];
    end
end
Arc = NewArc;

%returns list of points where the binary image equals 1
function White = FindWhite(Frame,Points)
[m,n] = size(Frame);
del = [];
for i = 1:length(Points)
    if Points(i,1) > n || Points(i,1) <= 0 || Points(i,2) > m ...
            || Points(i,2) <= 0
        del = [del,i];
    end
end
Points(del,:) = [];
White =[];
for Index = 1:size(Points,1)
    X=Points(Index,1);  %get the x and y values of the point
    Y=Points(Index,2);
    if Frame(Y,X)==1    %check if the x,y coordinate is white 
        White=[White;X Y];  %add the coordinate to the list
%         plot(X,Y,'ro');
    end
end
    
%Takes an arc through the fish, calculates the midpoint then sets a radius
%for a circle and cuts another arc. This works recursively getting a new
%arc over and over again until it runs out of fish. 
function Centers = FindMidPoint(FishImage,FishArc,Centers)
i = 1;
while (size(FishArc,1) > 1) && (i < 500)
    MidPoint = FishArc(floor(size(FishArc,1)/2),:);  %the midpoint
    Centers = [Centers;MidPoint];   %add the midpoint to the list

    NewRadius = RadFind(FishImage, MidPoint(1), MidPoint(2));

    %get the arc specified by the new midpoint and radius. The hemicircle
    %should start at right angles to a line between the two previous
    %centers
    CenterAngle = cart2pol(Centers(end,1)-Centers(end-1,1),Centers(end,2)-Centers(end-1,2));
    imshow(FishImage)
    %As in th initail step take a point along the line between the two
    %recent centers, 2* radius away from the tail
    %first find the temporary circle center
    [TempCenter(1), TempCenter(2)] = pol2cart(CenterAngle+pi,2*NewRadius);
    TempCenter(1) = TempCenter(1)+MidPoint(1);
    TempCenter(2) = TempCenter(2)+MidPoint(2);
    plot(TempCenter(1),TempCenter(2),'og');
    hold on
    RawArc = GetArc(TempCenter(1),TempCenter(2),NewRadius*3,CenterAngle-pi/2, CenterAngle+pi/2);
    RawArc = TrimArc(FishImage,RawArc);
    FishArc = FindWhite(FishImage,RawArc);   %find the intersection of the arc and the image
    plot(MidPoint(:,1),MidPoint(:,2),'.g'); drawnow; %plot a green circle on the midline
    plot(RawArc(:,1),RawArc(:,2),'r'); drawnow;%debug to show the arc
    %Centers = FindMidPoint(FishImage,FishArc,Centers);   %recursively look for midpoints
    i = i+1;
end

% blur and crop the image then invert and binary it
function FrameOut = ProcessImage(Frame,Level)
h = ones(5,5) / 25;     %blur the image to kill line artifacts
BlurredImage = imfilter(Frame,h);
%CroppedImage = double(BlurredImage);%(65:405,80:660); %this removes the borders from the images
                                             % assumes motionscope frame
%Level = graythresh(BlurredImage)*.8;         %set threshold a little darker than the auto computed one
FrameOut = ~im2bw(BlurredImage,Level);       %make image binary and invert it so fish is white
FrameOut = bwareaopen(FrameOut, 50);
FrameOut(1:3,:) = [];
FrameOut(end-2:end,:) = [];
FrameOut(:,1:3) = [];
FrameOut(:,end-2:end) = [];

%Smooth broken bits of fish
    se = strel('disk',10);
    FrameOut = imclose(FrameOut,se);


  
%this makes sure that none of the values in an arc point to locations that
%can't exist in the image.
function Arc = TrimArc(Image, Arc)
    Arc(Arc(:,1)>size(Image,2)-1,1) = size(Image,2)-1;
    Arc(Arc(:,2)>size(Image,1)-1,2) = size(Image,1)-1;
    Arc(Arc(:,1)<0,1) = 0;
    Arc(Arc(:,2)<0,1) = 0;
       
%find the white point furthest from the point clicked on by the user in 
%the direction away from the centroid    
function Nose = FindNose(Frame, X, Y, Angle)
    [XTemp,YTemp] = pol2cart(Angle,1);
    if X == 0
        X = X+1;
    end
    if Y == 0
        Y = Y+1;
    end
    while ~Frame(round(Y),round(X))==0
        X = X+(1.01*XTemp);
        Y = Y+(1.01*YTemp);
    end
    Nose = [X,Y];    
    
function rect = CropVideo(im)
    disp('Select the portion of the frame the fish swims through');
    choice = 0;
    while choice == 0
        imshow(im)
        rect = getrect;
        im2 = imcrop(im,rect);
        imshow(im2)
        choice = input('Does this look right? :');
    end
    
function [Back, Obj] = GetLevels(im)
    
    % Read in original RGB image.
    rgbImage = im;
    % Extract color channels.
    redChannel = rgbImage(:,:,1); % Red channel
    greenChannel = rgbImage(:,:,2); % Green channel
    blueChannel = rgbImage(:,:,3); % Blue channel
    % Create an all black channel.
    allBlack = zeros(size(rgbImage, 1), size(rgbImage, 2), 'uint8');
    % Create color versions of the individual color channels.
    just_red = cat(3, redChannel, allBlack, allBlack);
    just_green = cat(3, allBlack, greenChannel, allBlack);
    just_blue = cat(3, allBlack, allBlack, blueChannel);
    % Recombine the individual color channels to create the original RGB image again.
    recombinedRGBImage = cat(3, redChannel, greenChannel, blueChannel);
    % Display them all.
    subplot(3, 3, 2);
    imshow(rgbImage);
    fontSize = 20; title('Original RGB Image', 'FontSize', fontSize)
    subplot(3, 3, 4);
    imshow(just_red); title('Red Channel in Red', 'FontSize', fontSize)
    subplot(3, 3, 5);
    imshow(just_green); title('Green Channel in Green', 'FontSize', fontSize)
    subplot(3, 3, 6);
    imshow(just_blue); title('Blue Channel in Blue', 'FontSize', fontSize)
    subplot(3, 3, 8);
    imshow(recombinedRGBImage); title('Recombined to Form Original RGB Image Again', 'FontSize', fontSize)
    
    answer = questdlg('Which channel had the most contrast?', ...
        'Color Channels', ...
        'Red','Green','Blue','Blue');
    switch answer
        case 'Red'
            channel = 1;
        case 'Green'
            channel = 2;
        case 'Blue'
            channel = 3;
    end

    close all
    
    OBlu = []; BBlu = [];
    imshow(im); hold on
    disp('Get Fish Levels');
    [Xo Yo] = getpts;
    plot(Xo,Yo,'bo');
    hold on
    for i = 1:length(Xo)
        O = impixel(im,Xo(i),Yo(i));
        OBlu = [OBlu,O(channel)];
    end
    disp('Get background Levels');
    [Xb Yb] = getpts;
    plot(Xb,Yb,'ro');
    hold on
    for i = 1:length(Xb)
        B = impixel(im,Xb(i),Yb(i));
        BBlu = [BBlu,B(channel)];
    end
    
    
    
    MaxObj = max(OBlu); MinObj = min(OBlu);
    MaxBac = max(BBlu); MinBac = min(BBlu);
    
    % If the levels are overlapping, find the average
    % For now assuming that the background and fish are pretty different
    % so the only overlapping levels considered are as follows:
    % MaxFish > MaxBackground > MinFish > MinBackground
    % MaxBackground > MaxFish > MinBackground > MinFish
    % Looking to create an order that looks like one of the following:
    % MaxFish > MinFish > MaxBackground > MinBackground
    % MaxBackground > MinBackground > MaxFish > MinFish
    
    if MaxObj >= MinBac 
        if MaxBac >= MinObj 
            Avg = round(mean([MaxBac MinObj]));
            MinObj = Avg;
            MaxBac = Avg-1;
        end
    end
    if MaxBac >= MinObj 
        if MaxObj >= MinBac   
            Avg = round(mean([MaxObj MinBac]));
            MaxObj = Avg;
            MinBac = Avg+1;
        end
    end
    
    %Now I need to account for the fact that the user might not select the
    %full range of points on the fish (I did this and it leads to incorrect
    %D values which means incorrect wobble measurements. This is a cheat-y
    %fix but it works with my videos. This part of the code is less
    %adaptable for other users (especially if they are filming in a
    %background that is not super different from the fish). 
    %Assumes one of these two cases:
    % MinBackground < MaxBackground < MinFish < MaxFish 
    % MinFish < MaxFish < MinBackground < MaxBackground
    
    if MinObj > MaxBac
        MinObj = MinObj - ((MinObj-MaxBac)/2);
        MaxObj = MaxObj + ((MinObj-MaxBac)/2);
    end
    if MinBac > MaxObj
        MaxObj = MaxObj + ((MinBac-MaxObj)/2);
        MinObj = MinObj - ((MinBac-MaxObj)/2);
    end
    
    hold off
    Back = [MaxBac, MinBac]; Obj = [MaxObj, MinObj];
    
function [Back, Obj] = GetLevelsBW(im)   
    % Read in original RGB image.
    rgbImage = im;
        
    OBlu = []; BBlu = [];
    imshow(im); hold on
    disp('Get Fish Levels');
    [Xo Yo] = getpts;
    plot(Xo,Yo,'bo');
    hold on
    for i = 1:length(Xo)
        O = impixel(im,Xo(i),Yo(i));
        OBlu = [OBlu,O(1)];
    end
    disp('Get background Levels');
    [Xb Yb] = getpts;
    plot(Xb,Yb,'ro');
    hold on
    for i = 1:length(Xb)
        B = impixel(im,Xb(i),Yb(i));
        BBlu = [BBlu,B(1)];
    end
    
    MaxObj = max(OBlu); MinObj = min(OBlu);
    MaxBac = max(BBlu); MinBac = min(BBlu);
    
    % If the levels are overlapping, find the average
    % For now assuming that the background and fish are pretty different
    % so the only overlapping levels considered are as follows:
    % MaxFish > MaxBackground > MinFish > MinBackground
    % MaxBackground > MaxFish > MinBackground > MinFish
    % Looking to create an order that looks like one of the following:
    % MaxFish > MinFish > MaxBackground > MinBackground
    % MaxBackground > MinBackground > MaxFish > MinFish
    
    if MaxObj >= MinBac 
        if MaxBac >= MinObj 
            Avg = round(mean([MaxBac MinObj]));
            MinObj = Avg;
            MaxBac = Avg-1;
        end
    end
    if MaxBac >= MinObj 
        if MaxObj >= MinBac   
            Avg = round(mean([MaxObj MinBac]));
            MaxObj = Avg;
            MinBac = Avg+1;
        end
    end
    
    %Now I need to account for the fact that the user might not select the
    %full range of points on the fish (I did this and it leads to incorrect
    %D values which means incorrect wobble measurements. This is a cheat-y
    %fix but it works with my videos. This part of the code is less
    %adaptable for other users (especially if they are filming in a
    %background that is not super different from the fish). 
    %Assumes one of these two cases:
    % MinBackground < MaxBackground < MinFish < MaxFish 
    % MinFish < MaxFish < MinBackground < MaxBackground
    
    if MinObj > MaxBac
        MinObj = MinObj - ((MinObj-MaxBac)/2);
        MaxObj = MaxObj + ((MinObj-MaxBac)/2);
    end
    if MinBac > MaxObj
        MaxObj = MaxObj + ((MinBac-MaxObj)/2);
        MinObj = MinObj - ((MinBac-MaxObj)/2);
    end
    
    hold off
    Back = [MaxBac, MinBac]; Obj = [MaxObj, MinObj];
        