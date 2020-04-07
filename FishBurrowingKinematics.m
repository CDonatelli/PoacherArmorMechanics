function struct = FishBurrowingKinematics(structfile)
    load(structfile);
    clear 'structfile'
    vars = who;
    eval(['struct=',cell2mat(vars(1))]);
    %This imports the data

    prompt = {'Enter fish length in mm (leave as 0 if unknown)'};
    dlgtitle = 'Fish Length';
    dims = [1 35];
    definput = {'0'};
    length = inputdlg(prompt,dlgtitle,dims,definput);
    length = str2double(length);
    
    mids = struct.midLines;
    nfr = size(mids,2);
    % Calculate the scale of the video using the fish length
    % as the scale bar
        FishPixelLengths = [];
        % Loop through a few frames of the vide and calculate the
        % length of the midline at each frame
        for i = 1:round(nfr/15):nfr
            FishPixelLengths = [FishPixelLengths,...
                arclength(mids(i).MidLine(:,1),mids(i).MidLine(:,2))];
        end
        % Use the known length of the fish and the median of the 
        % midline measurements to calculate the scale of the video
        fishPixels = median(FishPixelLengths);
        
        struct.fishLength = length;
        VidScale = length/fishPixels;
        struct.VidScale = VidScale;
        
    % points to generate data for
    npts = 21;
    % initiating new variables
    x = []; y = []; tailPtCordsY = []; tailPtCordsX = [];
    for i = 1:nfr
        % Generate equation if the midline
        [pts, deriv, funct] = interparc(npts,  mids(i).MidLine(:,1),  ... 
                                        mids(i).MidLine(:,2), 'spline');
        % add those points to an array
        x = [x,pts(:,1)]; y = [y,pts(:,2)];
        
        % usee the above function to find the coordinates of the points
        % of interest in the tail region (to be used later
    end
    struct.X = x; struct.Y = y;
    % figure out time for each frame and make a vector of times
    
    prompt = {'Enter video frame rate.)'};
    dlgtitle = 'Frame Rate';
    dims = [1 35];
    definput = {'0'};
    fr = inputdlg(prompt,dlgtitle,dims,definput);
    fr = str2double(fr);
    

    total = nfr/fr; 
    struct.t = linspace(0,total,nfr)';
    s = linspace(1,length,npts);
    struct.s = s';

%%%% 2D Wave Kinematics
    tailY = smooth(struct.t, struct.Y(20,:));
    
    p = polyfit(struct.t, tailY,7);   % fit line for the tail wave
    yT = polyval(p, struct.t);        % y values for that line
    tailY = tailY - yT;               % subtract those y values from that
                                      % line to get actual amplitude
    
    %%%%% Mini-Peak Finder
        [p1,k1] = findpeaks(tailY,struct.t,'MinPeakProminence',3);
        [p2,k2] = findpeaks(-tailY,struct.t,'MinPeakProminence',3);
        p = [p1; -p2]; k = [k1; k2];
        
        % check if peakfinder did it's damn job
        [k,p] = correctPeakFinder(k,p,struct.t,tailY,['Tail peaks']);
        peaks = [k,p]; peaks = sortrows(peaks);
        k = peaks(:,1); p = peaks(:,2);
        tailPeaks = [k,p];
    %%%%%
    
    %%%%%Body Amplitude Finder
%     BodyAmps = [];
%     for i = 3:2:19
%         Y = (struct.Y(i,:));
%         [p1,k1] = findpeaks(Y,struct.t,'MinPeakProminence',5);
%         [p2,k2] = findpeaks(-Y,struct.t,'MinPeakProminence',5);
%         p = [p1';abs(p2')]; k = [k1';abs(k2')];
%         
%         % check if peakfinder did it's damn job
%         [k,p] = correctPeakFinder(k,p,struct.t,Y,['Body ',num2str(i),' peaks']);
%         peaks = [k,p]; peaks = sortrows(peaks);
%         k = peaks(:,1); p = peaks(:,2);
%         PositionPeaks = [k,p];
%         
%         p = polyfit(struct.t, Y,1);     % fit line for the tail wave
%         yT = p(1).*(PositionPeaks(:,1)) + p(2); % y values for that line
%         PositionAmps = PositionPeaks(:,2) - yT; % subtract those y values from that
%                                         % line to get actual amplitude
%         PositionAmps = abs(PositionAmps*struct.VidScale);
%         BodyAmps = [BodyAmps,median(PositionAmps)];
%         
%     end
    %%%%%   

    amps = tailPeaks(:,2);        
    tailAmps = abs(amps*struct.VidScale);
    wavenum = size(tailPeaks(:,2),1);

%     if mod(size(tailPeaks(:,1)),2) == 0
%         modEnd = 2;
%     else 
%         modEnd = 1;
%     end
%     wavelen = [];
%     for i = 1:2:size(tailPeaks(:,1),1)-modEnd
%         tailTime = tailPeaks(:,1);
%         wavelen = [wavelen; tailTime(i+2) - tailTime(i)];
%     end

    %%%%%%Wavelength calculator
%         [p1,k1] = findpeaks(y(20,:),'MinPeakProminence',5);
%         [p2,k2] = findpeaks(-y(20,:),'MinPeakProminence',5);
%         p = [p1';abs(p2')]; k = [k1';abs(k2')];
%         xValues = x(20,k);
%         % check if peakfinder did it's damn job
%         [xValues,p] = correctPeakFinder(xValues,p,x(20,:),y(20,:),['Wave peaks']);
%         wavePeaks = [xValues',p];
%         wavelen = [];
%         for i = 2:2:size(wavePeaks,1)
%             wavelen = [wavelen; ...
%                 pdist([wavePeaks(i-1,:);wavePeaks(i,:)])*struct.VidScale];
%         end
            
    %%%%%%%
    
    nose = [struct.midLines(1).MidLine(1,:);struct.midLines(end).MidLine(1,:)];
    distance = pdist(nose, 'euclidean');
    distance = distance.*struct.VidScale;                   % in mm
    struct.swimmingSpeed = (distance/struct.t(end));        % in mm/s
    struct.bendingFrequency = wavenum/2/struct.t(end);      % in hZ
    struct.bendingPeriod = 1/struct.bendingFrequency;       % in seconds
%     struct.wavelength = median(wavelen);
    struct.bendingStrideLength = distance/(wavenum/2);      % in mm
%     struct.bendingWS = struct.wavelength*struct.bendingFrequency;
    struct.bendingAmp = median(tailAmps);                   % in mm
    struct.bendingAmps = tailAmps;
%     struct.bodyAmps = BodyAmps;
    eval([cell2mat(vars(1)), '= struct'])
    save(cell2mat(vars(1)), cell2mat(vars(1)));
end

function [k,p] = correctPeakFinder(k,p,X,Y,Title)
    activeFigure = figure;
    plot(X,Y);
    hold on
    plot(k,p,'r*');
    title(Title);
    
    prompt = {'How Many False Peaks?', 'How Many Missing Peaks?'};
    BoxName = 'FindPeak Error correction';
    default = {'0','0'};
    answer = inputdlg(prompt, BoxName,1,default);
    answer = str2double(answer);
    % if the user needs to eliminate points
    if answer(1) ~= 0
        % to eliminate peaks
        for i = 1:answer(1)
            rect = getrect();
            elim = find(k>rect(1) & k<(rect(1)+rect(3)));
            k(elim) = []; p(elim) = [];
            plot(k,p,'o');
        end
    else
    end
    % if the user needs to specify points
    if answer(2) ~= 0
        [x, y] = getpts();
        p = [p;y];
        k = [k;x];
        plot(k,p,'o');
    else
    end
    close(activeFigure)
end
