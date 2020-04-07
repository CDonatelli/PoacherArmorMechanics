function struct = wobbleWaveKineMod(structfile)
    load(structfile);
    clear 'structfile'
    vars = who;
    eval(['struct=',cell2mat(vars(1))])
    % struct.fishLength = 85; <-- for debugging
    % Initiating variables
        % reading variables from structure

    prompt = {'Enter fish length in mm (leave as 0 if unknown)'};
    dlgtitle = 'Fish Length';
    dims = [1 35];
    definput = {'0'};
    length = inputdlg(prompt,dlgtitle,dims,definput);
    
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
        
        if length ~= 0
            struct.fishLength = length;
            VidScale = length/fishPixels;
        else
            struct.fishLength = 0;
            VidScale = 0;
        end
        
        struct.VidScale = VidScale;

    perc = 0.2:0.1:0.9;
    struct.tailPts = length.*perc;
    tail = struct.tailPts ./length; %for interparc
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
        for j = 1:size(tail,2)
            cordinate = funct(tail(j));
            tailPtCordsX(j,i) = cordinate(1);
            tailPtCordsY(j,i) = cordinate(2);
        end
    end
    struct.tailPtCordsX = tailPtCordsX;
    struct.tailPtCordsY = tailPtCordsY;
    struct.X = x; struct.Y = y;
    % figure out time for each frame and make a vector of times
    fr = 120;
    total = nfr/fr; 
    struct.t = linspace(0,total,nfr);
    s = linspace(1,length,npts);
    struct.s = s';

%%%% 2D Wave Kinematics
    tailY = struct.Y(19,:);
    
    %%%%% Mini-Peak Finder
        [p1,k1] = findpeaks(tailY,struct.t,'MinPeakProminence',5);
        [p2,k2] = findpeaks(-tailY,struct.t,'MinPeakProminence',5);
        p = [p1';abs(p2')]; k = [k1';abs(k2')];
        
        % check if peakfinder did it's damn job
        [k,p] = correctPeakFinder(k,p,struct.t,tailY,['Tail peaks']);
        peaks = [k,p]; peaks = sortrows(peaks);
        k = peaks(:,1); p = peaks(:,2);
        tailPeaks = [k,p];
    %%%%%
    
    %%%%%Body Amplitude Finder
    BodyAmps = [];
    for i = 3:2:19
        Y = (struct.Y(i,:));
        [p1,k1] = findpeaks(Y,struct.t,'MinPeakProminence',5);
        [p2,k2] = findpeaks(-Y,struct.t,'MinPeakProminence',5);
        p = [p1';abs(p2')]; k = [k1';abs(k2')];
        
        % check if peakfinder did it's damn job
        [k,p] = correctPeakFinder(k,p,struct.t,Y,['Body ',num2str(i),' peaks']);
        peaks = [k,p]; peaks = sortrows(peaks);
        k = peaks(:,1); p = peaks(:,2);
        PositionPeaks = [k,p];
        
        p = polyfit(struct.t, Y,1);     % fit line for the tail wave
        yT = p(1).*(PositionPeaks(:,1)) + p(2); % y values for that line
        PositionAmps = PositionPeaks(:,2) - yT; % subtract those y values from that
                                        % line to get actual amplitude
        PositionAmps = abs(PositionAmps*struct.VidScale);
        BodyAmps = [BodyAmps,median(PositionAmps)];
        
    end
    %%%%%
    
    
    p = polyfit(struct.t, tailY,1);     % fit line for the tail wave
    yT = p(1).*(tailPeaks(:,1)) + p(2); % y values for that line
    amps = tailPeaks(:,2) - yT;         % subtract those y values from that
                                        % line to get actual amplitude
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
        [p1,k1] = findpeaks(y(20,:),'MinPeakProminence',5);
        [p2,k2] = findpeaks(-y(20,:),'MinPeakProminence',5);
        p = [p1';abs(p2')]; k = [k1';abs(k2')];
        xValues = x(20,k);
        % check if peakfinder did it's damn job
        [xValues,p] = correctPeakFinder(xValues,p,x(20,:),y(20,:),['Wave peaks']);
        wavePeaks = [xValues',p];
        wavelen = [];
        for i = 2:2:size(wavePeaks,1)
            wavelen = [wavelen; ...
                pdist([wavePeaks(i-1,:);wavePeaks(i,:)])*struct.VidScale];
        end
            
    %%%%%%%
    
    nose = [struct.midLines(1).MidLine(1,:);struct.midLines(end).MidLine(1,:)];
    distance = pdist(nose, 'euclidean');
    distance = distance.*struct.VidScale;
    %speed in m/s
    struct.swimmingSpeed = (distance/struct.t(end))/1000;
    struct.bendingFrequency = wavenum/2/struct.t(end);
    struct.bendingPeriod = 1/struct.bendingFrequency;
    struct.wavelength = median(wavelen);
    struct.bendingStrideLength = distance/(wavenum/2);
    struct.bendingWS = struct.wavelength*struct.bendingFrequency;
    struct.bendingAmp = median(tailAmps);
    struct.bendingAmps = tailAmps;
    struct.bodyAmps = BodyAmps;
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
