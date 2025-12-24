
%Courtesy: Original code by Prof. Anand T. N. C. Associate Professor in the Department of Mechanical Engineering at IIT Madras

clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Input Variables                               %
%========================================================================%
fileName=''     ;% Enter filename here: if using one file only 
                              %this file is processed. If using many files, 
                              %this file is read to determine directory
numFiles='all'      ;% process 1 file or all files: 'all' or 'one'
filesLike='/User/clonedDirectory/SampleImages/Image/*.tif'  ;% if numfiles = 'all' give a filter to choose the files in the directory
binSize=5           ; %size of bin while plotting histograms
scale=10         ;% Enter magnification in micron per pixel. Enter 1 to get value in pixels
preFilter='none'    ;% preprocessing filter to remove speckles or noise etc:
                     % 'none', 'avg' ,'median', 'peak'; defaults to 'none'
filtSize=[3 3]      ;% size of preprocessing filter (typically [3 3]
bitDepth=8          ;% 8 bit or 16 bit images - can modify code to get this 
                     % info from arrays itself. Currenlty all images are
                     % converted to 8 bit before region growing. Will be
                     % updated in future versions
numWorkers=1        ;% number of threads to run in parallel. Use 1 for serial. distrib_computing_toolbox required to run in parallel (>1). Preferences of Parallel Computing Toolbox should support the number
                                  % when numWorkers=1, if matlab takes a long time create a pool due to the parfor statement in the prog, go to
                                  % Home->Preferences->Parallel Computing Toobox and uncheck the "Automatically
                                  % create a parallel pool..." option
deletePool=0        ;% 0 to keep parallel threads open (e.g. if you plan to run the program again soon. 1 to close the parallel pool (e.g. when running remotely)
%%%%%%%%%%%%%%%%%%% Droplet sizing options %%%%%%%%%%%%%%%%%%
% The values here determine the accuracy of the technique

subtMin = 0         ;% 0 = don't subtract; 1 = subtract minimum of image
refImage= 0         ;% 1 if you are supplying a recorded image, 0 to
                     % calculate a reference image. Calculating is common
refImFile=' '       ;% name of the file containing the reference image (Required only if refImage=1)
refImCalc='slidMax' ;% 'slidMax' , 'maxOfImage' , Filter used to calculate the reference image
filtxy= [170 170 ]  ;% Dimension of the filter - required for sliding maximum                     % filter. The filter size should be about 1.5 times the 
                     % largest particle dimension in pixels. Setting inmo
                     % to 4 below will bring up an imtool which can be used
                     % to find out this value
inmo=0        ;% Interactive mode: options 0 to 4. Interactive mode 
                     % will show many images (>2) and wait for user input (=4): that slows down 
                     % the code
normalize=1;% Normalize image after subtracting background. 0 to 
                     % turn off. Default is 1: normalizing takes into
                     % account some of the difference in laser intensity
                     % accross the image
globThreshPerc=48   ;% Global threshold in percentage of maximum value:
                     % Increase if the image is noisy, reduce or increase if droplets 
                     % are being split into parts and recognized as more 
                     % than one droplet: if multiple green circles come up 
                     % in the same region in the processed image. Reduce if
                     % droplets are not circled in Figure 
fillHoles=0         ;% Fill holes due to dark regions in center of drop. 1=on; 0=off.
                     % Filling holes slows the code, but may be required for large droplets                    
lowThresPerc = 42;    %Percentage of maximum value in the region. This value
                     % should be lower than globThreshPerc. Very low values
                     % will make droplets merge together. Shown in blue in
                     % Figure 
highThresPerc = 70  ;% Percentage of maximum value in the region. Shown in 
                     % green/red in figure
maxLowDist = 200    ;% maximum distance of pixel from centroid while 
                     % enlarging low threshold area. Same value is also used
                     % for high threshold area. Choose this based on size
                     % of largest droplet expected. Low values will
                     % artificially reduce droplet size. High values
                     % increase processing time

%%%%%%%%%%%%%%%%%%%%%%%% Post-processing options %%%%%%%%%%%%%%%%%%%%%%%%%%
maxAreaRatio = 1.25  ;% Maximum allowed ratio of low Theshold are to high
                     % threshold area. High values will allow more
                     % defocussed droplets to be sized. Droplets which meet
                     % this criterion are shown with highThresPerc area
                     % green. Droplets which do not meet are shown with red
figsave=0           ;% 0 or 1. 1=save figure showing droplets and histograms
                     % for each image in subdir. Saving figures slows down the code
                     % A single histogram of all droplets from all images
                     % is saved by default
subdir='processed'  ;% name of subfolder in which proessed images are saved
subsub='histo'      ;% name of subsubdirectory in which histograms are saved
datfile='dropletdatafile' ;% name of file in current directory in which
                           % droplet size data is saved. Can be used to
                           % post-process with different options quickly
dropletData(10000,21) = 0; % Size of array which has droplet data: set the first 
                           % value to expected number of droplets in a directory
                           % Second value is 21
smallestAllowed = 2      ;% Smallest allowed area in pixels
ellipticityAllowed = 0.6 ; % Range: 0 to 1; Allowed deviation from circle; Values close to 1
                           % will ignore non-circular droplets
                           
                           
                           
                           
                           
                           
                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Code starts here                              %
%========================================================================%
writtenon=date;
warning('off', 'Images:initSize:adjustingMag'); % Tells matlab not to display warning about image being too large to display at 100%
disp('************** Starting Droplet Sizing Algorithm *******************');
save('parameters.mat','writtenon','preFilter','fileName','numFiles','filesLike' ...
    ,'filtSize','bitDepth','subtMin','refImage','refImFile','refImCalc' ...
    ,'filtxy','normalize','globThreshPerc','lowThresPerc','highThresPerc','maxLowDist');
disp('Processing parameters saved in parameters.mat file');
fid=fopen(datfile,'w');
time0=tic;

% Create a parallel pool if none exists
if numWorkers>1
    if license('test','distrib_computing_toolbox')
        if isempty(gcp('nocreate'))
            if verLessThan('matlab', '8.2')
                matlabpool numWorkers;
            else
                pool = parpool(numWorkers);
            end
        else
            pool=gcp;
            if pool.NumWorkers~=numWorkers
                if verLessThan('matlab', '8.2')
                    matlabpool(close);
                    matlabpool numWorkers;
                else
                    delete(pool);
                    pool = parpool(numWorkers);
                end
            end
        end
    else
        disp('Not able to find parallel license');
    end
end

if(inmo)
    disp('**** Running in interactive mode.****');
else
    disp('**** Running in batch mode. Figures will not be displayed ****');
end
if(~figsave)
        disp('**** No figures will be saved ****');
end
switch numFiles
    case {'all'}
        direc = which(fileName);         %better to use fullfile?
        mkdir([fileparts(direc)],subdir);
        fileList = dir([fileparts(direc) filesep filesLike]);
        fileNames = {fileList.name}';
        cd(subdir);
        mkdir(subsub);
        cd ..;
    otherwise
        direc = which(fileName);
        mkdir([fileparts(direc)],subdir);
        fileNames{1}=fileName;
        cd(subdir);
        mkdir(subsub);
        cd ..;
        disp('You have chosen only 1 file to process');
end
% Creating mask to check if areas are within maxLowDist 
distMask=false(maxLowDist*2); % creating square matrix
x0=maxLowDist;
y0=maxLowDist;
for xi=1:size(distMask,1)
    for yi=1:size(distMask,2)
        distan=sqrt( (xi-x0)^2 + (yi-y0)^2 );
        if(distan<maxLowDist)
            distMask(xi,yi)=true;
        end
    end
end
% imshow(distMask)

h=1

for index=1:length(fileNames)
    time1=tic;
    clear lig1 f;
    close all
    I=imread(fileNames{index});
    disp(['**** Image ' num2str(index) ': Processing file: ',cell2mat(fileNames(index)),' ****']);
    I=I(:,:,1);
    [Ix Iy]=size(I);
% Preprocessing
    preProcStart=tic;
    switch preFilter
        case {'avg'}
            h = fspecial('average', filtSize);
            I1=filter2(h, I);
        case {'median'}
            I1=medfilt2(I,filtSize);
            I1=I1(ceil(0+filtSize(1)/2):ceil(size(I1,1)-filtSize(1)/2), floor(0+filtSize(2)/2):floor(size(I1,2)-filtSize(2)/2)); 
        case {'peak'}
            disp('Peak filter not coded yet.. Stopping!')
            break
        otherwise
            if(index==1)
            disp('No preprocessing filter chosen. Continuing')
            end
            I1=I;
    end
    if subtMin==1
        I2=I2-min(min(I1));
    else
        I2=I1;
    end
   preProcend=toc(preProcStart);
   disp(['Took ' num2str(preProcend) ' s to preprocess file ']);

    % Reference Image
    refImgCalcStart=tic;
    if(refImage==1)
        refImg=imread(refImFile);
    else
        switch(refImCalc)
            case {'slidMax'}
                refImg=ordfilt2(I1,filtxy(1)*filtxy(2),ones(filtxy(1),filtxy(2)));
            case {'maxOfImage'}
                refImg=double(max(max(I1)))*ones(Ix,Iy);
                
            otherwise
                disp('No background image chosen.. Stopping!')
                break
        end
    end
    I3=refImg-I2;
    if (normalize~=0)
        I4=double(I3)./double(refImg)*(2^bitDepth-1);
    else
        I4=I3;
    end
    I4disp=I4;
    refImgCalcEnd=toc(refImgCalcStart);
    disp(['Took ' num2str(refImgCalcEnd) ' s to create/process reference image']);
    %Global Thresholding  %% would it be better to subtract min and then have
    %as percentage of Max values instead of bitdepth?
    I5=I4;
    globThresh=globThreshPerc/100*(2^bitDepth-1);
    % I5(I5<globThresh)=0;
    % Thresholding and filling holes
    bw=im2bw(uint8(I5),globThresh/(2^bitDepth -1));
    if(fillHoles)
        bwFill=imfill(bw,'holes');
        diffFill=bwFill-bw;
        bw=bwFill;
        props = regionprops(bwFill,I4,'BoundingBox','MaxIntensity');
        if(~isempty(props))
            temp=I4;
            for i=1:size(props,1)
                bbox=props(i).BoundingBox;
                temp(max(1,floor(bbox(2))):min(size(temp,2),ceil(bbox(2)+bbox(4))),max(1,floor(bbox(1))):min(size(temp,1),ceil(bbox(1)+bbox(3))))=props(i).MaxIntensity;
            end
            I4(find(diffFill))=temp(find(diffFill));
        end
        if(inmo>3)
            figure();
            imshowpair(uint8(I4),uint8(temp),'montage');
            imshowpair(uint8(255*bw),uint8(temp),'montage');
        end
    end
    
    % Plotting Figures
    if(inmo)
        if(inmo>3)
            imtool(I);
            disp('Press any key to continue');
            pause;
        end
        if(inmo>2)
            procfig=figure('visible','off');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize window
            ll=0.05;bb=0.07;ww=0.1;hh=0.1; % changing left, bottom, width, height of subplots below
            h=subplot(2,3,1); imshow(I);
            title('Original Image','fontsize',12,'fontweight','bold');
            p = get(h, 'pos'); p(1)=p(1)-ll;p(2)=p(2)-bb; p(3)=p(3)+ww; p(4)=p(4)+hh;
            set(h, 'pos', p);
            h=subplot(2,3,2); imshow(I2);
            title('Preprocessed Image','fontsize',12,'fontweight','bold');
            p = get(h, 'pos'); p(1)=p(1)-ll;p(2)=p(2)-bb; p(3)=p(3)+ww; p(4)=p(4)+hh;
            set(h, 'pos', p);
            h=subplot(2,3,3); imshow(refImg);
            title('Reference Image','fontsize',12,'fontweight','bold');
            p = get(h, 'pos'); p(1)=p(1)-ll;p(2)=p(2)-bb; p(3)=p(3)+ww; p(4)=p(4)+hh;
            set(h, 'pos', p);
            h=subplot(2,3,4); imshow(uint8(I3));
            title('Subtracted Image','fontsize',12,'fontweight','bold');
            p = get(h, 'pos'); p(1)=p(1)-ll;p(2)=p(2)-bb; p(3)=p(3)+ww; p(4)=p(4)+hh;
            set(h, 'pos', p);
            h=subplot(2,3,5); imshow(uint8(I4disp));
            title('Normalized Subtracted Image','fontsize',12,'fontweight','bold');
            p = get(h, 'pos'); p(1)=p(1)-ll;p(2)=p(2)-bb; p(3)=p(3)+ww; p(4)=p(4)+hh;
            set(h, 'pos', p);
            h=subplot(2,3,6); imshow(uint8(bw*255));
            if(fillHoles)
                title('Thresholded and Filled Image','fontsize',12,'fontweight','bold');
            else
                title('Thresholded Image','fontsize',12,'fontweight','bold');
            end
            p = get(h, 'pos'); p(1)=p(1)-ll;p(2)=p(2)-bb; p(3)=p(3)+ww; p(4)=p(4)+hh;
            set(h, 'pos', p);
            figure(procfig);
            drawnow;
        end
    end
    

            
    centres  = regionprops(bw, 'centroid');
    seeds=cat(1,centres.Centroid);
    I6=uint8(I4);
    lastFilledRow=find(any(dropletData,2),1,'last');
    if(isempty(lastFilledRow))
        lastFilledRow=0;
    end
    
    lowThresVal=lowThresPerc/100*max(max(I6));
    
    
    disp(['Finding low and high threshold areas for ' num2str(size(seeds,1)) ' blobs. This could take a while..']);
    funcstart=tic;
    maxLowDistSq=maxLowDist^2;
    totCol=size(dropletData,2);
    parfor numBlobs=lastFilledRow+1:lastFilledRow+size(seeds,1)  %%% parfor
        localVar=zeros(1,totCol);
        initPos=round([seeds(numBlobs-lastFilledRow,2),seeds(numBlobs-lastFilledRow,1)]);  %%% y value, x value 
        [lowThRegion]=modRegionGrowing(I6, initPos,lowThresVal,maxLowDist,1,distMask);
        lowThCentreProp = regionprops(lowThRegion,I6,'Area','BoundingBox','Centroid', 'Eccentricity','EquivDiameter','MaxIntensity');
        if(~isempty(lowThCentreProp))% If low theshold does not give an area, skipping calculations
            highThresVal=highThresPerc/100*lowThCentreProp.MaxIntensity;
            [highThRegion]=modRegionGrowing(I6, initPos,highThresVal,maxLowDist,1,distMask);
            highThCentreProp = regionprops(highThRegion,'Area','BoundingBox','Centroid', 'Eccentricity','EquivDiameter');
            values=struct2array(highThCentreProp);   % Checking if high threshold gives an area
        else
            values=[];
        end
        if(isempty(values))
            continue
        else            % If high theshold does not give an area, skipping writing that line
            
            localVar(1)=index;
            localVar(2)=lowThCentreProp.Area;
            localVar(3:6)=lowThCentreProp.BoundingBox;
            localVar(7:8)=lowThCentreProp.Centroid;
            localVar(9)=lowThCentreProp.Eccentricity;
            
            colNum=10;
            localVar(colNum+1)=highThCentreProp.Area;
            localVar(colNum+2:colNum+5)=highThCentreProp.BoundingBox;
            localVar(colNum+6:colNum+7)=highThCentreProp.Centroid;
            localVar(colNum+8)=highThCentreProp.Eccentricity;
            
            if(localVar(2)/localVar(11))<maxAreaRatio
                localVar(20)=1;
                localVar(21)= sqrt((localVar(2)+localVar(colNum+1))/2*4/pi); % Calculating equivalent dia of circle from average area
            end
            dropletData(numBlobs,:)=localVar;
        end
        
    end
   funcend=toc(funcstart);
    
    
    
%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Post-processing                               %
%========================================================================%
% The portion above does not need to be rerun while making changes below    
if(inmo||figsave)
        f=figure('visible','off');       

        set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Maximize window
        ll=0.05;bb=0.07;ww=0.1;hh=0.1; % changing left, bottom, width, height of subplots below
        h=subplot(1,1,1); imshow(I);
        
        hold('on');
        for ii=lastFilledRow+1:lastFilledRow+size(seeds,1)
            if(dropletData(ii,1)==index)
                if(dropletData(ii,20))
                    h=rectangle('Curvature', [1 1],'Position',[dropletData(ii,(3)),dropletData(ii,(4)),dropletData(ii,(5)),dropletData(ii,(6))],...
                        'EdgeColor','b','LineWidth',2) ;
                    h=rectangle('Curvature', [1 1],'Position',[dropletData(ii,(12)),dropletData(ii,(13)),dropletData(ii,(14)),dropletData(ii,(15))],...
                        'EdgeColor','g','LineWidth',2) ;
                    diams1(ii)=dropletData(ii,21);
                elseif any(dropletData(ii,:))
                    h=rectangle('Curvature', [1 1],'Position',[dropletData(ii,(3)),dropletData(ii,(4)),dropletData(ii,(5)),dropletData(ii,(6))],...
                        'EdgeColor','b','LineWidth',2) ;
                    h=rectangle('Curvature', [1 1],'Position',[dropletData(ii,(12)),dropletData(ii,(13)),dropletData(ii,(14)),dropletData(ii,(15))],...
                        'EdgeColor','r','LineWidth',2) ;
                    
                end
            end
        end
        if(figsave)
            title('Number of droplets vs. diameter in pixel. Only maxAreaRatio applied. Not all droplets shown here are finally counted','fontsize',12,'fontweight','bold');
            cd(subdir); saveas(f,fileNames{index}); cd ..;
        end
        if inmo
            figure(f)
            title('Number of droplets vs. diameter in pixel. Only maxAreaRatio applied. Not all droplets shown here are finally counted','fontsize',12,'fontweight','bold');

        end
        clear f;
    end
    

    if(inmo>1||figsave)
        ff = figure('visible','off');
        
        if(exist('diams1'))
            mm=max(diams1);
            lastx=ceil(mm/binSize)*binSize;
            xx=(0:binSize:lastx);
            N=histc(diams1,xx);
            bar(xx,N,'histc');
            h = findobj(gca,'Type','patch');
            set(h,'FaceColor',[0 .3 .5],'EdgeColor','k');
            title('Number of droplets vs. diameter in pixel. Only maxAreaRatio applied','fontsize',12,'fontweight','bold');
            xlim([0 lastx]);
        end
        if(figsave)
            cd(subdir); cd(subsub); saveas(ff,fileNames{index}); cd ..; cd ..;
        end
        if(inmo)
            figure(ff)
        end
    end
    timeend=toc(time1);
    disp(['Took ' num2str(timeend) ' s to process file ' cell2mat(fileNames(index)) '. ' num2str(funcend) ' s of this time was spent in the region growing routine']);
   
end

cd(subdir);
disp('**** Saving array dropletData with all droplet data in the chosen subdir ****');
save dropletData;
cd ..;
%%
dropletData = dropletData(any(dropletData,2),:); % Removing empty rows
dropletData = unique(dropletData,'rows'); % Removing duplicate rows
d32=zeros(length(fileNames),1);
d2=zeros(length(fileNames),1);
d3=zeros(length(fileNames),1);
dropsPerImg=zeros(length(fileNames),1);
allDiams=zeros(size(dropletData,1),1);
for jj=1:size(dropletData,1)
        if(dropletData(jj,2)/dropletData(jj,11))<maxAreaRatio  % Col 20==1
            if((dropletData(jj,2)+dropletData(jj,11))/2>=smallestAllowed)
                ellRat=(dropletData(jj,5)+dropletData(jj,14))/(dropletData(jj,6)+dropletData(jj,15));
                if ellRat>1
                    ellRat=1/ellRat;
                end
                if ellRat > ellipticityAllowed
                    dropletData(jj,20)=4;% accepted droplets
                    dropletData(jj,21)= sqrt((dropletData(jj,2)+dropletData(jj,11))/2*4/pi); % Calculating equivalent dia of circle from average area
                    allDiams(jj)=dropletData(jj,21)*scale;
                else
                    dropletData(jj,20)=3;% rejected based on ellipticity ratio
                    
                end
            else
                dropletData(jj,20)=2; % rejected based on smallest size allowed
            end
        end
        im=dropletData(jj,1);
        if dropletData(jj,20)==4
            d2(im)=d2(im)+allDiams(jj).^2;
            d3(im)=d3(im)+allDiams(jj).^3;
            dropsPerImg(im)=dropsPerImg(im)+1;
        end
end



numRejAreaRatio=sum(dropletData(:,20)==0);
numRejSmallSize=sum(dropletData(:,20)==2);
numRejEllipt=sum(dropletData(:,20)==3);
numAccepted=sum(dropletData(:,20)==4);
% Calculating SMD
d32=d3./d2;
d32all=sum(d3)/sum(d2);
disp(['**** Post-processed ' num2str(size(dropletData,1)) ' droplets. Rejected ' num2str(numRejAreaRatio) ' based on area ratio, '])
disp(['**** ' num2str(numRejSmallSize) ' as they did not meet smallest acceptable size, and ' num2str(numRejEllipt) ' as they were not circular enough. ****' ]);
disp(['**** Accepted ' num2str(numAccepted) ' droplets out of ' num2str(size(dropletData,1)) '. ****']);
disp(['**** Average SMD of these droplets from all images is ' num2str(d32all) ' micron ****']);
allDiams=allDiams(any(allDiams,2),:);

% Plotting histogram of all droplets
allDrops = figure('visible','off');
mm=max(allDiams);
lastx=ceil(mm/binSize)*binSize;
xx=(0:binSize:lastx);
N=histc(allDiams,xx);
bar(xx,N,'histc');
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 .3 .5],'EdgeColor','k');
title('Number of droplets vs. diameter in micron','fontsize',12,'fontweight','bold');
if(isempty(lastx))
    lastx=0.001;
end
xlim([0 lastx]);
cd(subdir); cd(subsub);  
saveas(allDrops,'Histogram.png');cd ..; cd ..;
if(inmo)
    figure(allDrops);
end
% Plotting histogram of statistics per image
figNumDrops=figure('visible','off');
xx=(1:length(fileNames));
[AX,H1,H2] = plotyy(xx,dropsPerImg, xx,d32,'bar','plot');
xt=[1:1:length(fileNames)];
set(AX,'xtick',xt);
%set(H1,'Marker','*','MarkerSize',8);
set(H2,'Marker','o','Color','m','MarkerSize',8,'MarkerFaceColor',[1 0 0]);
set(AX(2),'ycolor','m')
grid on;
xlabel(AX(1),'Image number') % label x-axis
ylabel(AX(1),'Number of droplets accepted in image') % label left y-axis
ylabel(AX(2),'SMD from image') % label right y-axis
title('Number of droplets accepted vs. Image number','fontsize',12,'fontweight','bold');
cd(subdir); cd(subsub);
saveas(figNumDrops,'numberOfDroplets.png');cd ..; cd ..;
if(inmo)
    figure(figNumDrops);
end

%%
% Writing unscaled SMD values to file
fprintf(fid,'%s\t',d32all);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%s\t',d32);
fprintf(fid,'\n');

% SMD
% eccentricity
% size filter
 
fclose(fid);
if(deletePool)
    if ~isempty(gcp('nocreate'))
        if(~exist('pool'))
            pool=gcp;
        end
        if verLessThan('matlab', '8.2')
            matlabpool(close);
        else
            delete(pool); % stop parallel pool
        end
    end
end
progEnd=toc(time0);
disp(['**** Processing Complete: Took ' num2str(progEnd) ' s to process ' num2str(index) ' files ****']);
   

