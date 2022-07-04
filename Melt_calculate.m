%% Input Data
%
% I prefer using a script over function because loading the image data
% takes a lot of time, so working iteratively is very costly.
%
% M - merged image of the whole ice core length of for M x N, where M is
% the width and N is the length dimension of the core
%
% x - corresponding depth (N X 1)
%
% mp - match points or Age Depth points (I X 2), where column 1 is depth
% and column 2 is corresponding year

%% Load data and accordanize
startyear=2016; % Age of core top
melt_cutoff_perc=50; % Cut-off percentage for considering an ice lens as melt layer
melt_cutoff=1400*melt_cutoff_perc/100; % Cut-off in number of pixels. My images were 1400 pixel wide, chnage it accoring to your images


if ~exist('M') | ~exist('x') | ~exist('mp') % check for all existing data
    load('MeltData50m.mat','M','x','mp');
end


if mp(1,1)~=0; % Core top check
    mp=[0 0; mp];
    mp(:,2)=startyear:-1:1+startyear-length(mp);
end


%% Start working
Melt=table; % Initialize a table for results
Melt.Poly=zeros(length(mp),1);
Melt.Thick=zeros(length(mp),1);
Melt.Year=mp(:,2);


hhh=waitbar(0);
for i=2:length(mp)
    waitbar(i/length(mp)-1,hhh,['Calculating Melt for ' num2str(mp(i,2))]) % Waitbar
    
    id=x>mp(i-1,1) & x<=mp(i,1); % Find section for the annual layer
    AnLayer=M(:,id); % crop to annual layer
    
    
    [B,L] = bwboundaries(AnLayer,'noholes'); % Find melt layers
    
    if numel(B)~=0 % check if melt layer exist
        
        PolyArea = regionprops(AnLayer, 'Area'); % Find area of melt layers
        PolyArea=[PolyArea.Area];
        
        ThickArea=NaN(length(B),1); % empty variable to thickness area
        lab=[];
        for k = 1:length(B)  % loop through each melt layer
            
            boundary =[]; % initialize parameters
            r=[];
            binaryImage=[];
            heighth=[];
            maxWidth=[];
            
            boundary = B{k}; % select each melt layer
            r=range(boundary(:,1)); % find vertical range to ensure 50% coverage
            
            if r<melt_cutoff
                lab(k) = 0;
                
            elseif r>=melt_cutoff & r<=1400
                
                try
                    binaryImage = AnLayer(:,min(boundary(:,2)):max(boundary(:,2))); % crop out melt layer with 5 pixel buffer on each side
                    [maxWidth, rowOfMaxWidth] = max(sum(binaryImage, 2));
                    height=r;
                    ThickArea(k)=height*maxWidth;
                catch
                    
                end
                
                lab(k)=1;
                
            else
                lab(k)=NaN;
                
            end
            
        end
        
        idx=lab==1;
        
        ThickArea(~idx)=[];
        PolyArea(~idx)=[];
        
        Melt.Poly(i)=nansum(PolyArea)/(1400*sum(id));
        Melt.Thick(i)=nansum(ThickArea)/(1400*sum(id));
        
    end
    
end

delete(hhh)