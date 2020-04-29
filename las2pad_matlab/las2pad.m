
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this note is reproduced on each copy made. This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%   For more details See Arnqvist et al. Robust processing of airborne laser scans to plant area density profiles, 2020, biogeosciences.
% Please send a notification of any bugs to johan.arnqvist[at]geo.uu.se

%%%----------------------Start of program-------------------------------

%Set folder to collect .las data and destination folder for the PAD data

pathget='C:\Users\user\Documents\folderpath\';
pathend='C:\Users\user\Documents\folderdestinationpath\';

%%%%-------------Define grid sizes and constants

dx=20; %Grid size in z;
dy=20; %Grid size in y;
dz=5; %Gridsize in z-direction
hmax=40; %The height of the tallest allowed trees
k=0.5; %Extinction coefficient
%%%--------------------------------------------------------------

%Look for las data 
cd(pathget);
nms=ls;
inds=[];
for i=1:length(nms(:,1))
    if ~isempty(strfind(nms(i,:),'.las'))
        inds=[inds; i];
    end
end
nms=nms(inds,:);
inds=1:size(nms,1);

%Start the conversion
tic;
for il=1:length(inds)
    ind=inds(il);
    nme=nms(ind,:);
    cd(pathget);
    [X Y Z C Int Rnum Num phi]=LAStoMatrix([nme(1:20) '.las']);

    %Scan the data set and scale the values according to the return number and intensity
    Isum=Int;
    rmax=max(Num);
    disp ([num2str(length(find(Num==1))./length(Isum)*100) ' % of total shots with Rmax = 1']) 
    for i=2:rmax
        p=find(Rnum==i & Num==i);
        maxnum=length(p);
        for ii=1:(i-1) %Check that there is actually the right number of returns preceeding Num=i
            p=intersect(p,find((Rnum(i:end)-Rnum((i:length(Rnum))-ii))==ii)+(i-1));
        end
        disp ([num2str(length(find(Num==i))./length(Isum)*100) ' % of total shots with Rmax = ' num2str(i) ', the percentage of points in order is: ' num2str(100*length(p)./maxnum)]);
            
        if Num(1)~=1 %If the start is from outside the file
            p=p(2:end);
        end
        sum1=Int(p); 
        for ii=1:(i-1) %Find the total intensity for each shot/ray
            sum1=sum1+Int(p-ii);
        end
        Isum(p)=sum1; %The sum should be the total intensity
        for ii=1:(i-1) %Set the intensity sum for all returns of a shot/ray
            Isum(p-ii)=sum1;
        end
    end  
    
    
    Iscale=Int./Isum; %This scales according to share of the total intensity. For points that do not have the order of returns sorted, the scale will be 1.
    
    
    p=~isnan(Iscale);    %Remove points that become NaN
    X=X(p);
    Y=Y(p);
    Z=Z(p);
    C=C(p);
    Int=Int(p);
    Rnum=Rnum(p);
    Num=Num(p);
    phi=phi(p);
    Iscale=Iscale(p);

    ll(1)=floor(min(X)); %Lower left corner
    ll(2)=floor(min(Y));
    ur(1)=ceil(max(X)); %Upper right
    ur(2)=ceil(max(Y));


    sx=ceil((ur(1)-ll(1))/dx); %size in x (in grid points)
    sy=ceil((ur(2)-ll(2))/dy); %size in y (in grid points)
    sz=length(dz:dz:hmax);%size in z (in grid points)
    
    %Allocate matrixes
    gh=zeros(sy,sx); %Ground height
    th=gh;  % Tree height
    PAI=gh; %PAI
    wtrflag=gh; %Flag to indicate if there is water
    PAD=zeros(sy,sx,sz); %PAD
    
    grp=find(C==2); %Ground hits
    wtp=find(C==9); %Water hits
    ngrp=find(C~=2); %Non ground hits 
    rn1=find(Rnum==1);%Return number =1
    rnh=find(Rnum~=1);%Return number is higher than 1
    
    for i=1:sx %Step along the grid points in x
        tstamp=toc; %Start a counter to keep track of time
        disp(['row ' num2str(i) ' @ ' num2str(tstamp) ' seconds'])
        p1=find(ceil((X-ll(1))./dx)==i);
        for ii=1:sy %Step along the grid points in y
            p2=find(ceil((Y(p1)-ll(2))./dy)==ii);
            p2=p1(p2);
            p2gr=intersect(p2,grp); %Ground hits
            wtp2=intersect(p2,wtp); %Water hits
            if ~isempty(p2gr) %There must be values enough to determine the ground height 
                %The height of the ground is determined through the median
                %of ground points in the grid box
                zg=median(Z(p2gr)); %For ground height use only values flagged as ground hit.
                gh(ii,i)=zg;
                Zp2=Z(p2)-zg; %Store in order to speed up.
                th(ii,i)=max(Zp2); %The tree height is taken as the maximum within the dx*dy square. 

                I1=zeros(sz+1,1); %Intensity vector in z
                I1(1)=sum(Iscale(p2gr));
                for i3=1:sz
                    I1(i3+1)=sum(Iscale(p2((Zp2)<(i3*dz)))); %All reflections counted equally (scale takes into account the weight in the average)
                end
                P=(I1)./(I1(end)); %Ratio of reflections to total number of reflections
                %Calculcate PAD and PAI with Beer-Lambert law
                PAD(ii,i,:)=-diff(-mean(abs(cosd((phi(p2)))),'omitnan').*log(P)./k./dz);               
                PAI(ii,i)=-mean(abs(cosd((phi(p2)))),'omitnan').*log(P(1))./k;
            elseif ~isempty(wtp2)
                gh(ii,i)=median(Z(intersect(p2,wtp)));
                wtrflag(ii,i)=1;
            else
                gh(ii,i)=NaN;
                th(ii,i)=NaN;
                wtrflag(ii,i)=NaN;
                PAI(ii,i)=NaN;
                PAD(ii,i,:)=NaN;
            end
        end
    end
    cd(pathend)
    %save files
    evalc(['save grid' num2str(ll(1)) '_' num2str(ll(2)) ' gh th PAI PAD wtrflag dx dy dz'])
    tstamp=toc;
    disp(['grid' num2str(ll(1)) '_' num2str(ll(2)) ' is done @ ' num2str(tstamp) ' seconds'])
end
