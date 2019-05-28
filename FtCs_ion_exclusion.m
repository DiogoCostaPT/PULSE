

% Diogo Andre Pinho da Costa


%% EULER Explicit 

function c = FtCs_ion_exclusion(SWEi,tsteps,nlx,nly,D,c,deltx,delty,deltt,U,Katm)

atmi=[1:1:tsteps];
atmii=ceil(atmi*nly/tsteps);

dcdx=zeros(nlx,nly,tsteps);
dc2dx2=zeros(nlx,nly,tsteps);
dc2dy2=zeros(nlx,nly,tsteps);
dcdt=zeros(nlx,nly,tsteps);
laydepth=zeros(nlx,nly);

ldry=ones(nlx,nly);
ldry_prev=ones(nlx,nly);

upSnowlayer=zeros(tsteps); % variable that gets the upper snow cell
upSnowlayer_prev=zeros(tsteps); % variable that gets the upper snow cell

Acell=ones(nlx,nly)*deltx*delty;

h = waitbar(0,'Please wait...');

for t=2:tsteps  
      
     % check height of snowpack for numerical model
         upSnowlayer=ceil(SWEi(t)/delty);
         if upSnowlayer==0 % SWE =0
             c(:,:,t)=0.;
             continue
         else
             ldry(:,1:upSnowlayer)=0.;
             c(:,upSnowlayer+1:end,t)=0.;
            if upSnowlayer==1
                laydepth(:,2:end)=0;
                laydepth(:,1)=SWEi(t);
            else
                laydepth(:,1:upSnowlayer-1)=delty;
                laydepth(:,upSnowlayer)=SWEi(t)-(upSnowlayer-1)*delty;
            end
         end
    
       if SWEi(t)>=1
       c(2:end-1,upSnowlayer+1,t-1)=c(2:end-1,upSnowlayer+1,t-1)+ones(nlx-2,1)*Katm;      % mg/l
       end
       
    if upSnowlayer==1 % SWE =0
        continue
    end
       
    for yi=2:upSnowlayer+1; % Loop until the top of the snowpack
            
        for xi=2:nlx-1;         
         dcdy(xi,yi,t)=(c(xi,yi,t-1)-c(xi,yi-1,t-1))/delty;
         dc2dy2(xi,yi,t)=((c(xi,yi+1,t-1)-c(xi,yi,t-1))/delty-(c(xi,yi,t-1)-c(xi,yi-1,t-1))/delty)/(delty);         
         dc2dx2(xi,yi,t)=((c(xi+1,yi,t-1)-c(xi,yi,t-1))/deltx-(c(xi,yi,t-1)-c(xi-1,yi,t-1))/deltx)/(deltx);
         dcdt(xi,yi,t)=-U*dcdy(xi,yi,t)+D*(dc2dx2(xi,yi,t)+dc2dy2(xi,yi,t));
        end
         
        %end 
        
    end

    %dcdt(:,1,t)=0;
    %dcdt(:,end,t)=0;
    %dcdt(:,end-1,t)=0;
    
    c(:,:,t)=c(:,:,t-1)+dcdt(:,:,t)*deltt;
    
    
    % Von Neumann condition: dc/dx=0 in the river margins
    c(1,:,t)=c(2,:,t);  
    c(end,:,t)=c(end-1,:,t);   
    
    % Von Neumann condition: dc/dy=0 in the river BCs
    %c(:,1,t)=c(:,2,t);  
    %c(:,end,t)=c(:,end-1,t);   
    
    c(:,upSnowlayer+2:end,t)=0.;
    
    waitbar(t / tsteps)
    
end

close(h) 

