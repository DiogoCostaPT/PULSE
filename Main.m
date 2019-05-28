
function main

%% Inputs
% L - length to be analyzed, m
% B - width, m
% D - average depth, m
% deltx - to determine the number of layers, m - need to be multiple of L (to simplify for now)
% delty - to determine the number of layers, m - need to be multiple of B
% deltz - to determine the number of layers, m - need to be multiple of D
% M - mass of conservative contaminant, kg (e.g salt) - to be injected in
%        the river cell (assumption)
% CO - average concentration of the contaminant before the slug test
% s - position where the slug test takes palce, cross section cell (counting from the left margin)
% Tperd - simulation time interval, min
% deltt - time interval (s) - need to be multiple of t

% Outputs
% c - concentration, mg/l

 % Example 6: ion exclusion
D = 0       % m2/s
U = 0.005

Katm = 0.2     % kg


% example - lecture
B = 10        % m
H = 4         % m
L = 120      % m
deltx = 1    % m
delty = 1    % m
deltz = H    % m
C0 = 0        % mg/l
s = 12         % adm
%M = 300      % kg
%D = 0.009       % m2/s
Tperd = 90     % day 
deltt = 1    % h

cmaxlim = 10

A=B*H;


%%
tic

% Initializing
% 1) Defining the GRID
nly=L/delty;
nlx=B/deltx;
nlz=H/deltz;

if (nlx/round(nlx)|nly/round(nly)|nlz/round(nlz))~=1;
    display('Error: \Delta i, i=x,y,z, needs to be multiple of L,B and D, respectively')
    return    
end

% 4) Import SWE from CRHM
dirloc='D:\OneDrive\DI_PRF_CUR\UofS\7_Research_Sites\1_STC\4_CHRM_models\StepperTwin_STC';
CRHMrawALL=importdata([dirloc,'\CRHM_output_1.txt']);
timeCRHMall=datevec(datestr(CRHMrawALL.data(:,1)+693960));

SimSTART=[2008,10,1,1,0,0];
SimEND=[2009,09,30,0,0,0];

[CRHMraw timeCRHM]=GetCRHMresults(CRHMrawALL,timeCRHMall,SimSTART,SimEND);


% 1) Initializing the time steps
tsteps=numel(timeCRHM(:,1)); % hours

% 2) Initial concentration
c=ones(nlx,nly,tsteps);
c(:,:,:)=c*C0;
Vcell=deltx*delty*deltz;



HRUorderOutput=[1 2 3 4 5 6 7 8 9 10 11];
HRUsNo=11;
Years=[2008,2009,2010,2011];

i=8; % SWE -> used in Csnow and in Csoil (the later to determine when it's winter)
[colin colend]=findDataRaw(HRUsNo,i);
SWE=CRHMraw(:,colin:colend);
SWE=SWE(:,HRUorderOutput);

clear CRHMraw
clear CRHMrawALL

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
 surf(datenum(timeCRHM),[1:1:HRUsNo]',flipud(SWE'),'Parent',axes1,'EdgeColor','none','LineStyle','none','FaceLighting','phong');  
view(axes1,[-18 24.2]);
ylabel('hru')
xlabel('Time (month)')
datetick('x','mmm')
colorbar
grid(axes1,'on');
ylim([1 HRUsNo])

SWEi=SWE(:,4);


%% Calling numerical scheme
 c = FtCs_ion_exclusion(SWEi,tsteps,nlx,nly,D,c,deltx,delty,deltt,U,Katm);

toc
%% Plotting

% This loop is to plot the concentration change throught the test

figure(1)
filename = 'nada.gif';
% subplot(2,2,1)
% surf(c(:,:,1))
% drawnow
% 
% subplot(2,2,2)
% plot(c(s,:,1))
% drawnow
% 
% subplot(2,2,3)
% surf(c(:,:,1))
% view(90,0)
% drawnow
% 
% subplot(2,2,4)
% surf(c(:,:,1))
% view(0,90)
% drawnow




%cmaxlim = max(max(max(c)));



subplot(2,2,2)                  % Plant view
colormap(jet)
i2nan=find(c==0.);
c(i2nan)=NaN;

h1=surf(c(:,:,1))%,'EdgeColor','k','LineStyle','-','FaceLighting','phong');     
%shading interp   
caxis([0 cmaxlim])
ylim([1 B])
xlim([1 L])
grid on
ylabel('Z-Direction (snow depth) [m]')
ylabel('X-Direction [m]')
title('Initial Conditions: Snowpack ion concentrations (vertical profile)')
view(270,90)
drawnow; 
colorbar('location','eastoutside')


for ti=2:tsteps/24   
%for t=tsteps-1:tsteps 
%for t=1000:1500    
    
t=ti*24; % print per day

 %t=2 % this is the time step selected for       
         subplot(2,2,3)                 % 3D view
         colormap(jet)
         caxis([0 cmaxlim])
         shading interp
         surf(c(:,:,t-1),'EdgeColor','none','LineStyle','none','FaceLighting','phong');  
         caxis([0 cmaxlim])
         %axis equal
         %hold on
                  title('3D view')
         ylim([1 B])
         xlim([1 L])
        xlabel('Z-Direction (snow depth) [m]')
         ylabel('X-Direction [m]')
        zlabel('Concentration [mg/l]')
                title('3D view')
         zlim([0 cmaxlim])
         grid on
        %drawnow
        
        %colormap(hotflip)
        subplot(2,2,1)                  % Transverse view
        %plot(c(2,:,t),[1:1:nly])
        tplot=datenum(timeCRHM(1:t,:));
        plot(tplot,SWEi(1:t))
        %plotyy(tplot,SWE(1:t),tplot,c(ceil(end/2),:,t),
        %xlim([0 cmaxlim])
        ylim([1 max(SWEi)])
        xlim([datenum(timeCRHM(1,:)) datenum(timeCRHM(end,:))])
        datetick('x','yyyy-mm')
        %axis equal
        ylabel('Z-Direction (snow depth) [m]')
        xlabel('Time [day]')
        timeprint=['Date = ',num2str(timeCRHM(t,1:3)),' [yyyy mm dd]'];
         timestepprint= ['Timestep = ',num2str(t), ' h'];
           txt=[{'SWE [mm]',timeprint,timestepprint}];
           title(txt)
        drawnow;
        
% 
%         
%         subplot(2,3,2)                  % Cross-sectional view
%         colormap(jet)
%         surf(c(:,:,t-1),'EdgeColor','none','LineStyle','none','FaceLighting','phong');
%         caxis([0 cmaxlim])
%         shading interp
%          %colormap(hot)
%          zlim([0 cmaxlim])
%         ylabel('X-Direction [m]')
%         zlabel('Concentration [mg/l]')
%         title('Cross-sectional view')
%         view(90,0)
%         drawnow;   
        
        
        subplot(2,2,4)                  % Plant view
        colormap(jet)
        h2=surf(c(:,:,t-1))%,'EdgeColor','k','LineStyle','-','FaceLighting','phong');     
        shading interp   
        caxis([0 cmaxlim])
        ylim([1 B])
         xlim([1 L])
        grid on
        xlabel('Z-Direction (snow depth) [m]')
        ylabel('X-Direction [m]')
        title('Snowpack ion concentrations (vertical profile)')
        view(270,90)
        drawnow; 
        
%          
%frame = getframe(1);
%im = frame2im(frame);
%[imind,cm] = rgb2ind(im,256);

%    if t == 2;
 %   imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%    else
%    imwrite(imind,cm,filename,'gif','WriteMode','append');
 %   end
         
end
%end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary Functions

function [CRHMraw timeCRHM]=GetCRHMresults(CRHMrawALL,timeCRHMall,SimSTART,SimEND)
istart=find(datenum(timeCRHMall)==datenum(SimSTART));
iend=find(datenum(timeCRHMall)==datenum(SimEND));
CRHMraw=CRHMrawALL.data(istart:iend,:);
timeCRHM=timeCRHMall(istart:iend,:);
end

function [colin, colend]=findDataRaw(HRUsNo,i)
colin=1+(i-1)*HRUsNo+1;
colend=1+i*HRUsNo;
end
    