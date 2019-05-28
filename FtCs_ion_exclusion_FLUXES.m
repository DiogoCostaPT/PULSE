
function FtCs_ion_exclusion_FLUXES

% user input
Lh = 100;     % m
Lx = 3;        % m

dh=1;
dx=1;

Tperd = 90     % day 
deltt = 1    % h

katm=0.1;

hru=3;

SimSTART=[2008,10,1,1,0,0];
SimEND=[2009,09,30,0,0,0];


% Preliminar calc (don't change)
tsteps=(Tperd*24)/deltt; % hours


[timeCRHM SWEi]=Load_SWEi_CRHM(SimSTART,SimEND,hru) ;



ADEsolver(Lx,Lh,dx,dh,timeCRHM,SWEi,katm)


% Solver
 function ADEsolver(Lx,Lh,dx,dh,timeCRHM,SWEi,katm)
 
 qfcds=zeros(Lh);
 qfcdsW=zeros(Lh);
 csnow=zeros(Lx,Lh);
 cmax=zeros(Lx,Lh);
 cmin=zeros(Lx,Lh);
 h=zeros(Lx,Lh);
 h_prev=zeros(Lx,Lh);
 
 nx1=Lx/dx;
% nh1=Lh/dh;
 
fp = 0
fe = 0
fw = 0
fee= 0
fs= 0
fn= 0
fnn= 0
fem= 0
fnm= 0
cmax= 0
cmin= 0

pfce = 0
pfde= 0
qfcn= 0
qfdn = 0
pp = 0
qp = 0
area = 0
areae = 0
arean = 0
cvolrate = 0
cf = 0
cbilan = 0
cvolpot = 0
cvolrat = 0
ntp = 0
hp = 0
hne = 0
he = 0
hn = 0
hnue = 0
hnn = 0
hnew = 0
con = 0
itold = 1
ppw = 0
pps = 0 
wvch = 0
wvfws = 0
wvfen = 0
wvb = 0
swvb = 0
rat = 0
D_coef_x = 0
D_coef_y  = 0
 

nh=Lh/dh;
nx=Lx/dx;

tsteps=numel(timeCRHM);

ldry=zeros(Lx,Lh);
ldry_prev=zeros(Lx,Lh);

upSnowlayer=zeros(tsteps); % variable that gets the upper snow cell
upSnowlayer_prev=zeros(tsteps); % variable that gets the upper snow cell


%...    BOUNDARY CONDITIONS (snow-soil interface)
qfcds=zeros(nx,1);
pj = zeros(nh,1);

 for it=1:tsteps
 
% Initializing          
         pfe1=0;
         pfw1=0;
         qfn1=0;
         qfs1=0;
         dc1=0;
         
         % check height of snowpack for numerical model
         upSnowlayer=ceil(SWEi(it)/dh);
         if upSnowlayer==0 % SWE =0
             csnow(:,:)=0;
             continue
         else
             ldry(:,1:upSnowlayer)=1
            if upSnowlayer==1
                h(:,2:end)=0;
                h(:,1)=SWEi(it);
            else
                h(:,1:upSnowlayer-1)=dh;
                h(:,upSnowlayer)=SWEi(it)-(upSnowlayer-1)*dh;
            end
         end
         
%...    SET INITIAL CONDITIONS %DC
         if (it==1)
            for ih=1,nx
              for ix=1,nh
                    ldry_prev(ix,ih)=ldry(ix,ih);
                    if (h(ix,ih)>0)
                        csnow(ix,ih)=conc_c0;
                        h_prev(ix,ih)=h(ix,ih)  
                    else
                        csnow(ix,ih)=0;
                    end 
              end
            end
            return
         end

         
%...    ADJUST CONCENTRATION TO NEW DEPTH
        if(it>itold)
            itold=it;
                  for ih=1:nh
                    for ix=1:nx
                      hnew=h(ix,ih);
                      if(hnew>0)
                        %csnow(ix,ih)=csnow(ix,ih)*h_prev(ix,ih)/hnew;
                        %h_prev(ix,ih)=hnew
                      else
                        csnow(ix,ih)=0;
                      end
                    end
                  end
        end        

%...    POLLUTION SOURCES (atmospheric deposition)
        csnow(:,upSnowlayer)=csnow(ix,upSnowlayer)+katm;

        if upSnowlayer==1 % if only 1 cell don't apply the model
            continue
        end
            
%      for ih=1:nh
%         for ix=1:nx       
%             cmax(ix,ih)=max(csnow(ix-1,ih),csnow(ix+1,ih),csnow(ix,ih-1),csnow(ix,ih+1),csnow(ix,ih))
%             cmin(ix,ih)=min(csnow(ix-1,ih),csnow(ix+1,ih),csnow(ix,ih-1),csnow(ix,ih+1),csnow(ix,ih))  
%         end
%      end              
        
    
           
for ih=2:nh % the model start by calculating ih = 1 (and calculates all ix), and so on and so forth

%        2.1) initialization and BC
%          2.1.1) initializing
        %dx=dx0(ih)            % cell x-length
%              iprgx=1               % ???
        is=ih-1               % previous cell y-direction
        in=ih+1               % next cell y-direction
        inn=min(ih+2,nx1)     % next-next cell y-direction (limiter to make the model stop in a limit of the formain)
%        dxn=dx(in)           % ???
%     
%         2.1.2) BC (east flux into the formain - ix = 0)
        pfce=csnow(1,ih)*pj(0,ih)*dx     % convective flux at x = 0 ("c" refers to convective)
        hp=max(h(2,ih),hdry)                  % limiter to avoid instabilities (small values of h)
        he=max(h(3,ih),hdry)                  % limiter to avoid crashing 
        fp=csnow(1,ih)
        fe=csnow(2,ih)
        hne=0;sqrt(hp*nt(1,ih)*he*nt(2,ih))/sigc/abs(x(2)-x(1))*dx*D_coef_x  % [m3/s] - ?????
        pfde=0.                              % no diffusive flux over boundary
%              
        pfe1=pfce                             % now we can write the x-flux over the east boundary 
%
%         2.1.3) starts the x-direction calculation
        for ix=1:nx         
%            if (itmb.eq.10.and. ih.ge.63 .and. ix.ge.33) &      % ???       
          
            dn=dn0(ix)           % cell x-length
            iw=ix-1              % previous cell x-direction      
   		    ie=ix+1              % next cell x-direction
            iee=min(ix+2,nh1)    % next-next cell x-direction (limiter to make the model stop in a limit of the formain)       
            
 %           ppw=pp       % specific discharge at cell face (EAST)
 %           pps=qfcdsW(ix)
                      
         % ...       check if grid dry
        if(ldry(ix,ih))   % if true set everything to zero (flow2 model gives this statement)
            pfe1=0.
            qfcds(ix)=0
            qfcdsW(ix)=0
            if (MconcGW2SW(ix,ih,1,1)<=0)
                csnow(ix,ih)=0
            goto 90
            end
        end   
            
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % DC 
        % Calculation of Dispersion in rivers (%%%%%%%%%% WRONG)
%         if (D_coef.eq.0.) &
%            D_coef_x=phi_elder*h(ix,ih)*us(ix,ih)   % Lateral Dispersion (Elder, 1959) 
%            D_coef_y= 0.011*v(ix,ih)**2*B_average**2/(h(ix,ih)*us(ix,ih))
%        else         % defined by the user
            D_coef_x=D_coef
            D_coef_y=D_coef
%        end
       
    % ...       values in center of cells (area of the cells, water depth, flows, concentrations)
            area=arbase(ix,ih)      % area of the cell (being analyzed)
            areae=arbase(ie,ih)     % area of the adjacent cell - east (not yet calculated for the present time step)
            arean=arbase(ix,in)     % area of the adjacent cell - north
            ntp=nt(ix,ih)           % eddx viscosity [m/s]
            hp=h(ix,ih)             % water depth of the cell (being analyzed)
            he=h(ie,ih)             % water depth of the adjaent cell - east (not yet calculated for the present time step)
            hn=h(ix,in)             % water depth of the adjaent cell - west (not yet calculated for the present time step)
            fw=csnow(iw,ih)       % concentration at x-1 cell - west
            fp=csnow(ix,ih)       % concentration at cell (being analyzed)
            fe=csnow(ie,ih)       % concentration at x+1 cell - east
            fee=csnow(iee,ih)     % concentration at x+2 cell - east
            fs=csnow(ix,is)       % concentration at y-1 cell - south
            fn=csnow(ix,in)       % concentration at y+1 cell - north
            fnn=csnow(ix,inn)     % concentration at y+2 cell - north        


            % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 %.... limit the incoming fluxes to the possible arising from the increase/decrease in water level           
            
         
            pp=pj(ix,ih)
            qp=qj(ix,ih)
    
           
            
            
    %         2.1.4) determining the incoming at the cell in the X- and Y- direction (from the previous time-step)

    % ...       (X-direction) - fluxes over west face
            pfw1=pfe1       % convective + diffusive flux [m3/s] (mass balance) - coming from adjacent cells

    % ...       (Y-direction) - sum of diffusive+convective flux over south face 
            qfs1=qfcds(ix)             %[m3/s]

    %         2.1.5) determining the diffusion at the present time step X-direction (pfde, where d refers to diffusion)
 
    %   	x	location of x-lines	[m]
    %   	y	location of y-lines	[m]

    % ...       diffusive flux and mean concentration at east face
            if(ldry(ie,ih)==0)             % to check if the cell is not dry
                hnue=max(hp*nt(ix,ih)*he*nt(ie,ih),.0001)     % limiter
                hne=sqrt(hnue)/sigc/abs(x(ix)-x(ie))*dx*D_coef_x  % [m3/s] - x(ix)-x(ie) is the distance between the cell and the east-adjacent cell
                pfde=-hne*(fe-fp)                     % diffusive flux

                    if(pp>0)     % for the calculation of the advetive flux - the concentration average according to the direction of the flow (but forn't have experience with this method)
        %                fem=fp
                   %     if (.not.ldry(iw,ih)) &
                            fem=-.125*fw+.75*fp+.375*fe        % value at cell face
                   %      else
                    %    fem=0.5*fp+0.5*fe        % value at cell face
                    %     end
                         
                        else
                
                   %    if (.not.ldry(iee,ih)) &
                        fem=.375*fp+.75*fe-.125*fee
                   %     else
                    %     fem=0.5*fp+0.5*fe   
                    %    end
                     end

            else % if no flow, no concentration
            fem=0.
            pfde=0.
            end

%         2.1.6) if Boundary (overight the BC)
               
            fem=max(0.,fem)

            if(ix.eq.nh)   % if we are at the Boundary - overwrite with boundary condition
                fem=csnow(nh1,ih)
            end
   
    %         2.1.7) Calculating the advetive flux - X-direction     
    %            pfce=pj(ix,ih)*fem*dx          % [m3/s]
            pfce=pp*fem*dx          % [m3/s]

    %         2.1.8) Summing the advetive and diffusinve terms (X-direction)   
            pfe1=pfce+pfde                   % total flux = advective flux + diffusive

% ...       total flux at east face
%            pfce=pj(ix,ih)*fem*dx          % [m3/s]
            if(pfe1<0)         % inflwo from east cell
              if(ldry(ie,ih)==0)  % volume rate in east cell
                cvolrate=-(fe*he)*areae/dt_q 
                pfe1=max(pfe1,cvolrate*cflmb) %limit to available material
              else
                pfe1=0.
              end
            end
            
            
    %         2.1.5) determining the diffusion at the present time step Y-direction (pfde, where "d" refers to diffusion)

    %%%%%%%%%%%%%%%%%% DOUBTS (start) (AGAIN, the same forubt)
    % is "hne" the dispersion coefficient ?, in the case of the flow model it is the turbulent stress 
    % what is: "sigc" and "facdif"

            if(ldry(ix,in)==0)  
                hnue=max(.0001,hp*ntp*hn*nt(ix,in))
                hnn=sqrt(hnue)/sigc/dx*dn*D_coef_y              % [m3/s]
                qfdn=-hnn*(fn-fp)                    % diffusive flux
                if(qp>0)
                  %   if(.not.ldry(ix,is)) &
                    fnm=-.125*fs+.75*fp+.375*fn       % value at cell face
                  %   else
                  %   fnm=0.5*fp+.5*fn    
                  %   end
                 
                else
                %   if (.not.ldry(ix,inn)) &
                    fnm=.375*fp+.75*fn-.125*fnn
               % else
               %     fnm=.5*fp+.5*fn
               % end
                end
                
                
            else
                fnm=0.
                qfdn=0.
            end
            fnm=max(0.,fnm)

    % ???????????????????????????????????????????
    %%%%%%%%%%%%%%%%%% DOUBTS (end)

    %         2.1.6) if Boundary (overight the BC)
            if(ih.eq.nx) % overwrite with boundary condition
                fnm=csnow(ix,nx1)
            end

    %         2.1.7) Calculating the advetive flux - X-direction  
    %            qfcn=qj(ix,ij)*fnm*dn        % [m3/s]
            qfcn=qp*fnm*dn        % [m3/s]

    %         2.1.8) Summing the advetive and diffusinve terms (X-direction)   
            qfn1=qfcn+qfdn
            if(qfn1<0)              % inflwo from north cell
                if(ldry(ix,in)==0) 
        %          sourn=sour(i,ix,in)
                %cvolrate=-(fn*hn+sourn*dt_qes)*arean/dt_qes
                cvolrate=-(fn*hn)*arean/dt_q
                qfn1=max(qfn1,cvolrate*cflmb)      %limit to available material
                else
                qfn1=0.
                end
            end
      
    % ...       available volume rate [m3/s] in actual cell
            cvolpot=(fp*hp)*area % [m3] available material
            cvolrat=cvolpot/dt_q +(pfw1+qfs1) % inflow during actual time-step

    % ...       limit flux out of actual cell due to available material
            if (cvolrat>0)              % outflow is possible
                if(pfe1>0 & qfn1>0)   % both outflow
                if (pfe1+qfn1 > cvolrat)     % limit outflow to volrat
                    cf=qfn1/(pfe1+qfn1)
                    pfe1= (1.-cf)*cvolrat             % [m3/s]
                    qfn1=cf*cvolrat
                end
                elseif(pfe1 > 0)              % qfn1 is inflow
                pfe1=min(pfe1,(cvolrat-qfn1))
                elseif(qfn1>0)              % pfe1 is inflow
                qfn1=min(qfn1,(cvolrat-pfe1))         
                end
            else % bilance outflow with inflow
                if(pfe1>=0 & qfn1<0)  %restrict pfe1 to bilan
                cbilan=cvolrat-qfn1                
                if(cbilan>0)
                    pfe1=min(pfe1,cbilan)
                else
                    pfe1=0.
                end
                elseif(pfe1<0 & qfn1>=0) 
                cbilan=cvolrat-pfe1
                if(cbilan>0) 
                    qfn1=min(qfn1,cbilan)
                else
                    qfn1=0.
                end
                elseif(pfe1>=0. & qfn1>=0.)  
                pfe1=0.
                qfn1=0.
                end        
            end
 
            dc1=(pfw1-pfe1 + qfs1-qfn1)*dt_q/area  % [m]
            con=csnow(ix,ih)+dc1/h(ix,ih)           
               
           con=min(cmax(ix,ih),con)
           con=max(cmin(ix,ih),con)
           
%            if (pfw1>=0&qfs1>=0&pfe1<=0&qfn1>=0) con=0


            if (con<0.)  % DC - to be changed % limiter
                %write(nout, '(A23,I4,A7,I4,A11,f12.5,A6,f12.5,A5,f12.5,A5,f12.5,A6,f12.5,A7,f12.5,A7,f12.5,A7,f12.5,A7,f12.5,A8,f12.5,A8,f12.5,A7,f12.5,A7,f12.5,A8,f12.5,A8,f12.5,A6,f12.5,A6,f12.5,A6,f12.5,A6,f12.5)") 'Warning: C < 0 mg/l: ix= ',ix,' , ih= ',ih, ' ,t (sec)= ',tim,' ,C0= ',conc_SW(ix,ih,0),' ,C= ',con,' h= ', h(ix,ih),' h0= ', h0(ix,ih),' ,pfw= ', pfw1,' ,pfe= ', pfe1, ' ,qfs= ', qfs1, ' ,qfn= ', qfn1, ' ,pp(x)= ', pp, ' ,pq(y)= ', qp, ' ,ppw= ', ppw, ' ,pps= ', pps, ' ,u(x)= ', u(ix,ih), ' ,v(y)= ', v(ix,ih), ' ,he= ', h(ie,ih), ' ,hw= ', h(iw,ih), ' ,hs= ', h(ix,is), ' ,hn= ', h(ix,in)
                con=0
           elseif (con>4000) 
                %write(nout, '(A23,I4,A7,I4,A11,f12.5,A6,f12.5,A5,f12.5,A5,f12.5,A6,f12.5,A7,f12.5,A7,f12.5,A7,f12.5,A7,f12.5,A8,f12.5,A8,f12.5,A7,f12.5,A7,f12.5,A8,f12.5,A8,f12.5,A6,f12.5,A6,f12.5,A6,f12.5,A6,f12.5)") 'Warning: C > 2 mg/l: ix= ',ix,' , ih= ',ih, ' ,t (sec)= ',tim,' ,C0= ',conc_SW(ix,ih,0),' ,C= ',con,' h= ', h(ix,ih),' h0= ', h0(ix,ih),' ,pfw= ', pfw1,' ,pfe= ', pfe1, ' ,qfs= ', qfs1, ' ,qfn= ', qfn1, ' ,pp(x)= ', pp, ' ,pq(y)= ', qp, ' ,ppw= ', ppw, ' ,pps= ', pps, ' ,u(x)= ', u(ix,ih), ' ,v(y)= ', v(ix,ih), ' ,he= ', h(ie,ih), ' ,hw= ', h(iw,ih), ' ,hs= ', h(ix,is), ' ,hn= ', h(ix,in)
                %con=4000
           end  
            
            csnow(ix,ih)=con

            
    % ...       flux over south face
            qfcds(ix)=qfn1  % convective+diffusive flux


        end    %loop ix
%              if (ih=ihgr(iprgy))   
%                iprgy =iprgy+1
%              end

         end    %loop ih

    return
 end


 % Auxiliary function
function [timeCRHM SWEi]=Load_SWEi_CRHM(SimSTART,SimEND,hru)
         
         % 4) Import SWE from CRHM
dirloc='D:\OneDrive\DI_PRF_CUR\UofS\7_Research_Sites\1_STC\4_CHRM_models\StepperTwin_STC';
CRHMrawALL=importdata([dirloc,'\CRHM_output_1.txt']);
timeCRHMall=datevec(datestr(CRHMrawALL.data(:,1)+693960));

[CRHMraw timeCRHMvect]=GetCRHMresults(CRHMrawALL,timeCRHMall,SimSTART,SimEND);

HRUorderOutput=[1 2 3 4 5 6 7 8 9 10 11];
HRUsNo=11;
Years=[2008,2009,2010,2011];

i=8; % SWE -> used in Csnow and in Csoil (the later to determine when it's winter)
[colin colend]=findDataRaw(HRUsNo,i);
SWE=CRHMraw(:,colin:colend);
SWE=SWE(:,HRUorderOutput);

timeCRHM=datenum(timeCRHMvect);

clear CRHMraw
clear CRHMrawALL

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
 surf(timeCRHM,[1:1:HRUsNo]',flipud(SWE'),'Parent',axes1,'EdgeColor','none','LineStyle','none','FaceLighting','phong');  
view(axes1,[-18 24.2]);
ylabel('hru')
xlabel('Time (month)')
datetick('x','mmm')
colorbar
grid(axes1,'on');
ylim([1 HRUsNo])

SWEi=SWE(:,4);


function [CRHMraw timeCRHM]=GetCRHMresults(CRHMrawALL,timeCRHMall,SimSTART,SimEND)
istart=find(datenum(timeCRHMall)==datenum(SimSTART));
iend=find(datenum(timeCRHMall)==datenum(SimEND));
CRHMraw=CRHMrawALL.data(istart:iend,:);
timeCRHM=timeCRHMall(istart:iend,:);


function [colin, colend]=findDataRaw(HRUsNo,i)
colin=1+(i-1)*HRUsNo+1;
colend=1+i*HRUsNo;

