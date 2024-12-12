function [ cov] = Simulation3( step_size,time,rate,box_size,stick_prob)
%Simulation of spherical particle deposition
% %Input parameters:

% 1) step_size is the time step size of physics simulation (tested with values 
% from 0.001 to 0.01). Smaller values give better accuracy but takes more time
% 
% 2) time is the total simulated time
% 
% 3) rate is the number of generated particles per unit time
% 
% 4) box_size is the length of the box side with periodic boundary conditions
% 
% 5) stick_prob is the probability that during a collision of a moving particle
% with a deposit one, they stick together
% 


tic; %start timer


drawOn=0; %uncomment to draw each frame
calcCov=1; %1 to calculate surface coverage as a function of time

persistent particleRadius;
persistent xMax; %box corners
persistent xMin;
persistent yMin;
persistent yMax;
persistent N; %number of created particles
persistent movingIndices;
%coordinates of box corners
xMax=box_size;
xMin=-xMax;
yMax=xMax;
yMin=-yMax;

zMax=box_size; %height of the box

%data of particles:
positions=zeros(time*rate,3); %contains coordinates of every particle
deposed=zeros(time*rate); %contains information whether a particle was deposed onto surface
touchingInfo=zeros(time*rate,13); %contains information about particles that are in contact
notMovingLimit=3; %if a particle does not move this number of times in a row
%it is deposed to save computational time(only important when thermal motion is small) 

notMoving=zeros(1,time*rate); %contains how many times each particle does not move in a row

Nmoving=0;%number of moving particles

if(drawOn)
    figures=gobjects(time*rate); %contains graphical objects of spheres
end




%physical parameters
particleRadius=1e-6;
T=300; %K
kb=1.38e-23;
eta=8.9e-4;%viscocity
mass=1e-14;
dens=mass/(4/3*pi*particleRadius^3);
gforce=mass*(1-1000/dens)*[0,0,-9.8];%effective force due to gravity and buoyancy

stickProbability=stick_prob; %probability that two particles stick together after they collide

epsilon0=8.85e-12;
epsilonR=78;
charge=1.6e-19;
chargeConcentration=1e-1*6.02e26;%number density of ions in liquid
invDebLength=sqrt(charge^2*chargeConcentration/epsilon0/epsilonR/kb/T);
potential=30e-3;%surface potential of the particles

%constant of proportionality for electric double layer force:
elForceConst=64*pi*kb*T*chargeConcentration*particleRadius*potential^2/invDebLength;

ljs=particleRadius/5; %radius of Lennard Jones particles of which colloid particles are made of

Acc=1e-20; %Hamaker constant 
cutOffInnerEl=6.5e-3*particleRadius; %inner cut off distance for double layer force
cutOffInnerDisp=30^(-1/6)*ljs; %inner cut off distance for dispersion force
cutOffInnerLub=particleRadius*1e-3; %inner cut off distance for lubrication hydrodynamics
cutOffInnerBrown=30^(-1/6)*ljs; %inner cut off distance for lubrication brownian forces
%---------------

dirOfBrownF1=zeros(3,time*rate,time*rate); %contains information about pairwise brownian forces
dirOfBrownF2=zeros(3,time*rate); %contains the direction of isotropic brownian force for mid point algorithm


%spatial grid
%-----------------------------
%2d grid is used for finding nearest neighbours of each particle
nColumns=floor((xMax-xMin)/2/particleRadius); %number of columns in the grid
nRows=floor((yMax-yMin)/2/particleRadius);%number of rows in the grid
nLayers=floor(zMax/2/particleRadius); %maximum number of particles on top of each other
grid=zeros(nColumns,nRows,4*nLayers);
nInGrid=zeros(nColumns,nRows); %number of particles in each cell of the grid

inCells=zeros(4,2,time*rate); %every particle is in up to 4 cells of the grid
deltaX=(xMax-xMin)/nColumns; %size of cell
deltaY=(yMax-yMin)/nRows;


%-----------------------------


N=0; %number of created particles



simultColl=zeros(1,3); %contains information collision in the case of more than 2 particles colliding
simultShift=zeros(3,3);

%drawing figure in which spheres are drawn
%-------------------------
figure();
hold on;
plane1=[xMin,xMax,xMax,xMin];
plane2=[yMin,yMin,yMax,yMax];
plane3=[0,0,0,0];

fill3(plane1,plane2,plane3,'r');
alpha 0.1;
axis equal;
axis([2*xMin,2*xMax,2*yMin,2*yMax,0,zMax]);
xlabel('x');
ylabel('y');
%------------------------
cov=zeros(1,time); %surface coverage as a function of time

%the coverage is plotted against dimensionless reduced time,tau, which is  
%the time weighted by the number density of incoming particles
tau=pi*particleRadius^2/(xMax-xMin)/(yMax-yMin)*rate*(1:time); 
interval=1/rate;
nIntervals=time*rate;



for t=1:nIntervals
    
    
    N=N+1;
    
    Nmoving=Nmovig+1;
    movingIndices(Nmoving)=t;
    %physics simulation
    
    %making sure that a particle does not spawn inside other particle
    %-------------------------
    wrongPos=1;
    while(wrongPos)
        wrongPos=0;
        positions(t,:)=[xMin+rand()*(xMax-xMin);yMin+rand()*(yMax-yMin);zMax];
        cells=getCells(positions(t,:));
         for cel=1:4
             col=cells(cel,1);
             ro=cells(cel,2);
             
             tested=zeros(64,1);
             nTested=0;
             for it=1:nInGrid(col,ro)
                 p=grid(col,ro,it);
                 deltaR=positions(p,:)-positions(t,:);
                 deltaR=deltaR+getShift(deltaR);
                 if(abs(deltaR(3))<2*particleRadius && norm(deltaR)<2*particleRadius && ~any(tested(1:nTested)==p))
                     nTested=nTested+1;
                     tested(nTested)=p;
                     wrongPos=1;
                     break;
                 end
             end
             if(wrongPos)
                 break;
             end
         end
    end
    %--------------------------
    %adding particle to the grid
    inCells(:,:,N)=cells;
    for cel=1:4
         addToGrid(inCells(cel,1,N),inCells(cel,2,N),N);
    end
    
    
    physicsStep(interval);
   
    

    %calculates coverage
    if (mod(t,rate)==0 && calcCov)
        cov(t/rate)=coverage1();
    end
    fprintf('%.1f percent completed\n',t/nIntervals*100);
end

  
if(~drawOn)
    %draws all particles at the end of simulation if every frame is not
    %drawn
    drawAll();
end
hold off;
if (calcCov)
    %plots coverage
    xlabel('Reduced time');
    ylabel('Coverage');
    figure();
    plot(tau,cov);
end

fprintf('time %.3f s\n',toc);



  function addToGrid(column,row,i)
        %adds particle i to grid cell at given column and row
        nInGrid(column,row)=nInGrid(column,row)+1;
        grid(column,row,nInGrid(column,row))=i;
    end
    function removeFromGrid(column,row,i)
        %removes particle i from grid cell at given column and row
        sizeOfGrid=size(grid);
        indexArr=reshape(grid(column,row,:),[1,sizeOfGrid(3)]);
        indexArr=[indexArr(indexArr~=i),0];
        grid(column,row,:)=reshape(indexArr,[1,1,sizeOfGrid(3)]);
        nInGrid(column,row)=nInGrid(column,row)-1;
    end
    function addTouching(i,j)
        % adds information that particles are in contact
         touchingInfo(i,1)=touchingInfo(i,1)+1;
         touchingInfo(i,touchingInfo(i,1)+1)=j;
    end
    function removeTouching(i,j)
          % removes information that particles are in contact
        row=touchingInfo(i,2:end);
        row=row(row~=j);
        touchingInfo(i,2:end)=[row, 0];
        touchingInfo(i,1)=touchingInfo(i,1)-1;
    end
    function drawAll()
        %draws all spheres
        for i=1:N
            [x,y,z]=sphere;
            if(positions(i,3)<=particleRadius)
                c=[0 0 1];
            elseif(deposed(i))
                c=[0 1 0];
            else
                c=[1 0 0];
            end
            fig=surf(x*particleRadius+positions(i,1),y*particleRadius+positions(i,2),z*particleRadius+positions(i,3));
            set(fig, 'FaceColor', c);
            figures(i)=fig;
        end
        hold off;
    end




    function shift=periodicBoundary(position)
        %periodic boundary conditions for x and y directions - if a
        %particle leaves the box it is moved to the opposite side of the
        %box
        sideX=xMax-xMin;
        sideY=yMax-yMin;
        shift=[0,0,0];
        if(xMin>position(1))
            shift(1)=sideX*ceil(abs(xMin-position(1))/sideX);
        elseif (xMax<position(1))
            shift(1)=-sideX*ceil(abs(-xMax+position(1))/sideX);
        end
        
        if(yMin>position(2))
            shift(2)=sideY*ceil(abs(yMin-position(2))/sideY);
        elseif (yMax<position(2))
            shift(2)=-sideY*ceil(abs(-yMax+position(2))/sideY);
        end
        
    end
    function newMove=periodicShiftMove(i,move)
        %if a move of a particle is too large, it is reduced so that the
        %particle is always in the box
        newMove=move;
        sideX=xMax-xMin;
        sideY=yMax-yMin;
        if(positions(i,1)+move(1)>xMax+sideX)
            newMove(1)=move(1)-floor((positions(i,1)+move(1)-xMax)/sideX)*sideX;
        elseif(positions(i,1)+move(1)<xMin-sideX)
            newMove(1)=move(1)+floor((-positions(i,1)-move(1)+xMin)/sideX)*sideX;
        end
          if(positions(i,2)+move(2)>yMax+sideY)
            newMove(2)=move(2)-floor((positions(i,2)+move(2)-yMax)/sideY)*sideY;
        elseif(positions(i,2)+move(2)<yMin-sideY)
            newMove(2)=move(2)+floor((-positions(i,2)-move(2)+yMin)/sideY)*sideY;
          end
          if(positions(i,3)+move(3)>zMax*2)
              newMove(3)=move(3)-floor((positions(i,3)+move(3)-zMax)/zMax)*zMax+2*particleRadius;
              
          end
    end
    function c=coverage1()
        %calculates the fraction of surface area occupied by adsorbed particles 
        totalArea=(xMax-xMin)*(yMax-yMin);
        coveredArea=0;
        for i=1:N
            if (positions(i,3)<=particleRadius)
                coveredArea=coveredArea+pi*particleRadius^2;
            end
        end
        c=coveredArea/totalArea;
    end




%     function c=coverage(numbOfSteps)
%         totalArea=(xMax-xMin)*(yMax-yMin);
%         stepX=(xMax-xMin)/numbOfSteps;
%         stepY=(yMax-yMin)/numbOfSteps;
%         coveredArea=totalArea;
%         % g=zeros(numbOfSteps,numbOfSteps);
%         parfor i=1:numbOfSteps
%             
%             for j=1:numbOfSteps
%                 free=1;
%                 
%                 for par=1:N
%                     w=0;
%                     point=[xMin+i*stepX,yMin+j*stepY];
%                     if (norm([positions(par,1)-point(1);positions(par,2)-point(2)])<particleRadius)
%                         w=1;
%                     end
%                     if(w && deposed(par))
%                         free=0;
%                         break;
%                         %  g(mod(i,numbOfSteps)+1,floor(i/numbOfSteps)+1)=1;
%                     end
%                 end
%                 if(free)
%                     coveredArea=coveredArea-stepX*stepY;
%                 end
%             end
%         end
%         c=coveredArea/totalArea;
%         %         figure();
%         %         hold on;
%         %         for x=1:numbOfSteps
%         %             for y=1:numbOfSteps
%         %                 if (g(x,y)==1)
%         %                     fill3([(x-1)*stepX,x*stepX,x*stepX,(x-1)*stepX],[(y-1)*stepY,(y-1)*stepY,y*stepY,y*stepY],[0,0,0,0],'b');
%         %                 end
%         %             end
%         %         end
%         %         hold off;
%     end

    function fig=drawSphere(p,c)
        %draws a sphere of color specified by c 1x3 array at position of
        %p th particle
        if(~isempty(figures))
            try
                delete(figures(p));
            catch
            end
        end
        [x,y,z]=sphere;
        if (any(c))
            fig=surf(x*particleRadius+positions(p,1),y*particleRadius+positions(p,2),z*particleRadius+positions(p,3));
            set(fig, 'FaceColor', c);
        else
            fig=surf(x*particleRadius+positions(p,1),y*particleRadius+positions(p,2),z*particleRadius+positions(p,3));
            
        end
        
        drawnow;
        
        
    end

function shift=getShift(deltaR)
    %returns a periodic shift when two particles are in contact because of
    %periodic boundary conditions
    shift=zeros(1,3);
    if(deltaR(1)>deltaX*2.5)
        shift(1)=-xMax+xMin;
    elseif (deltaR(1)<-deltaX*2.5)
        shift(1)=xMax-xMin;
    end
    if(deltaR(2)>deltaY*2.5)
        shift(2)=-yMax+yMin;
    elseif (deltaR(2)<-deltaY*2.5)
        shift(2)=yMax-yMin;
    end
 end
    function v=volumeFraction(z,thickness)
        totalVol=(xMax-xMin)*(yMax-yMin)*thickness;
        occVol=0;
        for i=1:N
            if(deposed(i) && positions(i,3)<z+thickness+particleRadius  && positions(i,3)>z-particleRadius) %within the box
                h1=thickness-positions(i,3)+z;
                h2=positions(i,3)-z;
                capVol1=sphCapVolume(particleRadius,h1);
                capVol2=sphCapVolume(particleRadius,h2);
                occVol=occVol+4*pi/3*particleRadius^3-capVol1-capVol2;
            end
        end
        v=occVol/totalVol;
        %-------------------------------------
        %nested function to calculate the volume of spherical cap
        function v=sphCapVolume(R,h)
            if(h<R)
                v=pi*(2/3*R^3-R^2*h+h^3/3);
            else
                v=0;
            end
        end
        %-----------------------------------------
    end
    function f=dispElForce(r)
        %dispersion electric force for a given particle separation, r
        
        fDisp1=-32/3*Acc*particleRadius^6/r^3/(r^2-4*particleRadius^2)^2; %r^-6 term of LJ potential
     %   fDisp2=Acc/37800*ljs^6*((56*particleRadius*r^2+10*particleRadius^3-6*r^3)/(r-2*particleRadius)^8/r^2-(r^3+26*particleRadius*r^2+108*particleRadius^3+432*particleRadius^2*r)/(r+2*particleRadius)^8/r^2+12*(r^2-20*particleRadius^2)/r^9);
        f=fDisp1;
        %uncomment to include r^-12 term of LJ potential:
    %   fDisp2=Acc/37800*ljs^6*((56*particleRadius*r^2+10*particleRadius^3-6*r^3)/(r-2*particleRadius)^8/r^2-(r^3+26*particleRadius*r^2+108*particleRadius^3+432*particleRadius^2*r)/(r+2*particleRadius)^8/r^2+12*(r^2-20*particleRadius^2)/r^9);
        %f=f+fDisp2;
    end
    function w=smoothingFunc(x,b)
        %smoothing function to remove discontinuities due to cut offs
        w=1-(1-x).^b;
    end
 

    function [row,force]=getRMatrixRow(moveInd,keepDir)
        %returns a row of resistance matrix for particle of index moveInd
        %and total non-hydrodynamic force acting on it
        %keepDir is used to store random direction for mid-point algorithm
       i=movingIndices(moveInd);
        f=[0 0 0];
        row=zeros(3,3*Nmoving);
        %isotropic part of Brownian force
        if(any(dirOfBrownF2(:,i)))
            fBrown=dirOfBrownF2(:,i);
            dirOfBrownF2(:,i)=[0;0;0];
        else
            fBrown=randn(3,1);
            dirOfBrownF2(:,i)=fBrown;
        end
        
        
        
        ASumMatrix=zeros(3);
        cells=get16Cells(i);
        added=zeros(N,1);
        
        %the outer cut off of interactions is 4*particleRadius so all
        %nearest neighbours should be in 16 closest cells
        for c=1:16
            
            for itr=1:nInGrid(cells(c,1),cells(c,2))
                j=grid(cells(c,1),cells(c,2),itr);
                
                if(i~=j && abs(positions(i,3)-positions(j,3))<=4*particleRadius)
                    if(~added(j)) %storing particles that already added to the interactions
                  
                    added(j)=1;
                    deltaR=positions(j,:)-positions(i,:);
                     
                    %interactions through periodic boundaries
                    if(deltaR(1)>deltaX*2.5)
                        deltaR(1)=deltaR(1)-xMax+xMin;
                    elseif (deltaR(1)<-deltaX*2.5)
                        deltaR(1)=deltaR(1)+xMax-xMin;
                    end
                    if(deltaR(2)>deltaY*2.5)
                        deltaR(2)=deltaR(2)-yMax+yMin;
                    elseif (deltaR(2)<-deltaY*2.5)
                        deltaR(2)=deltaR(2)+yMax-yMin;
                    end

                    
                    r=norm(deltaR); %distance between particle centers
                    if(r<=4*particleRadius)
                        
                        
                        
                        
                        %fast lubrication dynamics (FLD)
                        %------------------------
                        dir=deltaR/r;
                        if(r<3*particleRadius)
                            del=r-2*particleRadius;
                             
                            if(del>cutOffInnerLub)
                                
                                Xpw=particleRadius/4/del+9/40*log(particleRadius/del);
                                Ypw=log(particleRadius/del)/6;
                                smooth=(smoothingFunc((3-r/particleRadius),4))^2;
                                Xpw=Xpw*smooth;
                                Ypw=Ypw*smooth;
                            else
                                
                                Xpw=particleRadius/4/cutOffInnerLub+9/40*log(particleRadius/cutOffInnerLub);
                                Ypw=log(particleRadius/cutOffInnerLub)/6;

                                
                            end
                            D=dir'*dir;
                            I=eye(3);
                            Apw=Xpw*D+Ypw*(I-D);
                            if(~deposed(j))
                                k=find(movingIndices==j);
                                row(:,k*3-2)=-Apw(:,1);
                                row(:,k*3-1)=-Apw(:,2);
                                row(:,k*3)=-Apw(:,3);
                            end
                            ASumMatrix=ASumMatrix+Apw;
                          
                            if(any(dirOfBrownF1(:,i,j)))
                                randVec=dirOfBrownF1(:,i,j);
                                if(~keepDir)
                                    dirOfBrownF1(:,i,j)=[0;0;0];
                                end
                            else
                                randVec=randn(3,1);
                                dirOfBrownF1(:,j,i)=-randVec;
                                dirOfBrownF1(:,i,j)=randVec;
                            end
                            
                            if(del>cutOffInnerBrown)
                                fBrown=fBrown+(sqrt(Xpw)*D+sqrt(Ypw)*(I-D))*randVec;
                            else
                                Xpw=particleRadius/4/cutOffInnerBrown+9/40*log(particleRadius/cutOffInnerBrown);
                                Ypw=log(particleRadius/cutOffInnerBrown)/6;
                                fBrown=fBrown+(sqrt(Xpw)*D+sqrt(Ypw)*(I-D))*randVec;
                            end
                        end
                         %------------------------
                         %dispersion force
                        if(r>2*particleRadius+cutOffInnerDisp) 
                            fDisp=dispElForce(r);
                            
                           
                        else
                            fDisp=dispElForce(2*particleRadius+cutOffInnerDisp);
                         
                        end
                        %double layer force
                        if(r>2*particleRadius+cutOffInnerEl) 
                          
                            fel=elForceConst*exp(-invDebLength*(r-2*particleRadius));
                        else
                            
                            fel=elForceConst*exp(-invDebLength*(cutOffInnerEl));
                        end
                        f=f-(fel+fDisp)*dir; %total electric force
                    end
                end
                end
            end
        end
        %in isotropic part of resistance matrix viscosity is a function of
        %volume fraction occupied by particles
        volFrac=N*4/3*pi*particleRadius^3/(xMax-xMin)/(yMax-yMin)/zMax;
        RCorrection=1+2.725*volFrac-6.583*volFrac^2;   
        
        %matrix row
        row(:,3*moveInd-2)=ASumMatrix(:,1)+[1;0;0]*RCorrection;
        row(:,3*moveInd-1)=ASumMatrix(:,2)+[0;1;0]*RCorrection;
        row(:,3*moveInd)=ASumMatrix(:,3)+[0;0;1]*RCorrection;
        row=row*6*pi*eta*particleRadius;
    
        %total force
        force=f'+fBrown*sqrt(12*pi*eta*particleRadius*kb*T/step_size)+gforce';
    end

    function u=getVelocities(keepDir)
        %returns vector containing all velocities of moving particles
        %uses FLD technique
        forces=zeros(3*Nmoving,1);
        RMatrix=zeros(3*Nmoving);
        
        for i=1:Nmoving
            
            [row,force]=getRMatrixRow(i,keepDir);
            forces(3*i-2)=force(1);
            forces(3*i-1)=force(2);
            forces(3*i)=force(3);
            RMatrix(3*i-2,:)=row(1,:);
            RMatrix(3*i-1,:)=row(2,:);
            RMatrix(3*i,:)=row(3,:);
        end
        u=RMatrix\forces;
    end
function shift=getShiftDeltaR(deltaR)
    %returns shift for particles in contact through periodic boundary
    shift=[0 0 0];
    if(deltaR(1)>2.1*particleRadius)
        shift(1)=-xMax+xMin;
    elseif (deltaR(1)<-2.1*particleRadius)
        shift(1)=xMax-xMin;
    end
    if(deltaR(2)>2.1*particleRadius)
        shift(2)=-yMax+yMin;
    elseif (deltaR(2)<-2.1*particleRadius)
        shift(2)=yMax-yMin;
    end
end
  function newVel=fixTouching(i,vel)
        %resting contacts algorithm for particles that are in contact at
        %the beginning of physics step
        
        %it prevents paritcle overlaps when they are in contact
        
        % uses Fast Contact Force Computation for Nonpenetrating Rigid Bodies
        %algorithm by David Baraff
        newVel=vel;
       
       n=touchingInfo(i,1);
       keys=touchingInfo(i,2:n+1);
       impulses=zeros(n,1);
       contactSpeed=zeros(n,1);
       contactNorms=zeros(n,3);
       collNotSolved=0;
       A=zeros(n);
       indices=[1:n];
       C=[];
       NC=[];
       d0=0;
       for c=1:n
           dir=normr(positions(keys(c),:)-positions(i,:)+getShiftDeltaR(positions(keys(c),:)-positions(i,:)));
           contactNorms(c,:)=-dir;
           contactSpeed(c)=-dir*vel;
           if(contactSpeed(c)<-1e-14)
                collNotSolved=1;
                if(d0==0)
                    d0=c;
                end
           

           end    
       end
       for a=1:n
           for b=a:n
               A(a,b)=contactNorms(a,:)*contactNorms(b,:)';
               A(b,a)=contactNorms(a,:)*contactNorms(b,:)';
           end
       end
       while (collNotSolved)
           collNotSolved=0;
           driveToZero(d0);
           newVel=vel;
           for c=1:n
               
               if(contactSpeed(c)<-1e-14 && ~any(C(C==c)))
                   collNotSolved=1;
                   d0=c;
                  break;
               end
               newVel=newVel+impulses(c)*contactNorms(c,:)';
           end
       end
       
        function driveToZero(d)
            deltaImpulse=fdirection(d);
            deltaContSpeed=A*deltaImpulse;
            [s,j]=maxStep(impulses,contactSpeed,deltaImpulse,deltaContSpeed,d);
            
            impulses=impulses+s*deltaImpulse;
            contactSpeed=contactSpeed+s*deltaContSpeed;
            
            if (any(C(C==j)))
                C=C(C~=j);
                NC(length(NC)+1)=j;
                driveToZero(d);
            elseif(any(NC(NC==j)))
                NC=NC(NC~=j);
                C(length(C)+1)=j;
                driveToZero(d);
            else
                if(d~=j)
                    error('error\n');
                end
                C(length(C)+1)=j;
            end
        end
        function deltaF=fdirection(d)
           deltaF=zeros(n,1);
           deltaF(d)=1;
           C=sort(C);
           if(~isempty(C))
           indsToRemove=indices;
           indsToRemove(C)=[];
           A11=A;
           A11(:,indsToRemove)=[];
           A11(indsToRemove,:)=[];
           v1=A(:,d);
           v1(indsToRemove)=[];
           
           x=-A11\v1;
           for u=1:length(C)
               deltaF(C(u))=x(u);
           end
           end
        end
        function [s,j]=maxStep(f,a,deltaF,deltaA,d)
            s=inf;
            j=-1;
            if(deltaA(d)>0)
                j=d;
                s=-a(d)/deltaA(d);
            end
            for q=1:length(C)
                if(deltaF(C(q))<0)
                    s1=-f(C(q))/deltaF(C(q));
                    if(s1<s-1e-14)
                        s=s1;
                        j=C(q);
                    end
                end
            end
            for q=1:length(NC)
                if(deltaA(NC(q))<0)
                    s1=-a(NC(q))/deltaA(NC(q));
                    if(s1<s-1e-14)
                        s=s1;
                        j=NC(q);
                    end
                end
            end
        end
    end
    function physicsStep(time)
        %updates positions of particles for a given time interval
        while (time>0)
            %obtaining velocities
            %--------------------------
            %velocities are calculated using mid point time stepping scheme
            u=getVelocities(1);
            velocities=reshape(u,[3,Nmoving]);
            initPos=positions;
            for m=1:Nmoving
                if(step_size>1e-6)
                     positions(m,:)=positions(m,:)+velocities(:,m)'*step_size;
                else
                    positions(m,:)=positions(m,:)+velocities(:,m)'*step_size/100;
                end
                positions(m,:)=positions(m,:)+periodicBoundary(positions(m,:));
            end
          
            u1=getVelocities(0);
            positions=initPos;
            if(step_size>1e-6) %derivative method is only valid for small step size
                velocities=reshape((u1+u)/2,[3,Nmoving]);
            else
                %for larger step size trapezoid rule is used
                velocities=reshape(u+50*(u1-u),[3,Nmoving]);
            end
            %------------------------------
            
            
            %now every particle is moved for time step_size until it 
            %collides with another particle
            %after the collision in can move in parallel direction to the
            %surface of particle to which it collided
            
            
            
            
            %continuous collision detection (the exact time of collision is
            %calculated for every possible collision)
            
            %-------------------------------------------------
            
            movingIndex=1;
            while (movingIndex<Nmoving+1)
                %loops over all moving particles
                i=movingIndices(movingIndex);
                  
                    vel=velocities(:,movingIndex);
                    move=(vel*step_size)'; %maximum move if there are no collisions
                    

                    deltaTime=step_size;
                    numberOfIter=0;
                    
                    while(deltaTime>0 && any(move) && ~deposed(i) && numberOfIter<10)

                            
                            %preventing particles to move in directions
                            %that would cause overlaps with particles that
                            %are already in contact
                            move=fixTouching(i,move')';
                            
                            %if move is too large reduces it according to
                            %periodic boundary conditions
                            move=periodicShiftMove(i,move);
                            
                            
                            %new velocity from this fixed move
                            vel=move'./deltaTime;
                        
                        if(any(move))
                            notMoving(i)=0;
                            if(positions(i,1)<xMin+2*particleRadius||positions(i,1)+move(1)<xMin+2*particleRadius) %possible collision with images through the wall at xMin
                                
                                shift(1)=xMax-xMin;
                            end
                            if(positions(i,1)>xMax-2*particleRadius||positions(i,1)+move(1)>xMax-2*particleRadius) %collision with the wall at xMax
                                
                                shift(1)=-xMax+xMin;
                            end
                            if(positions(i,2)<yMin+2*particleRadius||positions(i,2)+move(2)<yMin+2*particleRadius) %collision with the wall at yMin
                                
                                shift(2)=yMax-yMin;
                            end
                            if(positions(i,2)>yMax-2*particleRadius||positions(i,2)+move(2)>yMax-2*particleRadius) %collision with the wall at yMax
                                
                                shift(2)=-yMax+yMin;
                            end
                            
                            %stores particle that have been tested idices
                            %to prevent double checking
                            tested=zeros(16,1);
                            nTested=0;
                            
                            %gets cells through which particle moves
                            [cellsMove,numbC]=getCellsMoving(i,move);
                            
                            %detecting all possible collisions between the
                            %i th particle and all particles that are in its way 
                            for k=1:numbC
                                column=cellsMove(k,1);
                                row=cellsMove(k,2);
                                for itr=1:nInGrid(column,row)
                                    j=grid(column,row,itr);
                                    
                                    
                                    
                                    if(i~=j && abs(positions(i,3)-positions(j,3))<abs(move(3))+(2+1e-5)*particleRadius  && any(move) && ~any(tested(1:nTested)==j) &&~any(touchingInfo(i,2:touchingInfo(i,1)+1)-j==0))
                                        
                                        nTested=nTested+1;
                                        tested(nTested)=j;
                                        
                                        if(norm(positions(i,:)-positions(j,:))<norm(move)+particleRadius*2) % a collision might occur
                                            %more detailed calculations of
                                            %the time of collision
                                            delta_r=positions(j,:)-positions(i,:);
                                            discr=(delta_r*vel)^2-vel'*vel*(delta_r*delta_r'-(2*particleRadius)^2);
                                            if(discr>=0) %collision will occur
                                                
                                                dt=(delta_r*vel-sqrt(discr))/(vel'*vel);
                                                
                                                %now checks if the collision
                                                %occurs withing this step
                                                %and whether this is the first
                                                %collision
                                                
                                                if(dt>0 && dt<deltaTime&& dt<minT-1e-17)
                                                    
                                                    minT=dt;
                                                    coll=j;
                                                    fshift=[0,0,0];
                                                    nSimult=0;
                                                elseif(dt>0 && dt<deltaTime&& dt<minT+1e-17)
                                                    %collision with more
                                                    %than one sphere
                                                    nSimult=nSimult+1;
                                                    simultColl(nSimult)=j;
                                                    
                                                    
                                                end
                                                
                                            end
                                            
                                        end
                                        
                                        
                                        
                                        %if a particle is close to the
                                        %boundary it might collide wth
                                        %particles on the opposite side
                                        %of the box
                                        
                                        %the same checks as above are
                                        %repeated for particles in the
                                        %opposite side
                                        %--------------------------
                                        if(shift(1)~=0 && shift(2)~=0)
                                            delta_r=positions(j,:)-positions(i,:)-shift;
                                            disc=discriminant(delta_r,vel);
                                            if(disc>=0)
                                                dt=(delta_r*vel-sqrt(disc))/(vel'*vel);
                                                
                                                eif(dt>0&& dt<minT-1e-17)
                                                minT=dt;
                                                coll=j;
                                                fshift=shift;
                                                nSimult=0;
                                            elseif(dt>0&& dt<minT+1e-17)
                                                nSimult=nSimult+1;
                                                simultColl(nSimult)=j;
                                                
                                            end
                                        end
                                    end
                                    
                                    
                                    if(shift(1)~=0 )
                                        delta_r=positions(j,:)-positions(i,:)-[shift(1),0,0];
                                        disc=discriminant(delta_r,vel);
                                        if(disc>=0)
                                            dt=(delta_r*vel-sqrt(disc))/(vel'*vel);
                                            if(dt>0&& dt<minT-1e-17)
                                                minT=dt;
                                                coll=j;
                                                fshift=[shift(1),0,0];
                                                nSimult=0;
                                            elseif(dt>0&& dt<minT+1e-17)
                                                nSimult=nSimult+1;
                                                simultColl(nSimult)=j;
                                                
                                            end
                                        end
                                    end
                                    
                                    
                                    if(shift(2)~=0)
                                        delta_r=positions(j,:)-positions(i,:)-[0,shift(2),0];
                                        disc=discriminant(delta_r,vel);
                                        if(disc>=0)
                                            dt=(delta_r*vel-sqrt(disc))/(vel'*vel);
                                            if(dt>0&& dt<minT-1e-17)
                                                minT=dt;
                                                coll=j;
                                                fshift=[0,shift(2),0];
                                                nSimult=0;
                                            elseif(dt>0&& dt<minT+1e-17)
                                                nSimult=nSimult+1;
                                                simultColl(nSimult)=j;
                                                
                                            end
                                        end
                                        
                                    end
                                    %--------------------------
                                    %
                                end
                                
                            end
                            
                        
                        
                        
                        
                            
                            if(positions(i,3)+move(3)<particleRadius) %collision with the surface
                               
                                
                                dt=(particleRadius-positions(i,3))/vel(3);
                              
                                if(dt<minT && dt<deltaTime)
                                    minT=dt;
                                    coll=-1;
                                    deposed(i)=1;
                                    Nmoving=Nmoving-1;
                                    movingIndices(movingIndex)=[];
                                    movingIndex=movingIndex-1;
                                    
                                   
                                end
                            
                                
                            end
                            
                            
                            
                            if(coll~=0)
                                
                             
                                %updating the remaining time of the step and
                                %the position
                                deltaTime=deltaTime-minT;
                                positions(i,:)=positions(i,:)+vel'.*minT;
                                
                                if(coll>0)
                               %collsion with another particle
                                    
                                    %removes component of move in the
                                    %direction of line connecting centers
                                    %of colliding spheres
                                    move=(move-normr((positions(coll,:)-positions(i,:)-fshift)).*(normr(positions(coll,:)-positions(i,:)-fshift)*move'));
                                    
                                    
                                    addTouching(i,coll); %storing the indices of particles in contact
                                    if(~deposed(coll))
                                        addTouching(coll,i);
                                    end
                                    
                                    if(nSimult>0) %collision with more than 1 sphere
                                        for w=1:nSimult
                                            addTouching(i,simultColl(w));
                                            if(~deposed(simultColl(w)))
                                                addTouching(simultColl(w),i);
                                            end
                                        end
                                    end
                                    
                                    
                              
                              
                                end
                                
                            
                               
                            else
                                %no collision
                                positions(i,:)=positions(i,:)+move;
                                
                                deltaTime=0;

                            end
                            
                            %if particle is outside the box it is shifted
                            %to the opposite side
                            perShift=periodicBoundary(positions(i,:)); 
                            positions(i,:)=positions(i,:)+perShift;
                            if(drawOn)
                               figures(i)=drawSphere(i,[0 0 0]); %draws each frame
                            end
                 
                            
                            
                            %updating information of particles in contact
                            
                            keys=touchingInfo(i,2:touchingInfo(i,1)+1);


                            for q=1:length(keys)
                              
                                if(norm(positions(keys(q),:)-positions(i,:)+getShiftDeltaR(positions(keys(q),:)-positions(i,:)))>(2+1e-5)*particleRadius && coll~=keys(q))
                                    removeTouching(i,keys(q));
                                    if(~deposed(keys(q)))
                                        removeTouching(keys(q),i);
                                    end
                                end
                                
                            end
                            
                           
                            %-----------------space partition
                            %updating the grid
                            if(any(inCells(1,:,i)))
                                for c=1:4
                                   removeFromGrid(inCells(c,1,i),inCells(c,2,i),i);
                                end
                            end
                            
                            inCells(:,:,i)=getCells(positions(i,:));
                            for c=1:4
                               addToGrid(inCells(c,1,i),inCells(c,2,i),i);
                            end
                            %-----------------space partition
                            
                            %particle can stick to a deposed particle during
                            %a collision
                            if(coll>0 && deposed(coll) && stick())
                                deposed(i)=1;
                                Nmoving=Nmoving-1;
                                movingIndices(movingIndex)=[];
                                movingIndex=movingIndex-1;
                                break;
                            end
                        else
                            
                            notMoving(i)=notMoving(i)+1;
                            if(notMoving(i)>notMovingLimit)
                                deposed(i)=1;
                                Nmoving=Nmoving-1;
                                movingIndices(movingIndex)=[];
                                movingIndex=movingIndex-1;
                                break;
                            end
                            
                        end
                        
                    end
                    
                movingIndex=movingIndex+1;
            end
        
            time=time-step_size;
        end
        
    end

    function d=discriminant(delta_r,vel)
        d=(delta_r*vel)^2-vel'*vel*(delta_r*delta_r'-(2*particleRadius)^2);
    end
    function cells=get16Cells(i)
        %16 closest cells to a particle
        cells=zeros(16,2);
        cells4(1:4,:)=getCells(positions(i,:));
        minColumn=min(cells4(:,1));
        maxColumn=max(cells4(:,1));
        
        minRow=min(cells4(:,2));
        maxRow=max(cells4(:,2));
        cells(1:4,:)=cells4;
        if(minColumn<maxColumn-1)
          
            minColumn=maxColumn;
           
        end
        if(minRow<maxRow-1)
         
            minRow=maxRow;
           
        end
        m=minColumn-1;
        n=minRow-1;
        for i=1:12
            if(m<1)
                m=nColumns;
            elseif(m>nColumns)
                m=1;
            end
            if(n<1)
                n=nRows;
            elseif(n>nRows)
                n=1;
            end
            cells(4+i,:)=[m n];
             if(i<4)
                 n=n+1;
             elseif(i<7)
                 m=m+1;
             elseif(i<10)
                 n=n-1;
             else
                 m=m-1;
             end
        end
            
        
    end
    function cells=getCells(pos)
        %4 closest cells to a point at pos
        
        m=ceil((pos(1)-xMin)/deltaX);
        n=ceil((pos(2)-yMin)/deltaY);
    
        
        m1=round((pos(1)-xMin)./deltaX);
        n1=round((pos(2)-yMin)./deltaY);
       
        if(m1==m)
            m1=m+1;
        end
        if(n1==n)
            n1=n+1;
        end
        
        if(m1>nColumns)
            m1=1;
            
        end
        if(m1<1)
            m1=nColumns;
            
        end
        
        if(n1<1)
            n1=nRows;
            
        end
        if(n1>nRows)
            n1=1;
            
        end
        
        cells=[[m,n];[m,n1];[m1,n];[m1,n1]];
        
    end
    function [cells,numbOfCells]=getCellsMoving(i,move)
        %finds cells through which particle i moves by doing rasterisation
        %of move
        n=0;
        only1=0;
        numbOfCells=0;
        chkDist=[mod(positions(i,1)-xMin+deltaX/2,deltaX),mod(positions(i,2)-yMin+deltaY/2,deltaY)];
        if(chkDist(1)+move(1)<0 || chkDist(1)+move(1)>deltaX || chkDist(2)+move(2)<0 || chkDist(2)+move(2)>deltaY)
            grad=abs(move(2)/move(1));
            
            if(grad<1)
                bCoord=1;
                delta=deltaX;
                
                step=[delta*sign(move(1)), grad*delta*sign(move(2)), 0];
            else
                bCoord=2;
                delta=deltaY;
                
                step=[delta/grad*sign(move(1)), delta*sign(move(2)), 0];
            end
            if(delta>abs(move(bCoord)))
                
                delta=move(bCoord);
                step=[move(1),move(2),0];
            end
        else
            bCoord=1;
            delta=deltaX;
            only1=1;
            step=[0, 0, 0];
        end
        cells=zeros(1+floor(move(bCoord)/delta),2); %TODO FIX wrong size
        pos=-delta;
        while(abs(pos)<abs(move(bCoord)) || n==0 ||only1)
            position=positions(i,:)+n*step;
            currCells=getCells(position+periodicBoundary(position));
            numbOfAdded=0;
            for k=1:4
                if(~any(find(cells(:,1)==currCells(k,1) & cells(:,2)==currCells(k,2),1)))
                    numbOfCells=numbOfCells+1;
                    cells(numbOfCells,:)=currCells(k,:);
                    numbOfAdded=numbOfAdded+1;
                end
            end
            if(numbOfAdded>2 && n>0)
                
                %correction of rasterization of the center's trajectory
                
                der=step(1)/step(2);
                deltaXY=[ternary(step(1)>0,deltaX,0),ternary(step(2)>0,deltaY,0)];
                shiftR=0.5*[der*(deltaXY(2)-chkDist(2))+deltaXY(1)-chkDist(1),(deltaXY(1)-chkDist(1))/der+deltaXY(2)-chkDist(2),0];
                position=position-step+shiftR;
                currCells=getCells(position+periodicBoundary(position));
                for k=1:4
                    if(~any(find(cells(:,1)==currCells(k,1) & cells(:,2)==currCells(k,2),1)))
                        numbOfCells=numbOfCells+1;
                        cells(numbOfCells,:)=currCells(k,:);
                        numbOfAdded=numbOfAdded+1;
                    end
                end
            end
            n=n+1;
            pos=pos+delta;
            if(only1)
                only1=0;
            end
        end
        
        
    end

    function s=stick()
        s=rand()<stickProbability;
    end
    function x=ternary(condition,valueY,valueN)
        if(condition)
            x=valueY;
        else
            x=valueN;
        end
    end

end
