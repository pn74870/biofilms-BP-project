function [ cov] = SimulationPolydisp( step_size,time,rate)
%Simulation of spherical particle deposition

tic; %start timer





persistent particleRadius;
persistent xMax; %box corners
persistent xMin;
persistent yMin;
persistent yMax;
persistent N; %number of created particles
persistent movingIndices;
%coordinates of box corners
xMax=0.5e-5;
xMin=-xMax;
yMax=0.5e-5;
yMin=-yMax;

zMax=.5e-5; %height of the box

notMovingLimit=3;%if a particle does not move this number of times in a row it is deposed  

%rng(1);
%data of particles:
positions=zeros(2*time*rate,3);
positionsCenter=zeros(time*rate,3);
dirOfBrownF1=zeros(3,time*rate,time*rate);
dirOfBrownF2=zeros(3,time*rate);


omegas=zeros(3,time*rate);
deposed=zeros(1,time*rate*2); %contains information whether a particle was deposed onto surface
touchingInfo=cell(1,time*rate*2); %contains information about particles that are in contact
figures=gobjects(time*rate*2); %contains graphical objects of spheres
Nmoving=0;
velCenters=zeros(time*rate,3);
quaternions=zeros(4,time*rate);

%---------------
%physical parameters
particleRadius=1e-6;
T=300; %K
kb=1.38e-23;
d=2*particleRadius; %diameter
eta=8.9e-4;%viscocity
mass=1e-14;
dens=mass/(4/3*pi*particleRadius^3);
gforce=mass*(1-1000/dens)*[0,0,-9.8];%effective force due to gravity and buoyancy
stickProbability=0.1; %probability that two particles stick together after they collide
epsilon0=8.85e-12;
epsilonR=78;
charge=1.6e-19;
chargeConcentration=1e-4*6.02e29;%number density of ions in liquid
DebLength=sqrt(charge^2*chargeConcentration/epsilon0/epsilonR/kb/T);
potential=1e-3;%surface potential of the particles

%constant of proportionality for electric double layer force:
elForceConst=2*pi*epsilon0*epsilonR*DebLength*particleRadius*potential^2;

numbDensPart=3e26;
lje=1e3/6.02e23;
ljs=particleRadius/5; %radius of Lennard Jones particles of which colloid particles are made of
Acc=1e-20;%4*pi^2*lje*numbDensPart^2*ljs^6;

nSpheresInPart=2;
coframePos=[particleRadius,-particleRadius;0,0;0,0];
momOfInert0=particleRadius^2/5*[4,0,0;0,14,0;0,0,14];
%---------------

pairDir=containers.Map('KeyType','int32','ValueType','any');



virtFigures=gobjects(1,time*rate);


for o=1:time*rate*2
    touchingInfo{o}=containers.Map('KeyType','int32','ValueType','any');
end

notMoving=zeros(1,time*rate);

nColumns=floor((xMax-xMin)/2/particleRadius);
nRows=floor((yMax-yMin)/2/particleRadius);
nLayers=floor(zMax/2/particleRadius);
grid=zeros(nColumns,nRows,4*nLayers);
nInGrid=zeros(nColumns,nRows);
%grid=zeros(floor((xMax-xMin)/particleRadius),floor((yMax-yMin)/particleRadius),8);

inCells=zeros(4,2,time*rate*2); %every particle is in up to 4 cells of the grid
deltaX=(xMax-xMin)/nColumns;
deltaY=(yMax-yMin)/nRows;
deltaZ=zMax/nLayers;
tested=java.util.HashSet;
added=java.util.HashSet;

%velocity=[0;0;-50];


N=0; %number of created particles

calcCov=0; %1 to calculate surface coverage as a function of time

simultColl1=zeros(1,3);
simultShift1=zeros(3,3);

simultColl2=zeros(1,3);
simultShift2=zeros(3,3);


nOfZSteps=100;
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

cov=zeros(1,time);

interval=1/rate;
nIntervals=time*rate;

meanMove=zeros(1,time*rate);
sumOfMoves=0;
nMoves=0;
drawOn=1;
for t=1:nIntervals
    
    
    
    %sp(t,:)=positions(t,:);
    %positions(t,:)=sp(t,:);
    N=N+1;
    
    Nmoving=Nmoving+1;
    movingIndices(Nmoving)=t;
    %physics simulation
    wrongPos=1;
    while(wrongPos)
        %wrongPos=0;
        positionsCenter(t,:)=[xMin+rand()*(xMax-xMin);yMin+rand()*(yMax-yMin);zMax];
        quaternions(:,t)=randomQuaternion();
        [pos1,pos2]=getPositions(t);
        pos1=pos1+periodicBoundary(pos1);
        pos2=pos2+periodicBoundary(pos2);
        [wrongPos1,cells1]=checkColl(pos1);
        [wrongPos2,cells2]=checkColl(pos2);
        wrongPos= wrongPos1|| wrongPos2;
    end
    positions(2*t-1,:)=pos1;
    positions(2*t,:)=pos2;
    inCells(:,:,2*N-1)=cells1;
    inCells(:,:,2*N)=cells2;
    for cl=1:4
        addToGrid(inCells(cl,1,2*N),inCells(cl,2,2*N),2*N);
        addToGrid(inCells(cl,1,2*N-1),inCells(cl,2,2*N-1),2*N-1);
    end
    
%     nMoves=0;
%     sumOfMoves=0;
     physicsStep(interval);
%     meanMove(t)=sumOfMoves/nMoves;
    
    
    
    %drawSphere(particle.position);
    if (mod(t,rate)==0 && calcCov)
        cov(t/rate)=coverage1();
    end
    fprintf('%.1f percent completed\n',t/nIntervals*100);
end

drawAll();
hold off;
if (calcCov)
    figure();
    plot(cov);
end
figure();
plot(meanMove);
%volume fraction as a function of z
% zTop=max(positions(:,3));
% zStep=zTop/nOfZSteps;
% v=zeros(1,nOfZSteps);
% zVolF=linspace(0,zTop,nOfZSteps);
% for s=1:nOfZSteps
%     v(s)=volumeFraction(zStep*(s-1)+particleRadius,zStep);
% end
% 
% figure();
% plot(zVolF,v);

fprintf('time %.3f s\n',toc);

 function [wrongPos,cells]=checkColl(pos)
        wrongPos=0;
        cells=getCells(pos);
         for cel=1:4
             col=cells(cel,1);
             ro=cells(cel,2);
             
             tested.clear();
              for it=1:nInGrid(col,ro)
                 p=grid(col,ro,it);
                 deltaR=positions(p,:)-pos;
                 deltaR=deltaR+getShiftColl(deltaR);
                 if(abs(deltaR(3))<2*particleRadius && norm(deltaR)<2*particleRadius && tested.add(p))
                     wrongPos=1;
                     break;
                 end
             end
             if(wrongPos)
                 break;
             end
         end
    end
        function addToGrid(column,row,i)
        nInGrid(column,row)=nInGrid(column,row)+1;
        grid(column,row,nInGrid(column,row))=i;
    end
    function removeFromGrid(column,row,i)
        sizeOfGrid=size(grid);
        indexArr=reshape(grid(column,row,:),[1,sizeOfGrid(3)]);
        indexArr=[indexArr(indexArr~=i),0];
        grid(column,row,:)=reshape(indexArr,[1,1,sizeOfGrid(3)]);
        nInGrid(column,row)=nInGrid(column,row)-1;
    end
    function drawAll()
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
    function drawFrame()
        for i=1:N
            [figures(2*i-1),figures(2*i)]=drawSpheres(i,[0 0 0]);
        end
    end
function fig=drawSphere(p,c)
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


function fig=drawSphereTest(pos)
        [x,y,z]=sphere;
       
        fig=surf(x*particleRadius+pos(1),y*particleRadius+pos(2),z*particleRadius+pos(3));
       
        
        drawnow;
        
        
    end
    function q=multiplyQuaternions(q1, q2)
        v1=q1(2:4);
        v2=q2(2:4);
       q=[q1(1)*q2(1)-v1'*v2; q1(1)*v2+q2(1)*v1+cross(v1,v2)];
       q=q/norm(q');
        
    end
    function q=randomQuaternion()
        fi=2*pi*rand();
        theta=acos(2*rand()-1);
        q1=[cos(pi/4-theta/2);0;sin(-pi/4+theta/2);0];
        q2=[cos(fi/2);0;0;sin(fi/2)];
        q=multiplyQuaternions(q1,q2);
        
    end
    function R=quaternionToMatrix(q)
        q=q/norm(q');
        R=[1-2*q(3)^2-2*q(4)^2,2*q(2)*q(3)-2*q(1)*q(4),2*q(2)*q(4)+2*q(1)*q(3);
            2*q(2)*q(3)+2*q(1)*q(4), 1-2*q(2)^2-2*q(4)^2, 2*q(3)*q(4)-2*q(1)*q(2);
            2*q(2)*q(4)-2*q(1)*q(3), 2*q(3)*q(4)+2*q(1)*q(2), 1-2*q(2)^2-2*q(3)^2];
        
    end
    function [pos1,pos2]=getPositions(i)
        R=quaternionToMatrix(quaternions(:,i));
        pos1=(R*coframePos(:,1)+positionsCenter(i,:)')';
        pos2=(R*coframePos(:,2)+positionsCenter(i,:)')';
        pos1=pos1+periodicBoundary(pos1);
        pos2=pos2+periodicBoundary(pos2);
    end
    function shift=periodicBoundary(position)
        sideX=xMax-xMin;
        sideY=yMax-yMin;
        shift=[0,0,0];
        if(xMin>position(1))
            shift(1)=sideX;
        elseif (xMax<position(1))
            shift(1)=-sideX;
        end
        
        if(yMin>position(2))
            shift(2)=sideY;
        elseif (yMax<position(2))
            shift(2)=-sideY;
        end
        
    end

    function c=coverage1()
        totalArea=(xMax-xMin)*(yMax-yMin);
        coveredArea=0;
        % g=zeros(numbOfSteps,numbOfSteps);
        for i=1:N
            if (positions(i,3)<=particleRadius)
                coveredArea=coveredArea+pi*particleRadius^2;
            end
        end
        c=coveredArea/totalArea;
    end


function [ dx,dy,dz ] = BrownianMotion( stepSize,T,eta,d,force )
%generates random displacement
kb=1.38e-23;
D=kb*T/3/pi/eta/d;
k=sqrt(2*D*stepSize);
dx=k*randn()+D/kb/T*force(1)*stepSize;
dy=k*randn()+D/kb/T*force(2)*stepSize;
dz=k*randn()+D/kb/T*force(3)*stepSize;


end

    function c=coverage(numbOfSteps)
        totalArea=(xMax-xMin)*(yMax-yMin);
        stepX=(xMax-xMin)/numbOfSteps;
        stepY=(yMax-yMin)/numbOfSteps;
        coveredArea=totalArea;
        % g=zeros(numbOfSteps,numbOfSteps);
        parfor i=1:numbOfSteps
            
            for j=1:numbOfSteps
                free=1;
                
                for par=1:N
                    w=0;
                    point=[xMin+i*stepX,yMin+j*stepY];
                    if (norm([positions(par,1)-point(1);positions(par,2)-point(2)])<particleRadius)
                        w=1;
                    end
                    if(w && deposed(par))
                        free=0;
                        break;
                        %  g(mod(i,numbOfSteps)+1,floor(i/numbOfSteps)+1)=1;
                    end
                end
                if(free)
                    coveredArea=coveredArea-stepX*stepY;
                end
            end
        end
        c=coveredArea/totalArea;
        %         figure();
        %         hold on;
        %         for x=1:numbOfSteps
        %             for y=1:numbOfSteps
        %                 if (g(x,y)==1)
        %                     fill3([(x-1)*stepX,x*stepX,x*stepX,(x-1)*stepX],[(y-1)*stepY,(y-1)*stepY,y*stepY,y*stepY],[0,0,0,0],'b');
        %                 end
        %             end
        %         end
        %         hold off;
    end

    function [fig1,fig2]=drawSpheres(p,c)
        if(~isempty(figures))
            try
                delete(figures(2*p-1));
                delete(figures(2*p));
            catch
            end
        end
        [x,y,z]=sphere;
        [p1,p2]=getPositions(p);
        if (any(c))
            fig1=surf(x*particleRadius+p1(1),y*particleRadius+p1(2),z*particleRadius+p2(3));
            fig2=surf(x*particleRadius+p2(1),y*particleRadius+p2(2),z*particleRadius+p2(3));
            set(fig1, 'FaceColor', c);
            set(fig2, 'FaceColor', c);
        else
            fig1=surf(x*particleRadius+p1(1),y*particleRadius+p1(2),z*particleRadius+p1(3));
            fig2=surf(x*particleRadius+p2(1),y*particleRadius+p2(2),z*particleRadius+p2(3));
            
        end
        
        drawnow;
        
        
    end
    function fig=drawVirtualSphere(p,shift)
        if(~isempty(virtFigures))
            try
                delete(virtFigures(p));
            catch
            end
        end
        [x,y,z]=sphere;
        
        
        fig=surf(x*particleRadius+positions(p,1)+shift(1),y*particleRadius+positions(p,2)+shift(2),z*particleRadius+positions(p,3)+shift(3));
        set(fig, 'FaceColor', [1,1,1]);
        drawnow;
        
    end
function shift=getShiftColl(deltaR)
    shift=[0 0 0];
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

function shift=getShiftDeltaR(deltaR)
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
 function shift=getShift(pos,dist,move)
         shift=[0,0,0];
         if(pos(1)<xMin+dist||pos(1)+move(1)<xMin+dist) %possible collision with images through the wall at xMin
             
             shift(1)=xMax-xMin;
         end
         if(pos(1)>xMax-dist||pos(1)+move(1)>xMax-dist) %collision with the wall at xMax
             
             shift(1)=-xMax+xMin;
         end
         if(pos(2)<yMin+dist||pos(2)+move(2)<yMin+dist) %collision with the wall at yMin
             
             shift(2)=yMax-yMin;
         end
         if(pos(2)>yMax-dist||pos(2)+move(2)>yMax-dist) %collision with the wall at yMax
             
             shift(2)=-yMax+yMin;
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
        fDisp1=-32/3*Acc*particleRadius^6/r^3/(r^2-4*particleRadius^2)^2; %r^-6 term of LJ potential
        fDisp2=Acc/37800*ljs^6*((56*particleRadius*r^2+10*particleRadius^3-6*r^3)/(r-2*particleRadius)^8/r^2-(r^3+26*particleRadius*r^2+108*particleRadius^3+432*particleRadius^2*r)/(r+2*particleRadius)^8/r^2+12*(r^2-20*particleRadius^2)/r^9);
        f=fDisp1+fDisp2;
    end
    function w=smoothingFunc(x,b)
        w=1-(1-x).^b;
    end
    function r=random1Variance
        r=randn(3,1);
    end


    function [row,force]=getRMatrixRow(moveInd,keepDir)
       i=movingIndices(moveInd);
        f=[0 0 0];
        row=zeros(3,3*Nmoving);
        
        if(any(dirOfBrownF2(:,i)))
            fBrown=dirOfBrownF2(:,i);
            dirOfBrownF2(:,i)=[0;0;0];
        else
            fBrown=random1Variance();
            dirOfBrownF2(:,i)=fBrown;
        end
        
        
        nNear=0;
        ASumMatrix=zeros(3);
        cells=get16Cells(i);
        added.clear();
        for c=1:16
            
          for itr=1:nInGrid(cells(c,1),cells(c,2))
                j=grid(cells(c,1),cells(c,2),itr);
                if(i~=j && abs(positions(i,3)-positions(j,3))<=4*particleRadius&& added.add(j))
                   
                    deltaR=positions(j,:)-positions(i,:);
                    shift=getShiftColl(deltaR);
                    deltaR=deltaR+shift;
                    r=norm(deltaR);
                    if(r<=4*particleRadius)
                        
                        
                        
                        
                        
                        dir=deltaR/r;
                        if(r<3*particleRadius)
                            del=r-2*particleRadius;
                             del1=particleRadius/1000;
                            if(del>del1)
                                
                                Xpw=particleRadius/4/del+9/40*log(particleRadius/del);
                                Ypw=log(particleRadius/del)/6;
                                smooth=(smoothingFunc((3-r/particleRadius),4))^2;
                                Xpw=Xpw*smooth;
                                Ypw=Ypw*smooth;
                            else
                                
                                Xpw=particleRadius/4/del1+9/40*log(particleRadius/del1);
                                Ypw=log(particleRadius/del1)/6;

                                
                            end
                            D=dir'*dir;
                            I=eye(3);
                            Apw=Xpw*D+Ypw*(I-D);
                            if(~deposed(j))
                                k=find(movingIndices==j);
                                row(:,k*3-2)=Apw(:,1);
                                row(:,k*3-1)=Apw(:,2);
                                row(:,k*3)=Apw(:,3);
                            end
                            ASumMatrix=ASumMatrix+Apw;
                          
                            if(any(dirOfBrownF1(:,i,j)))
                                randVec=dirOfBrownF1(:,i,j);
                                if(~keepDir)
                                    dirOfBrownF1(:,i,j)=[0;0;0];
                                end
                            else
                                randVec=random1Variance();
                                dirOfBrownF1(:,j,i)=randVec;
                                dirOfBrownF1(:,i,j)=randVec;
                            end
                            
                            
                            fBrown=fBrown+(sqrt(Xpw)*D+sqrt(Ypw)*(I-D))*randVec;
                        end
                        fel=elForceConst/(1+exp(DebLength*(r-2*particleRadius)));
                        
                        if(r>2*particleRadius+30^(-1/6)*ljs) %recommended cut-off distance (potential is divergent at r=2*radius)
                            fDisp=dispElForce(r);
                            
                        else
                            fDisp=dispElForce(2*particleRadius+30^(-1/6)*ljs);
                        end
                        f=f-(fel+fDisp)*dir;
                    end
                end
            end
        end
        
        volFrac=N*4/3*pi*particleRadius^3/(xMax-xMin)/(yMax-yMin)/zMax;
        RCorrection=1+2.725*volFrac-6.583*volFrac^2;    
        row(:,3*moveInd-2)=ASumMatrix(:,1)+[1;0;0]*RCorrection;
        row(:,3*moveInd-1)=ASumMatrix(:,2)+[0;1;0]*RCorrection;
        row(:,3*moveInd)=ASumMatrix(:,3)+[0;0;1]*RCorrection;
        row=row*6*pi*eta*particleRadius;
    
        force=f'+fBrown*sqrt(12*pi*eta*particleRadius*kb*T/step_size)+gforce';
    end

    function u=getVelocities(keepDir)
        
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

   
 function newVel=fixTouching(partInd,vel)
%         function b=isCollSoved()
%            for c=1:n
%                if(
%         end
        newVel=vel;
        n=0;
       for b=1:length(partInd)
        n=n+length(touchingInfo{partInd(b)}.keys());
           
       end
       impulses=zeros(n,1);
       contactSpeed=zeros(n,1);
       contactNorms=zeros(n,3);
       collNotSolved=0;
       A=zeros(n);
       indices=[1:n];
       C=[];
       NC=[];
       d0=0;
        m=0;
       for i=1:length(partInd)
           keys=touchingInfo{partInd(i)}.keys();
           for c=1:length(touchingInfo{partInd(i)}.keys())
               
               m=m+1;
               dir=normr(positions(keys{c},:)-positions(partInd(i),:)-touchingInfo{partInd(i)}(keys{c}));
               contactNorms(m,:)=-dir;
               contactSpeed(m)=-dir*vel;
               if(contactSpeed(m)<0)
                   collNotSolved=1;
                   if(d0==0)
                       d0=m;
                   end
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
               
               if(contactSpeed(c)<0 && ~any(C(C==c)))
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
                    if(s1<s)
                        s=s1;
                        j=C(q);
                    end
                end
            end
            for q=1:length(NC)
                if(deltaA(NC(q))<0)
                    s1=-a(NC(q))/deltaA(NC(q));
                    if(s1<s)
                        s=s1;
                        j=NC(q);
                    end
                end
            end
        end
    end


    function physicsStep(time)
        %    initPos=positions;
        while (time>0)
%             u=getVelocities(1);
%             velocities=reshape(u,[3,Nmoving]);
%             initPos=positions;
%             for m=1:Nmoving
%                 positions(m,:)=positions(m,:)+velocities(:,m)'*step_size/100;
%                 positions(m,:)=positions(m,:)+periodicBoundary(positions(m,:));
%             end
%             %positions=positions+velocities'*step_size/100;
%             u1=getVelocities(0);
%             positions=initPos;
            
%             velocities=reshape(u+50*(u1-u),[3,Nmoving]);
            movingIndex=1;
            while (movingIndex<Nmoving+1)
                i=movingIndices(movingIndex);
                    %mid point algorithm
                    %----------------------
%                     u0=getVel(i,1);
%                     move=(u0*step_size)';
%                     initPos=positions(i,:);
%                     positions(i,:)=positions(i,:)+move;
%                     positions(i,:)=positions(i,:)+periodicBoundary(positions(i,:));
%                     u=getVel(i,0);
%                     positions(i,:)=initPos;
%                     velocities(:,i)=(u0+u)*.5;
         
%vel=velocities(:,movingIndex);
                    RotMatrixBodyWorld=quaternionToMatrix(quaternions(:,i));
                    rCenter=zeros(nSpheresInPart,3);
                    IMatrix=RotMatrixBodyWorld*momOfInert0/(RotMatrixBodyWorld);
                    sumV=[0;0;0];
                    L=[0;0;0];
                    for s=1:nSpheresInPart
                        [dx,dy,dz]=BrownianMotion(step_size,T,eta,d,gforce);%+elecForce(i));
                        move1=[dx,dy,dz];
                   % move=(vel*step_size)';
                       vel1=move1'/step_size;
                       sumV=sumV+vel1;
                       L=L+cross(positions(nSpheresInPart*i-2+s,:)'-positionsCenter(i,:)',vel1);
                       rCenter(s,:)=positions(nSpheresInPart*i-2+s,:)-positionsCenter(i,:);
                        rCenter(s,:)=rCenter(s,:)+getShiftDeltaR(rCenter);
                    end
                    
                
                    
                    
                    
                    omegas(:,i)=IMatrix\L;
                %    omegas(:,i)=[0;0;10];
                    vel=sumV/nSpheresInPart;
                    
                    move=vel'*step_size;
                    
                    
                    
                    
                    shift=getShift(positionsCenter(i,:),3*particleRadius,[0,0]);
                    
                    
                    nSimult2=0;
                    nSimult1=0;
                    
                    deltaTime=step_size;
                    numberOfIter=0;
                    coll=0;
                    
                    %TOUCHING ROTATIONS
                    r1Center=positions(2*i-1,:)-positionsCenter(i,:);
                    r1Center=r1Center+getShiftDeltaR(r1Center);
                    
                    r2Center=positions(2*i,:)-positionsCenter(i,:);
                    r2Center=r2Center+getShiftDeltaR(r2Center);
                    
                    instantVel=cross(omegas(:,i),r1Center');
                    wrongRot=0;
                    keys=touchingInfo{i*2-1}.keys();
                     for c=1:length(keys)
                         dir=positions(keys{c},:)-positions(2*i-1,:)-touchingInfo{2*i-1}(keys{c});
                         
                         if(dir*instantVel>1e-16)
                             omegas(:,i)=[0;0;0];
                             wrongRot=1;
                             break;
                         end
                     end
                     
                     if(~wrongRot)
                    keys=touchingInfo{i*2}.keys();
                     for c=1:length(keys)
                         dir=positions(keys{c},:)-positions(i*2,:)-touchingInfo{i*2}(keys{c});
                         
                         if(-dir*instantVel>1e-16)
                             omegas(:,i)=[0;0;0];
                             wrongRot=1;
                             break;
                         end
                     end
                     end
                        %------------------------------------------              
                        if(~wrongRot)
                            
                            minT=step_size;
                            omegaUnitVec=omegas(:,i)/norm(omegas(:,i)');
                            quaternion=[cos(norm(omegas(:,i))'*step_size/2);sin(norm(omegas(:,i)')*step_size/2)*omegaUnitVec]*step_size;
                            quaternion=quaternion/norm(quaternion');
                            moveRot=quaternionToMatrix(multiplyQuaternions(quaternion,quaternions(:,i)))*coframePos(:,1)-positions(2*i-1,:)'+positionsCenter(i,:)';
                            omegXY=omegaUnitVec-omegaUnitVec(3)*[0;0;1];
                            if(any(omegXY))
                                sinFi=normr(omegXY')*[0;1;0];
                                cosFi=normr(omegXY')*[1;0;0];
                            else
                                sinFi=0;
                                cosFi=1;
                            end
                            cosTheta=omegaUnitVec(3);
                            sinTheta=norm(cross(omegaUnitVec,[0;0;1]));
                            
                            RotMatrixWorldOmega=[cosTheta,0,-sinTheta;0,1,0;sinTheta,0,cosTheta]*[cosFi,sinFi,0;-sinFi,cosFi,0;0,0,1];
                            
%                             v=cross(omegaUnitVec,[0;0;1]);
%                             s=norm(v');
%                             c=omegaUnitVec'*[0;0;1];
%                             skewMatrix=[0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
%                             RotMatrixWorldOmega=eye(3)+skewMatrix+skewMatrix^2*(1-c)/s^2;
                            
                            
                            
                            tested.clear();
                            cells=get16Cells(positions(i,:));
                            
                            posInOmegaFrame1=RotMatrixWorldOmega*r1Center';
                            posInOmegaFrame2=RotMatrixWorldOmega*r2Center';
                            for ce=1:16
                                column=cells(ce,1);
                                row=cells(ce,2);
                               for itr=1:nInGrid(column,row)
                                    j=grid(column,row,itr);
                                    
                                    
                                    
                                    if(2*i-1~=j && 2*i~=j && abs(positionsCenter(i,3)-positions(j,3))<abs(moveRot(3))+abs(positionsCenter(i,3)-positions(2*i-1,3))+(2+1e-5)*particleRadius  && any(moveRot) && tested.add(j))
                                        
                                        if(~touchingInfo{2*i-1}.isKey(j))
                                            [dt,tshift]=getTimeToCollRot(i,j,posInOmegaFrame1,RotMatrixWorldOmega,shift);
                                            if(dt>0)
                                                if( dt<minT-1e-17)
                                                    minT=dt;
                                                    coll=j;
                                                    
                                                    partColl=2*i-1;
                                                    nSimult1=0;
                                                elseif(dt<minT+1e-17)
                                                    nSimult1=nSimult1+1;
                                                    simultColl1(nSimult1)=j;
                                                    
                                                end
                                            end
                                        end
                                        if(~touchingInfo{2*i}.isKey(j))
                                            [dt,tshift]=getTimeToCollRot(i,j,posInOmegaFrame2,RotMatrixWorldOmega,shift);
                                            if(dt>0)
                                                if( dt<minT-1e-17)
                                                    minT=dt;
                                                    coll=j;
                                                    
                                                    partColl=2*i;
                                                    nSimult2=0;
                                                elseif(dt<minT+1e-17)
                                                    nSimult2=nSimult2+1;
                                                    simultColl2(nSimult2)=j;
                                                    
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                                
                            
                            nVector=RotMatrixWorldOmega*[0;0;1];
                           
                            dFi=getAngleOfCollGround(posInOmegaFrame2,positionsCenter(i,3),nVector);
                            if(dFi>0)
                              
                                dt=dFi/norm(omegas(:,i)');
                                if(dt<minT)
                                    minT=dt;
                                    coll=-1;
                                    partColl=2*i;
                                end
                                
                            end
                           
                            
                            dFi=getAngleOfCollGround(posInOmegaFrame1,positionsCenter(i,3),nVector);
                            if(dFi>0)
                              
                                dt=dFi/norm(omegas(:,i)');
                                if(dt<minT)
                                    minT=dt;
                                    coll=-1;
                                    partColl=2*i-1;
                                end
                                
                            end
                            
                            if(minT>step_size)
                                coll=0;
                            end
                           
                            if(coll~=0)
                                angle=norm(omegas(:,i)')/2*minT;
                                quaternions(:,i)=multiplyQuaternions([cos(angle);sin(angle)*omegas(:,i)/norm(omegas(:,i)')],quaternions(:,i));
                                %quaternions(:,i)=quaternions(:,i)/norm(quaternions(:,i));
                                [positions(2*i-1,:),positions(2*i,:)]=getPositions(i);
                                delP=(quaternionToMatrix([cos(angle);sin(angle)*omegas(:,i)/norm(omegas(:,i)')])*r2Center')'+positionsCenter(i,:)+periodicBoundary( (quaternionToMatrix([cos(angle);sin(angle)*omegas(:,i)/norm(omegas(:,i)')])*r2Center')'+positionsCenter(i,:)) -positions(2*i,:);
                                
                                if(coll==-1)
                                    %deposed(partColl)=1; %TODO
                                    deposed(2*i)=1;
                                    deposed(2*i-1)=1;
                                    Nmoving=Nmoving-1;
                                    movingIndices(movingIndex)=[];
                                    movingIndex=movingIndex-1;
                                else
                                    fshift=getShiftDeltaR(positions(partColl,:)-positions(coll,:));
                                    touchingInfo{partColl}(coll)=fshift;
                                    if(~deposed(coll))
                                        touchingInfo{coll}(partColl)=-fshift;
                                    end
                                    if(nSimult1>0)
                                        for w=1:nSimult1
                                            fshift=getShiftDeltaR(positions(2*i-1,:)-positions(simultColl1(w),:));
                                            touchingInfo{2*i-1}(simultColl1(w))= fshift;
                                            if(~deposed(simultColl1(w)))
                                                touchingInfo{simultColl1(w)}(2*i-1)=-fshift;
                                            end
                                        end
                                    end
                                    if(nSimult2>0)
                                        for w=1:nSimult2
                                            fshift=getShiftDeltaR(positions(2*i,:)-positions(simultColl2(w),:));
                                            touchingInfo{2*i}(simultColl2(w))= fshift;
                                            if(~deposed(simultColl2(w)))
                                                touchingInfo{simultColl2(w)}(2*i)=-fshift;
                                            end
                                        end
                                    end
                                end
                            else
                                angle=norm(omegas(:,i)')/2*step_size;
                                quaternions(:,i)=multiplyQuaternions([cos(angle);sin(angle)*omegas(:,i)/norm(omegas(:,i)')],quaternions(:,i));
                               % quaternions(:,i)=quaternions(:,i)/norm(quaternions(:,i));
                                [positions(2*i-1,:),positions(2*i,:)]=getPositions(i);
                            end
                            onPositionChanged(2*i-1);
                            onPositionChanged(2*i);
                            if(drawOn)
                            figures(2*i-1)=drawSphere(2*i-1,[0 0 0]);   
                            figures(2*i)=drawSphere(2*i,[0 0 0]); 
                            end
                            for j=1:2*N
                                if(2*i~=j &&2*i-1~=j)
                                    p=2*i-1;
                                    if( ~touchingInfo{p}.isKey(j) && (norm(positions(p,:)-positions(j,:))-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+shift)-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+[shift(1),0,0])-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+[0,shift(2),0])-2*particleRadius<0))
                                        fprintf('critical error + %s\n',norm(positions(i,:)-positions(j,:))-2*particleRadius);
                                    end
                                    p=2*i;
                                    if( ~touchingInfo{p}.isKey(j) && (norm(positions(p,:)-positions(j,:))-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+shift)-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+[shift(1),0,0])-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+[0,shift(2),0])-2*particleRadius<0))
                                        fprintf('critical error + %s\n',norm(positions(i,:)-positions(j,:))-2*particleRadius);
                                    end
                                end
                            end
                        end       
        %------------------------------------------
                    
                    while(deltaTime>0 && any(move) && ~deposed(2*i) && numberOfIter<10)
                        
                        numberOfIter=numberOfIter+1;
                       shift1=getShift(positions(2*i-1,:),2*particleRadius,move);
                       shift2=getShift(positions(2*i,:),2*particleRadius,move);
                 
                        
                        %------------touching
                            oldVel=vel;
                            vel=fixTouching([2*i-1,2*i],vel);
%                             
%                             wrongMove=0;
                            keys=touchingInfo{2*i-1}.keys();
                            for c=1:length(keys)
                                dir1=positions(keys{c},:)-positions(2*i-1,:)-touchingInfo{2*i-1}(keys{c});
                                if(dir1*vel>1e-14)
                                    wrongMove=1;
                                    break;
                                end
                            end
                                keys=touchingInfo{2*i}.keys();
                            for c=1:length(keys)
                                dir1=positions(keys{c},:)-positions(2*i,:)-touchingInfo{2*i}(keys{c});
                                if(dir1*vel>1e-14)
                                    wrongMove=1;
                                    break;
                                end
                            end
%                             if(wrongMove)
%                                 move=[0,0,0];
%                             else
%                                 move=move2;
%                             end
%                             
                            move=vel'*deltaTime;
                            %------------------------------
                        
                        
                                                   
                                    if(any(move))
                                        
                                        minT=step_size;
                                        coll=0; %is there a collision
                                        fshift=[0,0,0];
                                        nSimult1=0;
                                        nSimult2=0;
                                  tested.clear();
                                  
                            [cellsMove,numbC]=getCellsMoving(2*i-1,move);
                            
                            for k=1:numbC
                                column=cellsMove(k,1);
                                row=cellsMove(k,2);
                                for itr=1:nInGrid(column,row)
                                    j=grid(column,row,itr);
                                    if(2*i-1~=j && 2*i~=j)
                                        if(abs(positions(2*i-1,3)-positions(j,3))<abs(move(3))+(2+1e-5)*particleRadius  && any(move) && tested.add(j) &&~touchingInfo{2*i-1}.isKey(j))
                                            [dt,tshift]=getTimeToColl(2*i-1,j,shift1,vel,move);
                                            
                                            if(dt>0 && dt<minT-1e-17)
                                                % f(j)=drawSphere(j,[0 1 0]);
                                                minT=dt;
                                                partColl=2*i-1;
                                                coll=j;
                                                fshift=tshift;
                                                nSimult1=0;
                                            elseif(dt>0 && dt<minT+1e-17)
                                                
                                                nSimult1=nSimult1+1;
                                                simultColl1(nSimult1)=j;
                                                simultShift1(nSimult1,:)=tshift;
                                                
                                            end
                                        end
                                        
                                      
                                    end
                                    %   for j=1:N
                                    %collision detection
                                    %test1=touchingInfo{i}.isKey(j);
                                    
                                    
                                    
                                end
                                
                            end
                            [cellsMove,numbC]=getCellsMoving(2*i,move);
                            tested.clear();
                            for k=1:numbC
                                column=cellsMove(k,1);
                                row=cellsMove(k,2);
                                for itr=1:nInGrid(column,row)
                                    j=grid(column,row,itr);
                                    if(2*i-1~=j && 2*i~=j)
                              if(abs(positions(2*i,3)-positions(j,3))<abs(move(3))+(2+1e-5)*particleRadius  && any(move) && tested.add(j) &&~touchingInfo{2*i}.isKey(j))
                                            [dt,tshift]=getTimeToColl(2*i,j,shift2,vel,move);
                                            
                                            if(dt>0 && dt<minT-1e-17)
                                                % f(j)=drawSphere(j,[0 1 0]);
                                                minT=dt;
                                                coll=j;
                                                partColl=2*i;
                                                fshift=tshift;
                                                nSimult2=0;
                                            elseif(dt>0 && dt<minT+1e-17)
                                                nSimult2=nSimult2+1;
                                                simultColl2(nSimult2)=j;
                                                simultShift2(nSimult2,:)=tshift;
                                                
                                            end
                              end
                                    end
                                end
                            end
                            %-----------------
                            if(positions(2*i,3)+move(3)<particleRadius) %collision with the surface
                                
                                
                                dt=(particleRadius-positions(2*i,3))/vel(3);
                                
                                if(dt<minT && dt<deltaTime)
                                    minT=dt;
                                    coll=-1;
                                    
                                end
                                
                                
                            end
                            if(positions(2*i-1,3)+move(3)<particleRadius) %collision with the surface
                                
                                
                                dt=(particleRadius-positions(2*i-1,3))/vel(3);
                                
                                if(dt<minT && dt<deltaTime)
                                    minT=dt;
                                    coll=-1;
                                    
                                end
                                
                                
                            end
                            
                            if(minT>deltaTime)
                                coll=0;
                            end
                            
                            if(coll~=0)
                                
                                % if(deposed(coll))
                                %    f(i)=drawSphere(i,[0 1 1]);
                                
                                deltaTime=deltaTime-minT;
                                positionsCenter(i,:)=positionsCenter(i,:)+vel'.*minT;
                               positionsCenter(i,:)=positionsCenter(i,:)+periodicBoundary(positionsCenter(i,:));
                                positions(i*2,:)=positions(i*2,:)+vel'.*minT;
                                positions(i*2-1,:)=positions(i*2-1,:)+vel'.*minT;
                                
                                sumOfMoves=sumOfMoves+vel'*vel.*minT/6;
                                nMoves=nMoves+1;
                                if(coll>0)
                                    %                                 numbOfColl(i)=numbOfColl(i)+1;
                                    %                                 if(numbOfColl(i)>collLimit)
                                    %                                    deposed(i)=1;
                                    %                                    break;
                                    %                                 end
                                    
                                    vel=(vel-normr((positions(coll,:)-positions(partColl,:)-fshift))'.*(normr(positions(coll,:)-positions(partColl,:)-fshift)*vel));
                                  %  vel=vel-(normr((positions(coll,:)-positions(partColl,:)-fshift)))'.*(normr(positions(coll,:)-positions(partColl,:)-fshift)*vel);
                                    
                                    touchingInfo{partColl}(coll)=fshift;
                                    if(~deposed(coll))
                                        touchingInfo{coll}(partColl)=-fshift;
                                    end
                                    
                                    if(nSimult1>0)
                                        for w=1:nSimult1
                                            touchingInfo{2*i-1}(simultColl1(w))= simultShift1(w,:);
                                            if(~deposed(simultColl1(w)))
                                                touchingInfo{simultColl1(w)}(2*i-1)=-simultShift1(w,:);
                                            end
                                        end
                                    end
                                    if(nSimult2>0)
                                        for w=1:nSimult2
                                            touchingInfo{2*i}(simultColl2(w))= simultShift2(w,:);
                                            if(~deposed(simultColl2(w)))
                                                touchingInfo{simultColl2(w)}(2*i)=-simultShift2(w,:);
                                            end
                                        end
                                    end
                                    
                                    %--------------------------sleep mode
                                    %                                 if(directionsX(i)~=sign(vel(1)))
                                    %                                     directionsX(i)=sign(vel(1));
                                    %                                     numbOfColl(i)=numbOfColl(i)+1;
                                    %                                  if(numbOfColl(i)>collLimit)
                                    %                                    deposed(i)=1;
                                    %                                    break;
                                    %                                  end
                                    %                                 end
                                    
                                    %                                 if(norm(initPos(i,:)-positions(i,:))<1e-16)
                                    %                                     deposed(i)=1;
                                    %                                     break;
                                    %                                 end
                                    
                                    %------------------------------------
                            
                                
                                else
                                    deposed(2*i)=1;
                                    deposed(2*i-1)=1;
                                    Nmoving=Nmoving-1;
                                    deltaTime=0;
                                    movingIndices(movingIndex)=[];
                                    movingIndex=movingIndex-1;
                                end
                                
                                
                                %           f(i)=drawSphere(i,[0 1 1]);
                                %   end
                            else
                                positions(2*i,:)=positions(2*i,:)+move;
                                positions(2*i-1,:)=positions(2*i-1,:)+move;
                                positionsCenter(i,:)=positionsCenter(i,:)+move;
                                positionsCenter(i,:)=positionsCenter(i,:)+periodicBoundary(positionsCenter(i,:));
                                sumOfMoves=sumOfMoves+move*move'/deltaTime/6;
                                nMoves=nMoves+1;
                                deltaTime=0;
                                %                         if(norm(move)<typicalMove*fractionOfTypical)
                                %                             notMoving(i)=notMoving(i)+1;
                                %                             if(notMoving(i)>notMovingLimit)
                                %                                 deposed(i)=1;
                                %                                 break;
                                %                             end
                                %                         else
                                %                             notMoving(i)=0;
                                %                         end
                            end
                            
                            
                            
                                 
                          
           
                            
%                
%                             if(coll>0 && deposed(coll) && stick())
%                                 deposed(i)=1;
%                                 break;
%                             end
%                         else
%
%                             notMoving(i)=notMoving(i)+1;
%                             if(notMoving(i)>notMovingLimit)
%                                 deposed(i)=1;
%                                 break;
%                             end
                    onPositionChanged(2*i-1);
                    onPositionChanged(2*i);
                    
                    for j=1:2*N
                        if(2*i~=j &&2*i-1~=j)
                            p=2*i-1;
                            if( ~touchingInfo{p}.isKey(j) && (norm(positions(p,:)-positions(j,:))-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+shift1)-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+[shift1(1),0,0])-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+[0,shift1(2),0])-2*particleRadius<0))
                                fprintf('critical error + %s\n',norm(positions(i,:)-positions(j,:))-2*particleRadius);
                            end
                             p=2*i;
                           if( ~touchingInfo{p}.isKey(j) && (norm(positions(p,:)-positions(j,:))-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+shift2)-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+[shift2(1),0,0])-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+[0,shift2(2),0])-2*particleRadius<0))
                                fprintf('critical error + %s\n',norm(positions(i,:)-positions(j,:))-2*particleRadius);
                            end
                        end
                    end

                    end
                       if(drawOn)             
                      figures(2*i-1)=drawSphere(2*i-1,[0 0 0]);   
                      figures(2*i)=drawSphere(2*i,[0 0 0]);              
                       end
                    end
                movingIndex=movingIndex+1;
            end
         %         drawFrame();
            time=time-step_size;
        end
        
            end
       
 

            function [dt,fshift]=getTimeToCollRot(i,j,posOmegaFrameI,RotMatrixWorldOmega,shift) %RotMatrixWorldOmega is problem probably
                
                dt=-1;
                fshift=[0,0,0];
               if(norm(positionsCenter(i,:)-positions(j,:))<particleRadius*3)
                                            posOmegaFrameJ=RotMatrixWorldOmega*(positions(j,:)-positionsCenter(i,:))'; 
                                            dFi=angleOfColl(posOmegaFrameJ,posOmegaFrameI);
                                            if(dFi>0)
                                            
                                                dt=dFi/norm(omegas(:,i)');
                                            
                                                
                                            end
                                        end
                                        if(any(shift))
                                            
                                            if(norm(positionsCenter(i,:)-positions(j,:)+shift)<particleRadius*3)
                                                posOmegaFrameJ=RotMatrixWorldOmega*(positions(j,:)-shift-positionsCenter(i,:))';
                                                dFi=angleOfColl(posOmegaFrameJ,posOmegaFrameI);
                                                if(dFi>0)
                                                   
                                                    dt1=dFi/norm(omegas(:,i)');
                                                    if(dt1<dt || dt<0)
                                                        dt=dt1;
                                                        fshift=shift;
                                                    end
                                                end
                                            end
                                            
                                            if(norm(positionsCenter(i,:)-positions(j,:)+[0,shift(2),0])<particleRadius*3)
                                                posOmegaFrameJ=RotMatrixWorldOmega*(positions(j,:)-[0,shift(2),0]-positionsCenter(i,:))';
                                                dFi=angleOfColl(posOmegaFrameJ,posOmegaFrameI);
                                                if(dFi>0)
                                                  dt1=dFi/norm(omegas(:,i)');
                                                    if(dt1<dt || dt<0)
                                                        dt=dt1;
                                                        fshift=[0,shift(2),0];
                                                    end
                                                end
                                            end
                                            
                                            if(norm(positionsCenter(i,:)-positions(j,:)+[shift(1),0,0])<particleRadius*3)
                                                 posOmegaFrameJ=RotMatrixWorldOmega*(positions(j,:)-[shift(1),0,0]-positionsCenter(i,:))';
                                                dFi=angleOfColl(posOmegaFrameJ,posOmegaFrameI);
                                                if(dFi>0)
                                                   dt1=dFi/norm(omegas(:,i)');
                                                    if(dt1<dt || dt<0)
                                                        dt=dt1;
                                                        fshift=[shift(1),0,0];
                                                    end
                                                end
                                            end 
                
                
                
                
            end
            end
        function onPositionChanged(i)
        
                       perShift=periodicBoundary(positions(i,:));
                            positions(i,:)=positions(i,:)+perShift;
                            keys=touchingInfo{i}.keys();
                            
                            
                            for q=1:length(keys)
                                if(any(perShift))
                                    
                                    touchingInfo{i}(keys{q})=touchingInfo{i}(keys{q})-perShift;
                                    if(~deposed(keys{q}))
                                        touchingInfo{keys{q}}(i)=touchingInfo{keys{q}}(i)+perShift;
                                    end
                                    
                                end
                                if(norm(positions(keys{q},:)-positions(i,:)-touchingInfo{i}(keys{q}))>(2+1e-5)*particleRadius)
                                    touchingInfo{i}.remove(keys{q});
                                    if(~deposed(keys{q}))
                                        touchingInfo{keys{q}}.remove(i);
                                    end
                                end
                                
                            end
                            %-----------------space partition
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
        
        end
        
       
        function dFi=angleOfColl(pos,pos0)%always returns positive angle if collisions occurs
            c=(pos(1)^2+pos(2)^2+pos0(1)^2+pos0(2)^2+(pos(3)-pos0(3))^2-4*particleRadius^2)/2/sqrt(pos0(1)^2+pos0(2)^2)/sqrt(pos(1)^2+pos(2)^2);
            dFi=atan2(pos(2),pos(1))-acos(c)-atan2(pos0(2),pos0(1));
            
            if(isreal(dFi))
               % fprintf('%f\n',norm([cos(dFi),-sin(dFi),0;sin(dFi),cos(dFi),0;0,0,1]*pos0-pos)-2*particleRadius);
                if(dFi<0)
                    dFi=dFi+2*pi;
                end
                
            else
                dFi=-1;
            end
        end
        function dFi=getAngleOfCollGround(pos0,zc,n)
           dFi=asin((particleRadius+zc)/sqrt(pos0(1)^2+pos0(2)^2)/sqrt(n(1)^2+n(2)^2))-atan2(n(1),n(2))-atan2(pos0(2),pos0(1));
           if(isreal(dFi))
              if(dFi<0)
                  dFi=dFi+2*pi; 
              end
           else
               dFi=-1;
           end
        end
      
        
        function [minT,fshift]=getTimeToColl(i,j,shift,vel,move)
            minT=-1;
          
                      fshift=[0 0 0];                  
                                        
                                       
                                        if(norm(positions(i,:)-positions(j,:))<norm(move)+particleRadius*2)
                                            delta_r=positions(j,:)-positions(i,:);%*1e5;
                                            %vel=vel;.*1e5;
                                            %particleRadius=particleRadius*1e5;
                                            
                                            discr=(delta_r*vel)^2-vel'*vel*(delta_r*delta_r'-(2*particleRadius)^2);
                                            if(discr>=0)
                                                dt=(delta_r*vel-sqrt(discr))/(vel'*vel);
                                                
                                                
                                               if(dt>0 && (dt<minT|| minT<0))
                                                    % f(j)=drawSphere(j,[0 1 0]);
                                                    minT=dt;
                                                  
                                                    fshift=[0,0,0];
                                                  
                                           
                                                    
                                                end
                                                
                                            end
                                            %vel=vel.*1e-5;
                                            % particleRadius=particleRadius*1e-5;
                                            
                                            %   if(deposed(j)) %collision with a fixed particle
                                            %      resolveCollision(i,j);
                                            %                             if(positions(j,3)<=particleRadius)
                                            %                                 resolveCollision(i,j);
                                            %                             else %collision with a moving particle
                                            %                                 dCenters=positions(j,:)-positions(i,:);
                                            %                                 positions(i,:)=positions(i,:)-0.5*(2*particleRadius-norm(dCenters))*normr(dCenters);
                                            %                                 positions(j,:)=positions(j,:)+0.5*(2*particleRadius-norm(dCenters))*normr(dCenters);
                                            %
                                            %                             end
                                            
                                            %TODO REPLACE ELSEIF WITH IF
                                        end
                                        
                                        
                                            
                                            
                                            
                                            if(shift(1)~=0 && shift(2)~=0 )
                                                delta_r=positions(j,:)-positions(i,:)-shift;
                                                disc=discriminant(delta_r,vel);
                                                if(disc>=0)
                                                dt=(delta_r*vel-sqrt(disc))/(vel'*vel);
                                               
                                                if(dt>0 && (dt<minT|| minT<0))
                                                    fshift=shift;
                                                    minT=dt;
                                                  
                                                end
                                                end
                                            end
                                            
                                            if(shift(1)~=0 )
                                                delta_r=positions(j,:)-positions(i,:)-[shift(1),0,0];
                                                disc=discriminant(delta_r,vel);
                                                if(disc>=0)
                                                    dt=(delta_r*vel-sqrt(disc))/(vel'*vel);
                                                    if(dt>0 && (dt<minT|| minT<0))
                                                        fshift=[shift(1),0,0];
                                                        minT=dt;
                                                    end
                                                end
                                            end
                                            
                                            
                                            if(shift(2)~=0)
                                                delta_r=positions(j,:)-positions(i,:)-[0,shift(2),0];
                                                disc=discriminant(delta_r,vel);
                                                if(disc>=0)
                                                    dt=(delta_r*vel-sqrt(disc))/(vel'*vel);
                                                    if(dt>0 && (dt<minT|| minT<0))
                                                        fshift=[0,shift(2),0];
                                                        minT=dt;
                                                    end
                                                end
                                                
                                            end
                                          
                                        
                                        
                                    end 
            
            
            
            
        
        
        
    function d=discriminant(delta_r,vel)
        d=(delta_r*vel)^2-vel'*vel*(delta_r*delta_r'-(2*particleRadius)^2);
    end

    function cells=get16Cells(pos)
        %16 closest cells to a particle
        cells=zeros(16,2);
        cells4(1:4,:)=getCells(pos);
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
        %     l=ceil((pos(3))/deltaZ);
        %              m1=m+1;
        %              m2=m-1;
        %              n1=n+1;
        %              n2=n-1;
        
        m1=round((pos(1)-xMin)./deltaX);
        n1=round((pos(2)-yMin)./deltaY);
        %   l1=round(pos(3)/deltaZ);
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

% function newMove=fixTouching (move,i)
% keys=touchingInfo{i}.keys();
% len=length(keys);
% newMove=move;
% if(len>0)
% 
%     dir=normr(positions(keys{1},:)-positions(i,:)-touchingInfo{i}(keys{1}));
%     if(dir*move'>0)
%         newMove=move-dir.*(dir*move');
%     else
%         newMove=move;
%     end
%     if(len>1)
%         moveWrong=0;
%         pairDir.remove(pairDir.keys());
%         rollMove=[0 0 0];
%         index=0;
%         for j=1:len
%             dir1=normr(positions(keys{j},:)-positions(i,:)-touchingInfo{i}(keys{j}));
%             %---rolling over two
%             %                         if(j==len && len>2)
%             %                             index=index+1;
%             %                             pairDir(index)=normr(cross(dir1,dir));
%             %                         else
%             for k=j+1:len
%                 index=index+1;
%                 dir2=positions(keys{k},:)-positions(i,:)-touchingInfo{i}(keys{k});
%                 pairDir(index)=normr(cross(dir1,dir2));
%             end
%             %  end
%             %----------------
%             %-----rolling over one
%             if(dir1*move'>0 )
%                 if(~any(rollMove))
%                     rollMove=move-dir1*(dir1*move');
%                 else
%                     rollMove=[0 0 0];
%                 end
%             end
%             %--------------------
%         end
% 
% 
%         for c=1:len
%             dir1=positions(keys{c},:)-positions(i,:)-touchingInfo{i}(keys{c});
%             dirKeys=pairDir.keys();
%             for k=1:length(dirKeys)
%                 checkDir=pairDir(dirKeys{k});
%                 if((checkDir*(move*checkDir'))*dir1'>1e-16)
%                     pairDir.remove(dirKeys(k));
% 
%                 end
%             end
% 
% 
%             if(any(rollMove) && rollMove*dir1'>1e-16)
%                 rollMove=[0 0 0];
%             end
% 
%             if(dir1*move'>1e-16)
%                 moveWrong=1;
%             end
% 
%         end
%         dirKeys=pairDir.keys();
%         if(moveWrong)
%             if(any(rollMove))
%                 newMove=rollMove;
%             elseif(~isempty(dirKeys))
%                 newMove=(pairDir(dirKeys{1})*move')*pairDir(dirKeys{1});
%             else
%                 newMove=[0 0 0];
%             end
%         else
%             newMove=move;
%         end
%     end
% 
% end
% end
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