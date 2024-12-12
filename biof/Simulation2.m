function [ cov] = Simulation2( step_size,time,rate)
%Simulation of spherical particle deposition

tic; %start timer





persistent particleRadius;
persistent xMax; %box corners
persistent xMin;
persistent yMin;
persistent yMax;
persistent N; %number of created particles

xMax=5;
xMin=-xMax;
yMax=5e-0;
yMin=-yMax;

zMax=10e-0;
zMin=zMax;

notMovingLimit=3;
numbOfColl=zeros(time*rate,3);
collLimit=200;
directionsX=zeros(1,time*rate);


positions=zeros(time*rate,3);


deposed=zeros(1,time*rate);

pairDir=containers.Map('KeyType','int32','ValueType','any');
particleRadius=1e-0;
typicalMove=0;
fractionOfTypical=0.001;
virtFigures=gobjects(1,time*rate);

touchingInfo=cell(1,time*rate);
for o=1:time*rate
    touchingInfo{o}=containers.Map('KeyType','int32','ValueType','any');
end

notMoving=zeros(1,time*rate);

nColumns=floor((xMax-xMin)/2/particleRadius);
nRows=floor((yMax-yMin)/2/particleRadius);

nLayers=floor(zMax/2/particleRadius);
% grid(nColumns,nRows)=java.util.HashSet;
% for a=1:nColumns
%     for b=1:nRows
%         grid(a,b)=java.util.HashSet;
%     end
% end

grid=zeros(nColumns,nRows,ceil(time*rate/nColumns/nRows*4));
nInGrid=zeros(nColumns,nRows);

%grid=zeros(floor((xMax-xMin)/particleRadius),floor((yMax-yMin)/particleRadius),8);

inCells=zeros(4,2,time*rate); %every particle is in up to 4 cells of the grid
deltaX=(xMax-xMin)/nColumns;
deltaY=(yMax-yMin)/nRows;
deltaZ=zMax/nLayers;
tested=java.util.HashSet;
added=java.util.HashSet;

%velocity=[0;0;-50];

f=gobjects(time*rate); %uncomment to draw each frame

N=0; %number of created particles

calcCov=1; %1 to calculate surface coverage as a function of time

simultColl=zeros(1,3);
simultShift=zeros(3,3);
%---------------
T=300; %K
kb=1.38e-23;
d=2*particleRadius;
eta=1e-3;
mass=1e-0;
dens=mass/(4/3*pi*particleRadius^3);
%gforce=mass*(1-1000/dens)*[0,0,-9.8];
gforce=mass*[0,0,-9.8];
stickProbability=0.0;
epsilon0=8.85e-12;
epsilonR=78;
charge=1.6e-19;
numberDensity=1e-7*6.02e29;
DebLength=sqrt(charge^2*numberDensity/epsilon0/epsilonR/kb/T);
potential=1e-3;
elForceConst=2*pi*epsilon0*epsilonR*DebLength*particleRadius*potential^2;

numbDensPart=3e26;
lje=1e3/6.02e23;
ljs=particleRadius/5;
Acc=1e-20;%4*pi^2*lje*numbDensPart^2*ljs^6;

%---------------
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

for t=1:nIntervals
    
    
    positions(t,:)=[xMin+rand()*(xMax-xMin);yMin+rand()*(yMax-yMin);zMin+rand()*(zMax-zMin)];
    %sp(t,:)=positions(t,:);
    %positions(t,:)=sp(t,:);
    N=N+1;
    %physics simulation
    
    
    
    inCells(:,:,N)=getCells(positions(N,:));
    for cel=1:4
        addToGrid(inCells(cel,1,N),inCells(cel,2,N),N);
    end
    
    
    physicsStep(interval);
    
    
    %             while (any(particle.velocity))
    %                 if(checkForCollision(particle))
    %                     particle.collide();
    %                     delLastParticle=0;
    %                 end
    %                 if (particle.position(3)<particle.radius)
    %                     particle.collide();
    %                     particle.deposed=1;
    %                     delLastParticle=0;
    %                 end
    %                 [dx,dy,dz]=BrownianMotion(step_size,T,eta,d,fluidVel);
    %                 particle.position=particle.position+[dx;dy;dz];
    %                 %  lastParticle=drawSphere(particle.position);
    %                 %  delLastParticle=1;
    %             end
    
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
%volume fraction as a function of z
zTop=max(positions(:,3));
zStep=zTop/nOfZSteps;
v=zeros(1,nOfZSteps);
zVolF=linspace(0,zTop,nOfZSteps);
for s=1:nOfZSteps
    v(s)=volumeFraction(zStep*(s-1)+particleRadius,zStep);
end

figure();
plot(zVolF,v);

fprintf('time %.3f s\n',toc);

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
            f(i)=fig;
        end
        hold off;
    end
    function drawFrame()
        for i=1:N
            f(i)=drawSphere(i,[0 0 0]);
        end
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

    function fig=drawSphere(p,c)
        if(~isempty(f))
            try
                delete(f(p));
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

    function f=elecForce(i)
        f=[0 0 0];
       
        cells=get16Cells(i);
        added.clear();
        for c=1:16
            
         %   iter=grid(cells(c,1),cells(c,2)).iterator();
            for it=1:nInGrid(cells(c,1),cells(c,2))
                j=grid(cells(c,1),cells(c,2),it);
                if(i~=j && abs(positions(i,3)-positions(j,3))<=4*particleRadius&& added.add(j))
                    shift=[0 0 0];
                    deltaR=positions(i,:)-positions(j,:);
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
                    deltaR=deltaR+shift;
                    if(norm(deltaR)<=4*particleRadius)
                        
                        
                        
                        
                        r=norm(deltaR);
                        dir=deltaR/r;
                        fel=elForceConst/(1+exp(DebLength*(r-2*particleRadius)))*dir;
                        
                        if(r>2*particleRadius+30^(-1/6)*ljs) %recommended cut-off distance (potential is divergent at r=2*radius)
                            fDisp1=-32/3*Acc*particleRadius^6/r^3/(r^2-4*particleRadius^2)^2; %r^-6 term of LJ potential
                            fDisp2=Acc/37800*ljs^6*((56*particleRadius*r^2+10*particleRadius^3-6*r^3)/(r-2*particleRadius)^8/r^2-(r^3+26*particleRadius*r^2+108*particleRadius^3+432*particleRadius^2*r)/(r+2*particleRadius)^8/r^2+12*(r^2-20*particleRadius^2)/r^9);
                            
                        else
                            fDisp1=0;
                            fDisp2=0;
                        end
                        f=f+fel+(fDisp1+fDisp2)*dir;
                    end
                end
            end
        end
    end
    
    function newVel=fixTouching(i,vel)
%         function b=isCollSoved()
%            for c=1:n
%                if(
%         end
        newVel=vel;
        keys=touchingInfo{i}.keys();
       n=length(keys);
       impulses=zeros(n,1);
       contactSpeed=zeros(n,1);
       contactNorms=zeros(n,3);
       collNotSolved=0;
       A=zeros(n);
       indices=[1:n];
       C=[];
       NC=[];
       for c=1:n
           dir=normr(positions(keys{c},:)-positions(i,:)-touchingInfo{i}(keys{c}));
           contactNorms(c,:)=-dir;
           contactSpeed(c)=dir*vel;
           if(contactSpeed(c)>0)
                collNotSolved=1;
                d0=c;
            end
       end
       for a=1:n
           for b=a:n
               A(a,b)=contactNorms(a)*contactNorms(b)';
               A(b,a)=contactNorms(a)*contactNorms(b)';
           end
       end
       while (collNotSolved)
           collNotSolved=0;
           driveToZero(d0);
           newVel=vel;
           for c=1:n
               
               if(contactSpeed(c)>0)
                   collNotSolved=1;
                   d0=c;
                  
               end
               newVel=newVel+impulses(n)*contactNorms(n,:)';
           end
       end
       
        function driveToZero(d)
            deltaImpulse=fdirection(d);
            deltaContSpeed=A*deltaImpulse;
            [s,j]=maxStep(impulses,contactSpeed,deltaImpulse,d);
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
           indsToRemove=indices;
           indsToRemove(C)=[];
           A11=A;
           A11(:,indsToRemove)=[];
           A11(indsToRemove,:)=[];
           v1=A(:,d);
           v1(indsToRemove)=[];
           x=-A11\v1;
           for u=1:length(C)
               deltaF(C(u))=x(C(u));
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
                if(deltaF(q)<0)
                    s1=-f(q)/deltaF(q);
                    if(s1<s)
                        s=s1;
                        j=q;
                    end
                end
            end
            for q=1:length(NC)
                if(deltaA(q)<0)
                    s1=-a(q)/deltaA(q);
                    if(s1<s)
                        s=s1;
                        j=q;
                    end
                end
            end
        end
    end

    function physicsStep(time)
        %    initPos=positions;
        while (time>0)
            
            for i=1:N
                if(~deposed(i))
                    [dx,dy,dz]=BrownianMotion(step_size,T,eta,d,gforce+elecForce(i));
                    move=[dx,dy,dz];
                    if(typicalMove==0)
                        typicalMove=norm(move);
                    end
                    
                    %  move=move+periodicBoundary(i);
                    
                    %  positions(i,:)=positions(i,:)+move;
                    
                    
                    
                    deltaTime=step_size;
                    numberOfIter=0;
                    
                    while(deltaTime>0 && any(move) && numberOfIter<10)% && minT>step_size/1e6)
                        
                        numberOfIter=numberOfIter+1;
                        vel=move'./deltaTime;
                        shift=[0,0,0];
                        minT=step_size;
                        coll=0; %is there a collision
                        fshift=[0,0,0];
                        nSimult=0;
                        %------------touching
%                         keys=touchingInfo{i}.keys();
%                         len=length(keys);
%                         if(len>0)
%                             
%                             dir=normr(positions(keys{1},:)-positions(i,:)-touchingInfo{i}(keys{1}));
%                             if(dir*move'>0)
%                                 newMove=move-dir.*(dir*move');
%                             else
%                                 newMove=move;
%                             end
%                             if(len>1)
%                                 moveWrong=0;
%                                 pairDir.remove(pairDir.keys());
%                                 rollMove=[0 0 0];
%                                 index=0;
%                                 for j=1:len
%                                     dir1=normr(positions(keys{j},:)-positions(i,:)-touchingInfo{i}(keys{j}));
%                                     %---rolling over two
%                                     %                         if(j==len && len>2)
%                                     %                             index=index+1;
%                                     %                             pairDir(index)=normr(cross(dir1,dir));
%                                     %                         else
%                                     for k=j+1:len
%                                         index=index+1;
%                                         dir2=positions(keys{k},:)-positions(i,:)-touchingInfo{i}(keys{k});
%                                         pairDir(index)=normr(cross(dir1,dir2));
%                                     end
%                                     %  end
%                                     %----------------
%                                     %-----rolling over one
%                                     if(dir1*move'>0 )
%                                         if(~any(rollMove))
%                                             rollMove=move-dir1*(dir1*move');
%                                         else
%                                             rollMove=[0 0 0];
%                                         end
%                                     end
%                                     %--------------------
%                                 end
%                                 
%                                 
%                                 for c=1:len
%                                     dir1=positions(keys{c},:)-positions(i,:)-touchingInfo{i}(keys{c});
%                                     dirKeys=pairDir.keys();
%                                     for k=1:length(dirKeys)
%                                         checkDir=pairDir(dirKeys{k});
%                                         if((checkDir*(move*checkDir'))*dir1'>1e-16)
%                                             pairDir.remove(dirKeys(k));
%                                             
%                                         end
%                                     end
%                                     
%                                     
%                                     if(any(rollMove) && rollMove*dir1'>1e-16)
%                                         rollMove=[0 0 0];
%                                     end
%                                     
%                                     if(dir1*move'>1e-16)
%                                         moveWrong=1;
%                                     end
%                                     
%                                 end
%                                 dirKeys=pairDir.keys();
%                                 if(moveWrong)
%                                     if(any(rollMove))
%                                         newMove=rollMove;
%                                     elseif(~isempty(dirKeys))
%                                         newMove=(pairDir(dirKeys{1})*move')*pairDir(dirKeys{1});
%                                     else
%                                         newMove=[0 0 0];
%                                     end
%                                 else
%                                     newMove=move;
%                                 end
%                             end
% %                             
% %                             
%                            
%                             %------------------------------
%                         end
                        
                         move=fixTouching(i,move')';
                            
                            keys=touchingInfo{i}.keys();
                            for c=1:length(keys)
                                dir1=positions(keys{c},:)-positions(i,:)-touchingInfo{i}(keys{c});
                                if(dir1*move'>1e-15)
                                    fprintf('asdf\n');
                                end
                            end
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
                            
                            %   plot3([positions(i,1),move(1)+positions(i,1)],[positions(i,2),move(2)+positions(i,2)],[positions(i,3),move(3)+positions(i,3)]);
                            tested.clear();
                            [cellsMove,numbC]=getCellsMoving(i,move);
                            
                            for k=1:numbC
                                column=cellsMove(k,1);
                                row=cellsMove(k,2);
                                for it=1:nInGrid(column,row)
                                   j=grid(column,row,it);
                                    
                                    %   for j=1:N
                                    %collision detection
                                    %test1=touchingInfo{i}.isKey(j);
                                    
                                    if(i~=j && abs(positions(i,3)-positions(j,3))<abs(move(3))+2*particleRadius  && any(move) && tested.add(j) &&~touchingInfo{i}.isKey(j))
                                        
                                        
                                        primChck=0; %primitive check of collisions through the boundaries
                                        if(norm(positions(i,:)-positions(j,:))<norm(move)+particleRadius*2)
                                            delta_r=positions(j,:)-positions(i,:);%*1e5;
                                            %vel=vel;.*1e5;
                                            %particleRadius=particleRadius*1e5;
                                            
                                            discr=(delta_r*vel)^2-vel'*vel*(delta_r*delta_r'-(2*particleRadius)^2);
                                            if(discr>=0)
                                                dt=(delta_r*vel-sqrt(discr))/(vel'*vel);
                                                
                                                if(dt<=0 && norm(positions(i,:)-positions(j,:))<=2*particleRadius &&(positions(j,:)-positions(i,:))*vel>0 )
                                                    fprintf('this should not have happened');
                                                elseif(dt>0 && dt<deltaTime&& dt<minT-1e-17)
                                                    % f(j)=drawSphere(j,[0 1 0]);
                                                    minT=dt;
                                                    coll=j;
                                                    fshift=[0,0,0];
                                                    nSimult=0;
                                                elseif(dt>0 && dt<deltaTime&& dt<minT+1e-17)
                                                    nSimult=nSimult+1;
                                                    simultColl(nSimult)=j;
                                                    simultShift(nSimult,:)=[0,0,0];
                                                    
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
                                        if(any(shift))
                                            delta_r=positions(j,:)-positions(i,:)-shift;
                                            disc=discriminant(delta_r,vel);
                                            
                                            
                                            if(shift(1)~=0 && shift(2)~=0 && disc>=0)
                                                dt=(delta_r*vel-sqrt(disc))/(vel'*vel);
                                                if(dt<=0 && norm(positions(i,:)-positions(j,:)+shift)<=2*particleRadius &&(positions(j,:)-positions(i,:)-shift)*vel>0 )
                                                    fprintf('this should not have happened');
                                                elseif(dt>0 && dt<deltaTime)
                                                    tShift=shift;
                                                    
                                                    primChck=1;
                                                end
                                            end
                                            
                                            if(~primChck && shift(1)~=0 )
                                                delta_r=positions(j,:)-positions(i,:)-[shift(1),0,0];
                                                disc=discriminant(delta_r,vel);
                                                if(disc>=0)
                                                    dt=(delta_r*vel-sqrt(disc))/(vel'*vel);
                                                    if(dt<=0 && norm(positions(i,:)-positions(j,:)+[shift(1),0,0])<=2*particleRadius &&(positions(j,:)-positions(i,:)-[shift(1),0,0])*vel>0 )
                                                        fprintf('this should not have happened');
                                                    elseif(dt>0 && dt<deltaTime)
                                                        tShift=[shift(1),0,0];
                                                        primChck=1;
                                                    end
                                                end
                                            end
                                            
                                            
                                            if(~primChck && shift(2)~=0)
                                                delta_r=positions(j,:)-positions(i,:)-[0,shift(2),0];
                                                disc=discriminant(delta_r,vel);
                                                if(disc>=0)
                                                    dt=(delta_r*vel-sqrt(disc))/(vel'*vel);
                                                    if(dt<=0 && norm(positions(i,:)-positions(j,:)+[0,shift(2),0])<=2*particleRadius &&(positions(j,:)-positions(i,:)-[0,shift(2),0])*vel>0 )
                                                        fprintf('this should not have happened');
                                                    elseif(dt>0 && dt<deltaTime )
                                                        tShift=[0,shift(2),0];
                                                        primChck=1;
                                                    end
                                                end
                                                
                                            end
                                            if (primChck)
                                                %                                     delta_r=positions(j,:)-positions(i,:)-tShift;
                                                %                                     discr=(delta_r*vel)^2-vel'*vel*(delta_r*delta_r'-(2*particleRadius)^2);
                                                
                                                
                                                if(dt<=0 && norm(positions(i,:)-positions(j,:)+tShift)<=2*particleRadius &&(positions(j,:)-positions(i,:)-tShift)*vel>0 )
                                                    fprintf('this should not have happened');
                                                    %TODO
                                                    
                                                elseif ( dt<minT-1e-17)
                                                    % f(j)=drawSphere(j,[0 1 0]);
                                                    minT=dt;
                                                    coll=j;
                                                    fshift=tShift;
                                                    nSimult=0;
                                                elseif(dt<minT+1e-17)
                                                    nSimult=nSimult+1;
                                                    simultColl(nSimult)=j;
                                                    simultShift(nSimult,:)=tShift;
                                                    
                                                end
                                                
                                            end
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                            if(positions(i,3)+move(3)<particleRadius) %collision with the surface
                                %  vel=vel.*1e5;
                                %  particleRadius=particleRadius*1e5;
                                
                                dt=(particleRadius-positions(i,3))/vel(3);
                                %   vel=vel.*1e-5;
                                %    particleRadius=particleRadius*1e-5;
                                if(dt<minT && dt<deltaTime)
                                    minT=dt;
                                    coll=-1;
                                    deposed(i)=1;
                                    deltaTime=0;
                                    %  deltaTime=0;
                                end
                                %  deposed(i)=1;
                                
                            end
                            
                            
                            
                            if(coll~=0)
                                
                                % if(deposed(coll))
                                %    f(i)=drawSphere(i,[0 1 1]);
                                
                                deltaTime=deltaTime-minT;
                                positions(i,:)=positions(i,:)+vel'.*minT;
                                
                                if(coll>0)
                                    %                                 numbOfColl(i)=numbOfColl(i)+1;
                                    %                                 if(numbOfColl(i)>collLimit)
                                    %                                    deposed(i)=1;
                                    %                                    break;
                                    %                                 end
                                    
                                    move=move-normr((positions(coll,:)-positions(i,:)-fshift)).*(normr(positions(coll,:)-positions(i,:)-fshift)*move');
                                    
                                    
                                    touchingInfo{i}(coll)=fshift;
                                    if(~deposed(coll))
                                        touchingInfo{coll}(i)=-fshift;
                                    end
                                    
                                    if(nSimult>0)
                                        for w=1:nSimult
                                            touchingInfo{i}(simultColl(w))= simultShift(w,:);
                                            if(~deposed(simultColl(w)))
                                                touchingInfo{simultColl(w)}(i)=-simultShift(w,:);
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
                                  
                                end
                                
                                
                                %           f(i)=drawSphere(i,[0 1 1]);
                                %   end
                            else
                                positions(i,:)=positions(i,:)+move;
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
                            perShift=periodicBoundary(positions(i,:));
                            positions(i,:)=positions(i,:)+perShift;
                            
                                  f(i)=drawSphere(i,[0 0 0]); %draw each frame
                            %                     for j=1:N
                            %                         if(i~=j && (norm(positions(i,:)-positions(j,:))-2*particleRadius<-particleRadius/100 || norm(positions(i,:)-positions(j,:)+shift)-2*particleRadius<-particleRadius/100 || norm(positions(i,:)-positions(j,:)+[shift(1),0,0])-2*particleRadius<-particleRadius/100 || norm(positions(i,:)-positions(j,:)+[0,shift(2),0])-2*particleRadius<-particleRadius/100))
                            %                             fprintf('critical error + %s\n',norm(positions(i,:)-positions(j,:))-2*particleRadius);
                            %                         end
                            %                     end
                            %
                            %                      for j=1:N
                            %                         if(i~=j && ~touchingInfo{i}.isKey(j) && (norm(positions(i,:)-positions(j,:))-2*particleRadius<0 || norm(positions(i,:)-positions(j,:)+shift)-2*particleRadius<0 || norm(positions(i,:)-positions(j,:)+[shift(1),0,0])-2*particleRadius<0 || norm(positions(i,:)-positions(j,:)+[0,shift(2),0])-2*particleRadius<0))
                            %                             fprintf('critical error + %s\n',norm(positions(i,:)-positions(j,:))-2*particleRadius);
                            %                         end
                            %                     end
                            keys=touchingInfo{i}.keys();
                            
                            
                            for q=1:length(keys)
                                if(any(perShift))
                                    
                                    touchingInfo{i}(keys{q})=touchingInfo{i}(keys{q})-perShift;
                                    if(~deposed(keys{q}))
                                        touchingInfo{keys{q}}(i)=touchingInfo{keys{q}}(i)+perShift;
                                    end
                                    
                                end
                                if(norm(positions(keys{q},:)-positions(i,:)-touchingInfo{i}(keys{q}))>(2+1e-5)*particleRadius && coll~=keys{q})
                                    touchingInfo{i}.remove(keys{q});
                                    if(~deposed(keys{q}))
                                        touchingInfo{keys{q}}.remove(i);
                                    end
                                end
                                
                            end
                            
                            %                     for q=1:numberOfTouch(i)
                            %                         if(any(perShift))
                            %                             touchingShift(q,:,i)=touchingShift(q,:,i)-perShift;
                            %                         end
                            %                         if(norm(positions(touchingInd(i,q),:)-positions(i,:)-touchingShift(q,:,i))>(2+1e-10)*particleRadius && coll~=touchingInd(i,q))
                            %                             numberOfTouch(i)=numberOfTouch(i)-1;
                            %                             touchingInd(i,q)=0;
                            %                             touchingShift(q,:,i)=[0,0,0];
                            %                             if(q>1 && q<12)
                            %                                 touchingInd(i,:)=[0 touchingInd(i,1:q-1) touchingInd(i,q+1:12)];
                            %                                 touchingShift(:,:,i)=[[0,0,0];touchingShift(1:q-1,:,i); touchingShift(q+1:12,:,i)];
                            %                             end
                            %                         end
                            %                     end
                            %                     while(touchingInd(i,1)==0 && any(touchingInd(i,:)))
                            %                         touchingInd(i,:)=circshift(touchingInd(i,:)',11)';
                            %                         touchingShift(:,:,i)=circshift(touchingShift(:,:,i),11);
                            %                     end
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
                            if(coll>0 && deposed(coll) && stick())
                                deposed(i)=1;
                                break;
                            end
                        else
                            
                            notMoving(i)=notMoving(i)+1;
                            if(notMoving(i)>notMovingLimit)
                                deposed(i)=1;
                                break;
                            end
                            
                        end
                        
                    end
                    
                end
            end
                 % drawFrame();
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