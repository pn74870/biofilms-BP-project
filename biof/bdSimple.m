function cov=bdSimple(time,rate)
%Simulation of spherical particle deposition

tic; %start timer
v=0;





persistent particleRadius;
persistent xMax;
persistent xMin;
persistent yMin;
persistent yMax;
persistent N;
persistent pos;

xMax=2e-0;
xMin=-xMax;
yMax=2e-0;
yMin=-yMax;

zMax=5e-1;
zMin=5e-1;



positions=zeros(time*rate,2);
deposed=zeros(1,time*rate);



particleRadius=1e-1;
virtFigures=gobjects(1,1);
%--------del---------------
% touchingShift=zeros(12,3,time*rate);
% numberOfTouch=zeros(1,time*rate);
% touchingInd=zeros(time*rate,12);
%------------------------
touchingInfo=cell(1,time*rate);
for o=1:time*rate
    touchingInfo{o}=containers.Map('KeyType','int32','ValueType','any');
end



nColumns=floor((xMax-xMin)/2/particleRadius);
nRows=floor((yMax-yMin)/2/particleRadius);
grid(nColumns,nRows)=java.util.HashSet;
for a=1:nColumns
    for b=1:nRows
        grid(a,b)=java.util.HashSet;
    end
end
pairDir=containers.Map('KeyType','int32','ValueType','any');

%grid=zeros(floor((xMax-xMin)/particleRadius),floor((yMax-yMin)/particleRadius),8);

inCells=zeros(4,2,time*rate); %every particle is in up to 4 cells of the grid
deltaX=(xMax-xMin)/nColumns;
deltaY=(yMax-yMin)/nRows;

tested=java.util.HashSet;
 

%velocity=[0;0;-50];

f=gobjects(time*rate); %uncomment to draw each frame

N=0; %number of created particles

calcCov=1; %1 to calculate surface coverage as a function of time

simultColl=zeros(1,3);
simultShift=zeros(3,3);
%---------------
T=300; %K
d=2*particleRadius;
eta=1e-0;
mass=1e-0;
gforce=mass*[0,0,-9.8];
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
coveredArea=0;
totalArea=(xMax-xMin)*(yMax-yMin);
interval=1/rate;
nIntervals=time*rate;

fixed=containers.Map('KeyType','int32','ValueType','any');
for t=1:nIntervals
    
    pos=[xMin+rand()*(xMax-xMin),yMin+rand()*(yMax-yMin)];
    


    fixed.remove(fixed.keys());

    if(resolveCollisions())
        N=N+1;
        positions(N,:)=pos;
        coveredArea=coveredArea+pi*particleRadius^2;
        f(N)=drawPart(N);
    end
         %sp(t,:)=positions(t,:);
        %positions(t,:)=sp(t,:);
        
        %physics simulation
        
      
        
        
        
       % g=zeros(numbOfSteps,numbOfSteps);
 
        
%         
%      
%                     
%                     inCells(:,:,N)=getCells(positions(N,:));
%                     for cel=1:4
%                          grid(inCells(cel,1,N),inCells(cel,2,N)).add(N);
%                     end
%         
%         
%         physicsStep(interval);
            
            
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
        cov(t/rate)=coveredArea/totalArea;
    end
    fprintf('%.1f percent completed\n',t/nIntervals*100);
end
delete(virtFigures);
drawAll();
hold off;
if (calcCov)
    figure();
    plot(cov);
end
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

    function b=isObtuseTriangle(a,b,c)
        r1=a-b;
        r2=c-b;
        r3=a-c;
        r4=-r2;
        r5=-r1;
        r6=-r3;
        
        b=r1*r2'<0 || r3*r4'<0 || r5*r6'<0;
        
    end
 
    function deposed=resolveCollisions()
    pShift=periodicBoundary(pos);
    pos=pos+pShift;
    keys=fixed.keys();
    numbOfFix=length(keys);
    if(any(pShift) && numbOfFix>0)
        for i=1:numbOfFix
            fixed(keys{i})=fixed(keys{i})-pShift;
        end
    end
   shift=[0,0];
   shiftsColl=zeros(6,2);
   if(pos(1)<xMin+2*particleRadius)
       shift(1)=xMax-xMin;
   
   elseif(pos(1)>xMax-2*particleRadius)
       shift(1)=-xMax+xMin;
   end

    if(pos(2)<yMin+2*particleRadius)
       shift(1)=yMax-yMin;
   
   elseif(pos(2)>yMax-2*particleRadius)
       shift(1)=-yMax+yMin;
    end
        
    virtFigures=drawSphere(pos);
    collisions=zeros(1,6);
    numbOfColl=numbOfFix;
    
    
    for i=1:numbOfFix
        collisions(i)=keys{i};
        shiftsColl(i,:)=fixed(keys{i});
    end
    for i=1:N
        if(~any(collisions==i))
        if((norm(pos-positions(i,:))-2*particleRadius<0))
            numbOfColl=numbOfColl+1;
            collisions(numbOfColl)=i;
            shiftsColl(numbOfColl,:)=[0 0];
        elseif(any(shift) && (norm(pos+shift-positions(i,:))-2*particleRadius<0))
            numbOfColl=numbOfColl+1;
            collisions(numbOfColl)=i;
            shiftsColl(numbOfColl,:)=shift;
        elseif (any(shift) &&(norm(pos+[shift(1),0]-positions(i,:))-2*particleRadius<0))
            numbOfColl=numbOfColl+1;
            collisions(numbOfColl)=i;
            shiftsColl(numbOfColl,:)=[shift(1),0];
            
        elseif (any(shift) && (norm(pos+[0, shift(2)]-positions(i,:))-2*particleRadius<0))
            numbOfColl=numbOfColl+1;
            collisions(numbOfColl)=i;
            shiftsColl(numbOfColl)=[0, shift(2)];
            
        end
        end
    end 
    if(numbOfColl==0 || (numbOfColl==3 && fixed.isKey(collisions(1)) && fixed.isKey(collisions(2)) && fixed.isKey(collisions(3))) || (numbOfColl==2 && fixed.isKey(collisions(1)) && fixed.isKey(collisions(2))) || (numbOfColl==1 && fixed.isKey(collisions(1)) ))
        deposed=1;
    else
        deposed=-1;
        switch numbOfColl
            case 1
                pos=positions(collisions(1),:)-shiftsColl(1,:)+normr(pos+shiftsColl(1,:)-positions(collisions(1),:))*2*particleRadius;
               
                if(~fixed.isKey(collisions(1)))
                    fixed(collisions(1))=shiftsColl(1,:);
                end
            case 2
                rollOverSpheres(collisions(1),collisions(2),shiftsColl(1,:),shiftsColl(2,:),pos,0);
                 if(~fixed.isKey(collisions(1)))
                fixed(collisions(1))=shiftsColl(1,:);
                 end
                  if(~fixed.isKey(collisions(2)))
                
                fixed(collisions(2))=shiftsColl(2,:);
                  end
            case 3
                
                
            
                if(isObtuseTriangle(positions(collisions(1),:)-shiftsColl(1,:),positions(collisions(2),:)-shiftsColl(2,:),positions(collisions(3),:)-shiftsColl(3,:)))
                    rollOver1=collisions(1);
                    rollOver2=collisions(2);
                    rollShift1=shiftsColl(1,:);
                    rollShift2=shiftsColl(2,:);
                    other=3;
                    side=norm(positions(rollOver1,:)-rollShift1-positions(rollOver2,:)+rollShift2);
                    if(side<norm(positions(collisions(3),:)-shiftsColl(3,:)-positions(rollOver1,:)+rollShift1))
                        rollOver2=collisions(3);
                        rollShift2=shiftsColl(3,:);
                        side=norm(positions(rollOver1,:)-rollShift1+rollShift2-positions(rollOver2,:));
                        other=2;
                    end
                    if(side<norm(positions(collisions(3),:)-shiftsColl(3,:)-positions(collisions(2),:)+shiftsColl(2,:)))
                        rollOver2=collisions(3);
                        rollShift2=shiftsColl(3,:);
                        rollOver1=collisions(2);
                        rollShift1=shiftsColl(2,:);
                        other=1;
                    end
                    rollOverSpheres(rollOver1,rollOver2,rollShift1,rollShift2,positions(collisions(other),:),1);
                    
                    if(fixed.isKey(rollOver1) && fixed.isKey(rollOver2))
                        deposed=0;
                    else
                    if(~fixed.isKey(rollOver1))
    
                    fixed(rollOver1)=rollShift1;
                     end
                      if(~fixed.isKey(rollOver2))
                  
                    fixed(rollOver2)=rollShift2;
                    
                      end
                    end
                      if(fixed.isKey(other))
                        fixed.remove(other);
                      end
                else
                    deposed=0;
                end
            otherwise
                deposed=0;
        end
        if(deposed==-1)
           deposed=resolveCollisions(); 
        end
        
    end

end
        function rollOverSpheres(coll1, coll2,shiftsColl1,shiftsColl2,point,inv)
            
            deltaDir=-normr(positions(coll1,:)-shiftsColl1-positions(coll2,:)+shiftsColl2);
            scal=sqrt(4*particleRadius^2-(norm(positions(coll1,:)-shiftsColl1-positions(coll2,:)+shiftsColl2)/2)^2);
            if(inv)
                pos=(positions(coll1,:)-shiftsColl1+positions(coll2,:)-shiftsColl2)/2-normr(point-positions(coll1,:)+shiftsColl1-deltaDir*((point-positions(coll1,:)+shiftsColl1)*deltaDir'))*scal;
            else
                pos=(positions(coll1,:)-shiftsColl1+positions(coll2,:)-shiftsColl2)/2+normr(point-positions(coll1,:)+shiftsColl1-deltaDir*((point-positions(coll1,:)+shiftsColl1)*deltaDir'))*scal;

            end
            if(~isreal(pos))
                fprintf('asdf');
            end
        end
        
        
    function drawAll()
        [x,y,z]=sphere;
        for i=1:N
            
          
                if(deposed(i))
                    c=[0 1 0];
                else
                    c=[1 0 0];
                end
                fig=surf(x*particleRadius+positions(i,1),y*particleRadius+positions(i,2),z*particleRadius);
                set(fig, 'FaceColor', c);
            end
        
        hold off;
        
    end
%    function c=coverage()
%        totalArea=(xMax-xMin)*(yMax-yMin);
%        overlapArea=0;
%        nIn=0;
%        inside(1,N)=Particle;
%        for j=1:N
%            within(xMin,xMax,particles(j).position(1));
%            within(yMin,yMax,particles(j).position(2));
%            if(within(xMin,xMax,particles(j).position(1)) && within(yMin,yMax,particles(j).position(2)))
%                nIn=nIn+1;
%                inside(nIn)=particles(nIn);
%               % text(inside(nIn).position(1),inside(nIn).position(2),inside(nIn).position(3)+particleRadius*1.5,'in');
%            end
%        end
%        
%        for i=1:nIn
%           for k=1:nIn
%               if (i~=k)
%                 dist=norm([inside(i).position(1);inside(i).position(2)]-[inside(k).position(1);inside(k).position(2)]);
%                 q=dist/2/particleRadius;
%                 if(q<=1)
%                     overlapArea=overlapArea+particleRadius.^2*(acos(q)-q*sqrt(1-q.^2));
%                 end
%               end
%           end
%        end
%        c=(nIn*pi*particleRadius^2-overlapArea)/totalArea;
%        end
%    

   
 function shift=periodicBoundary(position)
        sideX=xMax-xMin;
        sideY=yMax-yMin;
        shift=[0,0];
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
            

            
            

    

function fig=drawSphere(pos)
    try
        delete(virtFigures);
    catch
    end
        [x,y,z]=sphere;
    
            fig=surf(x*particleRadius+pos(1),y*particleRadius+pos(2),z*particleRadius);
            
        

        drawnow;

        
end


function fig=drawPart(i)
        try
        delete(f(i));
        catch
        end
        [x,y,z]=sphere;
    
            fig=surf(x*particleRadius+pos(1),y*particleRadius+pos(2),z*particleRadius);
            
        

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








    function cells=getCells(pos)
             %cell of ith particle
             m=ceil((pos(1)-xMin)./deltaX);
             n=ceil((pos(2)-yMin)./deltaY);
%              m1=m+1;
%              m2=m-1;
%              n1=n+1;
%              n2=n-1;

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
        
        n=0;
        only1=0;
        numbOfCells=0;
        chkDistX=mod(positions(i,1)-xMin+deltaX/2,deltaX)+move(1);
        chkDistY=mod(positions(i,2)-yMin+deltaY/2,deltaY)+move(2);
        if(chkDistX<0 || chkDistX>deltaX || chkDistY<0 || chkDistY>deltaY)
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
        else
            bCoord=1;
            delta=deltaX;
            only1=1;
            step=[0, 0, 0];
        end
        cells=zeros(1+floor(move(bCoord)/delta),2);
        pos=-delta;
        while(abs(pos)<abs(move(bCoord)) || n==0 ||only1)
             currCells=getCells(positions(i,:)+n*step+periodicBoundary(positions(i,:)+n*step));       
             for k=1:4
                 if(~any(find(cells(:,1)==currCells(k,1) & cells(:,2)==currCells(k,2),1)))
                     numbOfCells=numbOfCells+1;
                     cells(numbOfCells,:)=currCells(k,:);
                 end
             end
             n=n+1;
             pos=pos+delta;
             if(only1)
                 only1=0;
             end
        end
        
        
        end
                       
                        
 
 
end