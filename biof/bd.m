function [ cov] = bd( step_size,time,rate)
%Simulation of spherical particle deposition

tic; %start timer
v=0;





persistent particleRadius;
persistent xMax;
persistent xMin;
persistent yMin;
persistent yMax;
persistent N;

xMax=0.8e-0;
xMin=-xMax;
yMax=0.8e-0;
yMin=-yMax;

zMax=4.5e-1;
zMin=4.5e-1;



positions=zeros(time*rate,3);
deposed=zeros(1,time*rate);
removed=zeros(1,time*rate);
velocities=zeros(1,time*rate);

particleRadius=1e-1;
virtFigures=gobjects(1,time*rate);
%--------del---------------
% touchingShift=zeros(12,3,time*rate);
% numberOfTouch=zeros(1,time*rate);
% touchingInd=zeros(time*rate,12);
%------------------------
touchingInfo=cell(1,time*rate);
for o=1:time*rate
    touchingInfo{o}=containers.Map('KeyType','int32','ValueType','any');
end

notMoving=zeros(1,time*rate);

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
                         grid(inCells(cel,1,N),inCells(cel,2,N)).add(N);
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
% if (calcCov)
%     figure();
%     plot(cov);
% end
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

    function drawAll()
        [x,y,z]=sphere;
        for i=1:N
            
            if(~removed(i))
                if(deposed(i))
                    c=[0 1 0];
                else
                    c=[1 0 0];
                end
                fig=surf(x*particleRadius+positions(i,1),y*particleRadius+positions(i,2),z*particleRadius+positions(i,3));
                set(fig, 'FaceColor', c);
            end
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


    function physicsStep(time)
        while (time>0)
            for i=1:N
                if(~deposed(i) &&~removed(i))
                    [dx,dy,dz]=BrownianMotion(step_size,T,eta,d,gforce);
                    move=[dx,dy,dz];
                    
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
               keys=touchingInfo{i}.keys();
               len=length(keys);
               if(len>0)
                   
                   dir=normr(positions(keys{1},:)-positions(i,:)-touchingInfo{i}(keys{1}));
                   if(dir*move'>0)
                        newMove=move-dir.*(dir*move');
                   else
                       newMove=move;
                   end
                   
                   coplanar=1;
                    remove(pairDir,pairDir.keys());
                    index=0;
                    for j=1:len
                        dir1=positions(keys{j},:)-positions(i,:)-touchingInfo{i}(keys{j});
                        if(j==len)
                            index=index+1;
                            pairDir(index)=normr(cross(dir1,dir));
                        else
                            for k=j+1:len
                                 index=index+1;
                                 dir2=positions(keys{k},:)-positions(i,:)-touchingInfo{i}(keys{k});
                                 pairDir(index)=normr(cross(dir1,dir2));
                            end
                        end
                    end
                   dirKeys=pairDir.keys();
                   for c=1:len
                           
                       dir1=positions(keys{c},:)-positions(i,:)-touchingInfo{i}(keys{c});
%                        if(~any(dir2))
%                                dir2=normr(cross(dir1,dir));
%                        elseif(dir2*dir1'>1e-16)
%                            coplanar=0;
%                        end
                           
                       if(dir1*newMove'>0)
                       
                          
                           
                           
                            newMove=dir2*(dir2*move');
                           if(c>2 && (dir1*newMove'>1e-16 || ~coplanar)) %try to replace with  dir1*newMove'>1e-16
                               newMove=[0,0,0];
                               break;
                           end
                          
                 
                       end   
                  
                   end
                   
                  
                    move=newMove; 
                 
                    for c=1:len
                        dir1=positions(keys{c},:)-positions(i,:)-touchingInfo{i}(keys{c});
                        if(dir1*move'>1e-16)
                   %         error('wrong direction');
                        end
                    end
                    vel=move'./deltaTime;
                   %------------------------------
               end    
                  if(any(move))
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
                        iter=grid(column,row).iterator();
                        while (iter.hasNext())
                            j=iter.next();
                           
                         %   for j=1:N
                        %collision detection
                        test1=touchingInfo{i}.isKey(j);
                        test2=tested.contains(j);
                        if(i~=j &&~removed(j) && ~test1 && ~test2)
                            tested.add(j);
   
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
                                    elseif(dt>0 && dt<deltaTime&& dt<minT)
                                        % f(j)=drawSphere(j,[0 1 0]);
                                        minT=dt;
                                        coll=j;
                                        fshift=[0,0,0];
                                        nSimult=0;
                                    elseif(dt>0 && dt<deltaTime&& dt==minT)
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
                                    if(dt>0 && dt<deltaTime)
                                        tShift=shift;
                                        
                                        primChck=1;
                                    elseif(dt<=0 && norm(positions(i,:)-positions(j,:)+shift)<=2*particleRadius &&(positions(j,:)-positions(i,:)-shift)*vel>0 )
                                        fprintf('this should not have happened');
                                    end
                                end 
       
                                if(~primChck && shift(1)~=0 )
                                    delta_r=positions(j,:)-positions(i,:)-[shift(1),0,0];
                                    disc=discriminant(delta_r,vel);
                                    if(disc>=0)
                                    dt=(delta_r*vel-sqrt(disc))/(vel'*vel);
                                    if(dt>0 && dt<deltaTime)
                                        tShift=[shift(1),0,0];
                                        primChck=1;
                                    elseif(dt<=0 && norm(positions(i,:)-positions(j,:)+[shift(1),0,0])<=2*particleRadius &&(positions(j,:)-positions(i,:)-[shift(1),0,0])*vel>0 )
                                       fprintf('this should not have happened');
                                    end
                                    end
                                end   
                                 
                                
                                if(~primChck && shift(2)~=0)
                                    delta_r=positions(j,:)-positions(i,:)-[0,shift(2),0];
                                    disc=discriminant(delta_r,vel);
                                  if(disc>=0)
                                    dt=(delta_r*vel-sqrt(disc))/(vel'*vel);
                                if(dt>0 && dt<deltaTime )
                                    tShift=[0,shift(2),0];
                                    primChck=1;
                                elseif(dt<=0 && norm(positions(i,:)-positions(j,:)+[0,shift(2),0])<=2*particleRadius &&(positions(j,:)-positions(i,:)-[0,shift(2),0])*vel>0 )
                                        fprintf('this should not have happened');
                                end 
                                  end     
                                    
                                end
                                if (primChck)
                                    %                                     delta_r=positions(j,:)-positions(i,:)-tShift;
                                    %                                     discr=(delta_r*vel)^2-vel'*vel*(delta_r*delta_r'-(2*particleRadius)^2);
                                    
                                    
                                    if(dt<=0 && norm(positions(i,:)-positions(j,:)+tShift)<=2*particleRadius &&(positions(j,:)-positions(i,:)-tShift)*vel>0 )
                                        fprintf('this should not have happened');
                                        %TODO
                                        
                                    elseif ( dt<minT)
                                        % f(j)=drawSphere(j,[0 1 0]);
                                        minT=dt;
                                        coll=j;
                                        fshift=tShift;
                                        nSimult=0;
                                    elseif(dt==minT)
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
                                 
                                 %---------del------------
%                                  numberOfTouch(i)=numberOfTouch(i)+1;
%                                  touchingShift(numberOfTouch(i),:,i)=fshift;
%                                  touchingInd(i,numberOfTouch(i))=coll;
%                                      
%                                      if(~deposed(coll) && ~any(touchingInd(coll,:)==i) )
%                                          numberOfTouch(coll)=numberOfTouch(coll)+1;
%                                          touchingShift(numberOfTouch(coll),:,coll)=-fshift;
%                                          touchingInd(coll,numberOfTouch(coll))=i;
%                                      end
                                 %--------del-------------------
                            end
                            
                           
                 %           f(i)=drawSphere(i,[0 1 1]);
                     %   end
                    else
                        positions(i,:)=positions(i,:)+move;
                        deltaTime=0;
                
                    end
                    perShift=periodicBoundary(positions(i,:));
                    positions(i,:)=positions(i,:)+perShift;
  
               f(i)=drawSphere(i,[0 0 0]); %draw each frame
%                     for j=1:N
%                         if(i~=j && ~removed(j) && (norm(positions(i,:)-positions(j,:))-2*particleRadius<-particleRadius/100 || norm(positions(i,:)-positions(j,:)+shift)-2*particleRadius<-particleRadius/100 || norm(positions(i,:)-positions(j,:)+[shift(1),0,0])-2*particleRadius<-particleRadius/100 || norm(positions(i,:)-positions(j,:)+[0,shift(2),0])-2*particleRadius<-particleRadius/100))
%                             error('critical error + %s\n',norm(positions(i,:)-positions(j,:))-2*particleRadius);
%                         end
%                     end
%                     
                     keys=touchingInfo{i}.keys();
                     
                         
                     for q=1:length(keys)
                         if(any(perShift))
                            
                                 touchingInfo{i}(keys{q})=touchingInfo{i}(keys{q})-perShift;
                                 if(~deposed(keys{q}) && ~removed(keys{q}))
                                     touchingInfo{keys{q}}(i)=touchingInfo{keys{q}}(i)+perShift;
                                 end
                             
                         end
%                         if(norm(positions(keys{q},:)-positions(i,:)-touchingInfo{i}(keys{q}))>(2+1e-10)*particleRadius && coll~=keys{q}) 
%                            touchingInfo{i}.remove(keys{q});
%                            keyTest=sprintf('%d removed %d',i,keys{q});
%                            removedCollisions(keyTest)=keys{q};
%                             if(~deposed(keys{q}))
%                                 touchingInfo{keys{q}}.remove(i);
%                             end
%                         end 
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
                        grid(inCells(c,1,i),inCells(c,2,i)).remove(i);
                    end
                    end
                    
                    inCells(:,:,i)=getCells(positions(i,:));
                    for c=1:4
                         grid(inCells(c,1,i),inCells(c,2,i)).add(i);
                    end
              %-----------------space partition        

                    
                  else
                      
                      removed(i)=1;
                      keys=touchingInfo{i}.keys();
                      for w=1:length(keys)
                          if(~deposed(i))
                          touchingInfo{keys{w}}.remove(i);
                          end
                      end
                    if(any(inCells(1,:,i)))
                    for c=1:4
                        grid(inCells(c,1,i),inCells(c,2,i)).remove(i);
                    end
                    end
                    break;
                      
                  end  
                 end

                end
            end
            
            time=time-step_size;
        end
      
    end

    function d=discriminant(delta_r,vel)
            d=(delta_r*vel)^2-vel'*vel*(delta_r*delta_r'-(2*particleRadius)^2);
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
                       
                        
 
    function overlap=checkCol(i)
        overlap=0;
        for k=1:N
            %collision detection
            if(norm(positions(k,:)-positions(i,:))<2*particleRadius && i~=k) %overlapping spheres
                overlap=1;
                break;
            end
        end
    end
end