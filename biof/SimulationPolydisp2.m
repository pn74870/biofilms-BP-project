function [ cov] = SimulationPolydisp2( step_size,time,rate,box_size,stick_prob)
%Simulation of polydisperse particle deposition

% %Input parameters:

% 1) step_size is the time step size of physics simulation (tested with values
% from 0.001 to 0.01). Smaller values give better accuracy but the simulation then
% takes more time to finish
%
% 2) time is the total simulated time
%
% 3) rate is the number of generated particles per unit time
%
% 4) box_size is the length of the box side with periodic boundary conditions
%
% 5) stick_prob is the probability that during a collision between a moving particle
% and a deposited one, they stick together
%

%returns an array of the surface coverage evaluated at every unit of time

tic; %start timer





persistent particleRadius;
persistent xMax; %box corners
persistent xMin;
persistent yMin;
persistent yMax;
persistent N; %number of created particles
persistent movingIndices;
%coordinates of box corners
xMax=box_size/2;
xMin=-xMax;
yMax=xMax;
yMin=-yMax;

zMax=1e-5; %height of the box


maxSpheresInPart=3; %largest number of spheres in one particle

%data of particles:
positions=zeros(maxSpheresInPart*time*rate,3); %positions of indvidual spheres
positionsCenter=zeros(time*rate,3); %positions of centers of the rigid bodies



omegas=zeros(3,time*rate);
deposited=zeros(1,time*rate*maxSpheresInPart); %contains information whether a particle has been deposited onto surface
touchingInfo=cell(1,time*rate*maxSpheresInPart); %contains information about particles that are in contact
figures=gobjects(time*rate*maxSpheresInPart); %contains graphical objects of spheres
Nmoving=0;

quaternions=zeros(4,time*rate);
debugOn=0;
%---------------
%physical parameters
particleRadius=1e-6;
T=300; %K
kb=1.38e-23;
d=2*particleRadius; %diameter
eta=8.9e-4;%viscocity
mass=1e-14;
dens=mass/(4/3*pi*particleRadius^3);
gforce=mass*(1-1000/dens)*[0;0;-9.8];%effective force due to gravity and buoyancy
stickProbability=stick_prob; %probability that two particles stick together after they collide



restFramePos=zeros(3,3,3);
restFramePos(1,1,2)=particleRadius;
restFramePos(1,2,2)=-particleRadius;
restFramePos(1,:,3)=[0,-particleRadius,particleRadius];
restFramePos(2,:,3)=[2,-1,-1]*particleRadius/sqrt(3);

momOfInert0(:,:,1)=particleRadius^2*2/5*[1,0,0;0,1,0;0,0,1];
momOfInert0(:,:,2)=particleRadius^2/5*[4,0,0;0,14,0;0,0,14];
momOfInert0(:,:,3)=particleRadius^2/5*[16,0,0;0,16,0;0,0,26];

mobMatrices(:,:,1)=getMobilityMatrix(1);
mobMatrices(:,:,2)=getMobilityMatrix(2);
mobMatrices(:,:,3)=getMobilityMatrix(3);

startEndInd=zeros(time*rate,2);
types=zeros(time*rate,1);
%---------------

pairDir=containers.Map('KeyType','int32','ValueType','any');



virtFigures=gobjects(1,time*rate);


for o=1:time*rate*maxSpheresInPart
    touchingInfo{o}=containers.Map('KeyType','int32','ValueType','any');
end

notMoving=zeros(1,time*rate);

nColumns=floor((xMax-xMin)/2/particleRadius);
nRows=floor((yMax-yMin)/2/particleRadius);
nLayers=floor(zMax/2/particleRadius);
grid=zeros(nColumns,nRows,4*nLayers);
nInGrid=zeros(nColumns,nRows);


inCells=zeros(4,2,time*rate*2); %every particle is in up to 4 cells of the grid
deltaX=(xMax-xMin)/nColumns;
deltaY=(yMax-yMin)/nRows;

tested=java.util.HashSet;


%velocity=[0;0;-50];


N=0; %number of created particles
nSpheres=0;
calcCov=1; %1 to calculate surface coverage as a function of time




figure();
hold on;
plane1=[xMin,xMax,xMax,xMin];
plane2=[yMin,yMin,yMax,yMax];
plane3=[0,0,0,0];

fill3(plane1,plane2,plane3,'r');
alpha 0.1;
axis equal;
axis([2*xMin,2*xMax,2*yMin,2*yMax,0,1.5*zMax]);
xlabel('x');
ylabel('y');

cov=zeros(1,time);

interval=1/rate;
nIntervals=time*rate;

meanMove=zeros(1,time*rate);
sumOfMoves=0;
nMoves=0;
drawOn=1;
mColor=[1 1 0]; %color of the spheres drawn at the end of the simulation
for t=1:nIntervals
    
    
    

    N=N+1;
    
    Nmoving=Nmoving+1;
    movingIndices(Nmoving)=t;
    type=randi(3,1);
    types(N)=type;
    %physics simulation
    wrongPos=1;
    startEndInd(N,:)=[nSpheres+1 , nSpheres+type];
    while(wrongPos)
        %wrongPos=0;
        allCells=zeros(4,2,type);
        positionsCenter(t,:)=[xMin+rand()*(xMax-xMin);yMin+rand()*(yMax-yMin);zMax];
        quaternions(:,t)=randomQuaternion();
        posit=getPositions(t);
        for parts=1:type
            [wrongPos,cells1]=checkColl(posit(parts,:));
            
            if(wrongPos)
                break;
            end
            allCells(:,:,parts)=cells1;
        end
         
        
    end
    for part=1:type
        nSpheres=nSpheres+1;
        positions(nSpheres,:)=posit(part,:);
        inCells(:,:,nSpheres)=allCells(:,:,part);
        for cl=1:4
            addToGrid(inCells(cl,1,nSpheres),inCells(cl,2,nSpheres),nSpheres);      
        end
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
                 if(abs(deltaR(3))<2*particleRadius  && tested.add(p))
                     if(norm(deltaR)<2*particleRadius)
                         wrongPos=1;
                         break;
                     end
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
        for i=1:nSpheres
            [x,y,z]=sphere;

            fig=surf(x*particleRadius+positions(i,1),y*particleRadius+positions(i,2),z*particleRadius+positions(i,3));
            set(fig, 'FaceColor', mColor);
            figures(i)=fig;
        end
        hold off;
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
    function pos=getPositions(i)
        pos=zeros(types(i),3);
        R=quaternionToMatrix(quaternions(:,i));
        for p=1:types(i)
            pos(p,:)=(R*restFramePos(:,p,types(i))+positionsCenter(i,:)')';
            pos(p,:)=pos(p,:)+periodicBoundary(pos(p,:));
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
        for i=1:nSpheres
            if (positions(i,3)<=particleRadius)
                coveredArea=coveredArea+pi*particleRadius^2;
            end
        end
        c=coveredArea/totalArea;
    end
    function ep=leviCivita()
        ep=zeros(3,3,3);
       for i=1:3
           for j=1:3
               for k=1:3
                   if(i~=j && j~=k && i~=k)
                       
                       P=[1,2,3;circshift([1 2 3],[0 1]);circshift([1 2 3],[0 2])];
                       if(any( ismember(P,[i j k],'rows')))
                           ep(i,j,k)=1;
                       else
                           ep(i,j,k)=-1;
                       end
                       
                 
                   end
               end
           end
       end
    end

    function q=omegaToQuaternion(omega,time)
        q=[sin(norm(omega')*time/2);cos(norm(omega')*time/2*omega/norm(omega'))];
    end
    
    function A=getCrossProdMatrix(r)
         A=[0, r(3), -r(2); -r(3), 0, r(1); r(2), -r(1), 0]; 
    end
    function N=getMobilityMatrix(type) %matrix for the equation [u;omega]=N*[Force;Torque]
        %for detailed descriptio see 
        %Improved hydrodynamic interaction in macromolecular bead models
        %B. Carrasco and J. Garc??a de la Torre
        
        %most variables are named the same as in the article
        %Rtt,Rtr... are the resistance matrices for the whole body
        
        I=eye(3);%identity matrix
        utt=zeros(type*3);
        urt=zeros(type*3);
        urr=zeros(type*3);
        
        for i=1:type
            utt(3*i-2:3*i,3*i-2:3*i)=1/6*I;
            urr(3*i-2:3*i,3*i-2:3*i)=1/8/particleRadius^2*I;
            for j=i+1:type
                rij=restFramePos(:,j,type)-restFramePos(:,i,type);
                rijNorm=norm(rij);
                rijDir=rij/rijNorm;
                epsilonRij=getCrossProdMatrix(rijDir); 
                Pij=rij*rij'/(rij'*rij);
                utt(3*i-2:3*i,3*j-2:3*j)=1/8*particleRadius*(1/rijNorm*(I+Pij)+2*particleRadius^2/rijNorm^3*(I-3*Pij)); %modified Oseen tensor
                utt(3*j-2:3*j,3*i-2:3*i)=utt(3*i-2:3*i,3*j-2:3*j);
                
                urt(3*i-2:3*i,3*j-2:3*j)=-1/8*particleRadius/rijNorm^2*epsilonRij;
                urt(3*j-2:3*j,3*i-2:3*i)=-urt(3*i-2:3*i,3*j-2:3*j);
                
                urr(3*i-2:3*i,3*j-2:3*j)=1/16*particleRadius/rijNorm^3*(3*Pij-I);
                urr(3*j-2:3*j,3*i-2:3*i)=urr(3*i-2:3*i,3*j-2:3*j);
            end
        end
        utr=urt';
        
       % zetaMatrix=inv([utt,utr;urt,urr]);
        
        
        zetaTT=inv(utt-utr*inv(urr)*urt);
        zetaRR=inv(urr-urt*inv(utt)*utr);
        zetaTR=-inv(urr)*utr*zetaTT;
        zetaRT=zetaTR';
        
        
        Rtt=zeros(3);
        Rtr=zeros(3);
        Rrr=zeros(3);
       
        for i=1:type
            ri=restFramePos(:,i,type);
            Ai=getCrossProdMatrix(ri);
            for j=1:type
                rj=restFramePos(:,j,type);
                Aj=getCrossProdMatrix(rj);
                Rtt=Rtt+zetaTT(3*i-2:3*i,3*j-2:3*j);
                Rtr=Rtr-zetaTT(3*i-2:3*i,3*j-2:3*j)*Aj+zetaTR(3*i-2:3*i,3*j-2:3*j);
                Rrr=Rrr+zetaRR(3*i-2:3*i,3*j-2:3*j)-zetaRT(3*i-2:3*i,3*j-2:3*j)*Aj+Ai*zetaTR(3*i-2:3*i,3*j-2:3*j)-Ai*zetaTT(3*i-2:3*i,3*j-2:3*j)*Aj;
                
            end
        end
        Rrt=Rtr';
        multiplier=1/pi/eta/particleRadius;
        N=inv([Rtt,Rtr;Rrt,Rrr])*multiplier;
    end
  function N=getMobilityMatrixTest(type,k) %matrix for the equation [u;omega]=N*[Force;Torque]
        %for detailed descriptio see 
        %Improved hydrodynamic interaction in macromolecular bead models
        %B. Carrasco and J. Garc??a de la Torre
        
        %most variables are named the same as in the article
        %Rtt,Rtr... are the resistance matrices for the whole body
        
        I=eye(3);%identity matrix
        utt=zeros(type*3);
        urt=zeros(type*3);
        urr=zeros(type*3);
        
        for i=1:type
            utt(3*i-2:3*i,3*i-2:3*i)=1/6*I;
            urr(3*i-2:3*i,3*i-2:3*i)=1/8/particleRadius^2*I;
            for j=i+1:type
                rij=positions(startEndInd(k,1)+i-1,:)'-positions(startEndInd(k,1)+j-1,:)';
                rij=rij+getShiftDeltaR(rij')';
                rijNorm=norm(rij);
                rijDir=rij/rijNorm;
                epsilonRij=getCrossProdMatrix(rijDir); 
                Pij=rij*rij'/(rij'*rij);
                utt(3*i-2:3*i,3*j-2:3*j)=1/8*particleRadius*(1/rijNorm*(I+Pij)+2*particleRadius^2/rijNorm^3*(I-3*Pij)); %modified Oseen tensor
                utt(3*j-2:3*j,3*i-2:3*i)=utt(3*i-2:3*i,3*j-2:3*j);
                
                urt(3*i-2:3*i,3*j-2:3*j)=-1/8*particleRadius/rijNorm^2*epsilonRij;
                urt(3*j-2:3*j,3*i-2:3*i)=-urt(3*i-2:3*i,3*j-2:3*j);
                
                urr(3*i-2:3*i,3*j-2:3*j)=1/16*particleRadius/rijNorm^3*(3*Pij-I);
                urr(3*j-2:3*j,3*i-2:3*i)=urr(3*i-2:3*i,3*j-2:3*j);
            end
        end
        utr=urt';
        
       % zetaMatrix=inv([utt,utr;urt,urr]);
        
        
        zetaTT=inv(utt-utr*inv(urr)*urt);
        zetaRR=inv(urr-urt*inv(utt)*utr);
        zetaRT=-inv(urr)*urt*zetaTT;
        zetaTR=zetaRT';
        
        
        Rtt=zeros(3);
        Rtr=zeros(3);
        Rrr=zeros(3);
       
        for i=1:type
            ri=positions(startEndInd(k,1)+i-1,:)'-positionsCenter(k,:)';
            Ai=getCrossProdMatrix(ri);
            for j=1:type
                rj=positions(startEndInd(k,1)+j-1,:)'-positionsCenter(k,:)';
                Aj=getCrossProdMatrix(rj);
                Rtt=Rtt+zetaTT(3*i-2:3*i,3*j-2:3*j);
                Rtr=Rtr-zetaTT(3*i-2:3*i,3*j-2:3*j)*Aj+zetaTR(3*i-2:3*i,3*j-2:3*j);
                Rrr=Rrr+zetaRR(3*i-2:3*i,3*j-2:3*j)-zetaRT(3*i-2:3*i,3*j-2:3*j)*Aj+Ai*zetaTR(3*i-2:3*i,3*j-2:3*j)-Ai*zetaTT(3*i-2:3*i,3*j-2:3*j)*Aj;
                
            end
        end
        Rrt=Rtr';
        multiplier=1/pi/eta/particleRadius;
        N=inv([Rtt,Rtr;Rrt,Rrr])*multiplier;
    end
function [ velocity,omega ] = BrownianMotion( i,force,torque )
%generates random linear and angular velocity of rigid body
% the method is described in 
% Brownian Dynamics of Confined Rigid Bodies
% by Steven Delong,Florencio Balboa Usabiaga,and Aleksandar Donev



rotMatrix1=quaternionToMatrix(quaternions(:,i)); %this matrix rotates points
rotMatrix6DT1=[rotMatrix1',zeros(3);zeros(3),rotMatrix1';]; %this matrix transforms 6D mobility tensor
rotMatrix6D1=[rotMatrix1,zeros(3);zeros(3),rotMatrix1;]; %transpose of rotMatrix6D1
W1=randn(6,1);

Nmatrix=rotMatrix6D1*mobMatrices(:,:,types(i))*rotMatrix6DT1;

v=Nmatrix*[force;torque]+sqrt(4*kb*T/step_size)*chol(Nmatrix)*W1;

velocity=v(1:3);
omega=v(4:6);



end


% function [ move,omega ] = BrownianMotion( i,force )
% %generates random displacement
% D=kb*T/3/pi/eta/d/types(i);
% k=sqrt(2*D*step_size);
% 
% move=randn(1,3)*k+D/kb/T*force*step_size;
% W1=randn(3,1);
% W2=randn(3,1);
% rotMatrix1=quaternionToMatrix(quaternions(:,i));
% Mwt1=1/6/types(i)/pi/eta/particleRadius*inv(rotMatrix1*momOfInert0(:,:,types(i))*rotMatrix1');
% N1=chol(Mwt1);
% omega=sqrt(2*kb*T/step_size*2)*N1*W1;
% 
% % rotMatrix2=quaternionToMatrix(multiplyQuaternions(omegaToQuaternion(omega1,step_size/2),quaternions(:,i)));
% % Mwt2=1/6/types(i)/pi/eta/particleRadius*inv(rotMatrix2*momOfInert0(:,:,types(i))*rotMatrix2');
% % omega=sqrt(kb*T/step_size)*Mwt2/N1*(W1+W2);
% % P=leviC*quaternions(2:4,i);
% % psi=1/2*[-quaternions(2:4,i)';quaternions(1,i)*eye(3)-P];
% % 
% % delQuater=psi*(chol(Mwt)*W)-kb*T/4*W'*Mwt*W*quaternions(:,i);
% 
% 
% end

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
                    if(w && deposited(par))
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

    function drawSpheres(p,c)
       
        [x,y,z]=sphere; 
        for i=startEndInd(p,1):startEndInd(p,2)
            
            if(~isempty(figures))
                try
                    delete(figures(i));
                    
                catch
                end
            end
        
        
        if (any(c))
            figures(i)=surf(x*particleRadius+positions(i,1),y*particleRadius+positions(i,2),z*particleRadius+positions(i,3));
            set(figures(i), 'FaceColor', c);
        else
            figures(i)=surf(x*particleRadius+positions(i,1),y*particleRadius+positions(i,2),z*particleRadius+positions(i,3));
            
        end
        
        drawnow;
        end
        
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
            if(deposited(i) && positions(i,3)<z+thickness+particleRadius  && positions(i,3)>z-particleRadius) %within the box
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


   
 function newVel=fixContacts(partInd,vel)
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
               if(contactSpeed(m)<-1e-15)
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
               
               if(contactSpeed(c)<-1e-15 && ~any(C(C==c)))
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

                    rCenter=zeros(types(i),3);
%                     IMatrix=RotMatrixBodyWorld*momOfInert0(:,:,types(i))/(RotMatrixBodyWorld);
%                     sumV=[0;0;0];
%                     L=[0;0;0];
                    nPar=0;
                     [vel,omegas(:,i)]=BrownianMotion(i,gforce*types(i),[0;0;0]);%+elecForce(i));
                   for par=startEndInd(i,1):startEndInd(i,2)
                       
%                         move1=[dx,dy,dz];
                   % move=(vel*step_size)';
%                        vel1=move1'/step_size;
%                        sumV=sumV+vel1;
%                        L=L+cross(positions(par,:)'-positionsCenter(i,:)',vel1);
                       nPar=nPar+1;
                       rCenter(nPar,:)=positions(par,:)-positionsCenter(i,:);
                       rCenter(nPar,:)=rCenter(nPar,:)+getShiftDeltaR(rCenter(nPar,:));
                    end
                    
                
                    
                    
                    
                %    omegas(:,i)=IMatrix\L;
                %    omegas(:,i)=[0;0;10];
                 %   vel=sumV/maxSpheresInPart;
                    
                 %   move=vel'*step_size;
                    move=vel'*step_size;
                    
                    
                    
                    shift=getShift(positionsCenter(i,:),3*particleRadius,[0,0]);
                    
                    if(types(i)>1)
                 
                    nSimult=zeros(types(i),1);
                    simultColl=zeros(types(i),7);
                    
                    coll=0;
                    
                    %TOUCHING ROTATIONS
                    
                    wrongRot=0;
                    nPar=0;
                    for par=startEndInd(i,1):startEndInd(i,2)
                        nPar=nPar+1;
                        instantVel=cross(omegas(:,i),rCenter(nPar,:)');
                        
                        keys=touchingInfo{par}.keys();
                        for c=1:length(keys)
                            dir=positions(keys{c},:)-positions(par,:)-touchingInfo{par}(keys{c});
                            
                            if(dir*instantVel>1e-16)
                                omegas(:,i)=[0;0;0];
                                wrongRot=1;
                                break;
                            end
                        end
                        if(wrongRot)
                            break;
                        end
                    end
                    
                        %------------------------------------------              
                        if(~wrongRot)
                            
                            minT=step_size;
                            omegaUnitVec=omegas(:,i)/norm(omegas(:,i)');
                            quaternion=[cos(norm(omegas(:,i))'*step_size/2);sin(norm(omegas(:,i)')*step_size/2)*omegaUnitVec]*step_size;
                            quaternion=quaternion/norm(quaternion');
                            moveRot=quaternionToMatrix(multiplyQuaternions(quaternion,quaternions(:,i)))*restFramePos(:,1,types(i))-positions(startEndInd(i,1),:)'+positionsCenter(i,:)';
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
                            
                            
                            
                            
                            cells=get16Cells(positionsCenter(i,:));
                            nPar=0;
                            nVector=RotMatrixWorldOmega*[0;0;1];
                            
                            
                           
                           for sph=startEndInd(i,1):startEndInd(i,2)
                               nPar=nPar+1;
                               posInOmegaFrame=RotMatrixWorldOmega*rCenter(nPar,:)';
                               tested.clear();
                            for ce=1:16
                                column=cells(ce,1);
                                row=cells(ce,2);
                               for itr=1:nInGrid(column,row)
                                    j=grid(column,row,itr);
                                    
                                    
                                    
                                    if((j<startEndInd(i,1) || j>startEndInd(i,2)) && abs(positionsCenter(i,3)-positions(j,3))<abs(moveRot(3))+abs(positionsCenter(i,3)-positions(sph,3))+(2+1e-5)*particleRadius  && any(moveRot) && tested.add(j))
                                        
                                        if(~touchingInfo{sph}.isKey(j))
                                            [dt,tshift]=getTimeToCollRot(i,j,posInOmegaFrame,RotMatrixWorldOmega,shift);
                                            if(dt>0)
                                                if( dt<minT-1e-17)
                                                    minT=dt;
                                                    coll=j;
                                                    
                                                    partColl=sph;
                                                    nSimult(nPar)=0;
                                                elseif(dt<minT+1e-17)
                                                    nSimult(nPar)=nSimult(nPar)+1;
                                                    simultColl(nPar,nSimult(nPar))=j;
                                                    
                                                end
                                            end
                                        end
                                   
                                    end
                                end
                            end
                                
                            
                            
                       
                        
                            dFi=getAngleOfCollGround(posInOmegaFrame,positionsCenter(i,3),nVector);
                            if(dFi>0)
                              
                                dt=dFi/norm(omegas(:,i)');
                                if(dt<minT)
                                    minT=dt;
                                    coll=-1;
                                    partColl=sph;
                                end
                                
                            end
                           end
                            
                            if(minT>step_size)
                                coll=0;
                            end
                           
                            if(coll~=0)
                                angle=norm(omegas(:,i)')/2*minT;
                                quaternions(:,i)=multiplyQuaternions([cos(angle);sin(angle)*omegas(:,i)/norm(omegas(:,i)')],quaternions(:,i));
                                %quaternions(:,i)=quaternions(:,i)/norm(quaternions(:,i));
                                positions(startEndInd(i,1):startEndInd(i,2),:)=getPositions(i);
                                
                                if(coll==-1)
                                    %deposed(partColl)=1; %TODO
                                    depose(i);
                                    movingIndices(movingIndex)=[];
                                    movingIndex=movingIndex-1;
                                else
                                    fshift=getShiftDeltaR(positions(partColl,:)-positions(coll,:));
                                    touchingInfo{partColl}(coll)=fshift;
                                    if(~deposited(coll))
                                        touchingInfo{coll}(partColl)=-fshift;
                                    elseif(stick())
                                        depose(i);
                                        movingIndices(movingIndex)=[];
                                        movingIndex=movingIndex-1;
                                    end
                                    for sp=1:types(i)
                                        sph=startEndInd(i,1)+sp-1;
                                        if(nSimult(sp)>0)
                                            for w=1:nSimult(sp)
                                                fshift=getShiftDeltaR(positions(sph,:)-positions(simultColl(sp,w),:));
                                                touchingInfo{sph}(simultColl(sp,w))= fshift;
                                                if(~deposited(simultColl(sp,w)))
                                                    touchingInfo{simultColl(sp,w)}(sp)=-fshift;
                                                end
                                            end
                                        end
                                    end
                                  
                                end
                            else
                                angle=norm(omegas(:,i)')/2*step_size;
                                quaternions(:,i)=multiplyQuaternions([cos(angle);sin(angle)*omegas(:,i)/norm(omegas(:,i)')],quaternions(:,i));
                               % quaternions(:,i)=quaternions(:,i)/norm(quaternions(:,i));
                                positions(startEndInd(i,1):startEndInd(i,2),:)=getPositions(i);
                            end
                           
                            for sph=startEndInd(i,1):startEndInd(i,2)
                                onPositionChanged(sph);
                                if(debugOn)
                                for j=1:nSpheres
                                    if(startEndInd(i,1)>j || startEndInd(i,2)<j)
                                       
                                        p=sph;
                                        if( ~touchingInfo{p}.isKey(j) && (norm(positions(p,:)-positions(j,:))-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+shift)-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+[shift(1),0,0])-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+[0,shift(2),0])-2*particleRadius<0))
                                            fprintf('critical error + %s\n',norm(positions(i,:)-positions(j,:))-2*particleRadius);
                                        end

                                    end
                                end
                                end
                            end
                            
                         
                        end       
        %------------------------------------------
                    end
                deltaTime=step_size;
                numberOfIter=0;
                    while(deltaTime>0 && any(move) && ~deposited(startEndInd(i,1)) && numberOfIter<10)
                        
                        numberOfIter=numberOfIter+1;
                      
                 
                        
                        %------------touching
                            oldVel=vel;
                            vel=fixContacts(startEndInd(i,1):startEndInd(i,2),vel);
%                             
%                             wrongMove=0;
                            for sph=startEndInd(i,1):startEndInd(i,2)
                                keys=touchingInfo{sph}.keys();
                                for c=1:length(keys)
                                    dir1=positions(keys{c},:)-positions(sph,:)-touchingInfo{sph}(keys{c});
                                    if(dir1*vel>1e-14)
                                        wrongMove=1;
                                        break;
                                    end
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
                        
                        
                                                   
                                    if(any(move>1e-15))
                                        
                                        minT=deltaTime;
                                        coll=0; %is there a collision
                                        fshift=[0,0,0];
                                        nSimult=zeros(types(i),1);
                                        simultColl=zeros(types(i),7);
                                        nPar=0;
                                        for sph=startEndInd(i,1):startEndInd(i,2)
                                            nPar=nPar+1;
                                            tested.clear();
                                            
                                            [cellsMove,numbC]=getCellsMoving(sph,move);
                                            shift=getShift(positions(sph,:),2*particleRadius,move);
                                           
                                            for k=1:numbC
                                                column=cellsMove(k,1);
                                                row=cellsMove(k,2);
                                                for itr=1:nInGrid(column,row)
                                                    j=grid(column,row,itr);
                                                    if(startEndInd(i,1)>j || startEndInd(i,2)<j)
                                                        if(abs(positions(sph,3)-positions(j,3))<abs(move(3))+(2+1e-5)*particleRadius  && any(move) && tested.add(j) &&~touchingInfo{sph}.isKey(j))
                                                            [dt,tshift]=getTimeToColl(sph,j,shift,vel,move);
                                                            
                                                            if(dt>0 && dt<minT-1e-17)
                                                                % f(j)=drawSphere(j,[0 1 0]);
                                                                minT=dt;
                                                                partColl=sph;
                                                                coll=j;
                                                                fshift=tshift;
                                                                nSimult(nPar)=0;
                                                            elseif(dt>0 && dt<minT+1e-17)
                                                                
                                                                 nSimult(nPar)=nSimult(nPar)+1;
                                                                 simultColl(nPar,nSimult(nPar))=j;
                                                                
                                                            end
                                                        end
                                                        
                                                        
                                                    end
                                                    %   for j=1:N
                                                    %collision detection
                                                    %test1=touchingInfo{i}.isKey(j);
                                                    
                                                    
                                                    
                                                end
                                                
                                            end
                                           
                                            %-----------------
                                            if(positions(sph,3)+move(3)<particleRadius) %collision with the surface
                                                
                                                
                                                dt=(particleRadius-positions(sph,3))/vel(3);
                                                
                                                if(dt<minT && dt<deltaTime)
                                                    minT=dt;
                                                    coll=-1;
                                                    partColl=sph;
                                                end
                                                
                                                
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
                                
                                for sph=startEndInd(i,1):startEndInd(i,2)
                                    positions(sph,:)=positions(sph,:)+vel'.*minT;
                                    
                                end
                               
                                
                                sumOfMoves=sumOfMoves+vel'*vel.*minT/6;
                                nMoves=nMoves+1;
                                if(coll>0)
                                    %                                 numbOfColl(i)=numbOfColl(i)+1;
                                    %                                 if(numbOfColl(i)>collLimit)
                                    %                                    deposed(i)=1;
                                    %                                    break;
                                    %                                 end
                                    fshift=getShiftDeltaR(positions(partColl,:)-positions(coll,:));
                                    vel=(vel-normr((positions(coll,:)-positions(partColl,:)-fshift))'.*(normr(positions(coll,:)-positions(partColl,:)-fshift)*vel));
                                  %  vel=vel-(normr((positions(coll,:)-positions(partColl,:)-fshift)))'.*(normr(positions(coll,:)-positions(partColl,:)-fshift)*vel);
                                    
                                    touchingInfo{partColl}(coll)=fshift;
                                    if(~deposited(coll))
                                        touchingInfo{coll}(partColl)=-fshift;
                                    elseif(stick())
                                        depose(i);
                                        deltaTime=0;
                                        movingIndices(movingIndex)=[];
                                        movingIndex=movingIndex-1;
                                    end
                                    for sp=1:types(i)
                                        sph=startEndInd(i,1)+sp-1;
                                        if(nSimult(sp)>0)
                                            for w=1:nSimult(sp)
                                                fshift=getShiftDeltaR(positions(sph,:)-positions(simultColl(sp,w),:));
                                                touchingInfo{sph}(simultColl(sp,w))= fshift;
                                                if(~deposited(simultColl(sp,w)))
                                                    touchingInfo{simultColl(sp,w)}(sp)=-fshift;
                                                end
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
                                    depose(i);
                                    deltaTime=0;
                                    movingIndices(movingIndex)=[];
                                    movingIndex=movingIndex-1;
                                end
                                
                                
                                %           f(i)=drawSphere(i,[0 1 1]);
                                %   end
                            else
                               positionsCenter(i,:)=positionsCenter(i,:)+move;
                               positionsCenter(i,:)=positionsCenter(i,:)+periodicBoundary(positionsCenter(i,:));
                                
                                for sph=startEndInd(i,1):startEndInd(i,2)
                                    positions(sph,:)=positions(sph,:)+move;
                                    
                                end
                                
                                
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
                   
                    
                    
                   for sph=startEndInd(i,1):startEndInd(i,2)
                                onPositionChanged(sph);
                                if(debugOn)
                                for j=1:nSpheres
                                    if(startEndInd(i,1)>j || startEndInd(i,2)<j)
                                       
                                        p=sph;
                                        if( ~touchingInfo{p}.isKey(j) && (norm(positions(p,:)-positions(j,:))-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+shift)-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+[shift(1),0,0])-2*particleRadius<0 || norm(positions(p,:)-positions(j,:)+[0,shift(2),0])-2*particleRadius<0))
                                            fprintf('critical error + %s\n',norm(positions(i,:)-positions(j,:))-2*particleRadius);
                                        end

                                    end
                                end
                                end
                            end

                    end
                     
                    end
                movingIndex=movingIndex+1;
                if(drawOn)
                    drawSpheres(i,mColor);
                end
            end
         %         drawFrame();
            time=time-step_size;
        end
        
            end
       
 
    function depose(i)
        deposited(startEndInd(i,1):startEndInd(i,2))=ones(1,types(i));
        Nmoving=Nmoving-1;
    end
            function [dt,fshift]=getTimeToCollRot(i,j,posOmegaFrameI,RotMatrixWorldOmega,shift) %RotMatrixWorldOmega is problem probably
                
                dt=-1;
                fshift=[0,0,0];
                
               if(norm(positionsCenter(i,:)-positions(j,:))<2.001*particleRadius+norm(posOmegaFrameI))
                                            posOmegaFrameJ=RotMatrixWorldOmega*(positions(j,:)-positionsCenter(i,:))'; 
                                            dFi=angleOfColl(posOmegaFrameJ,posOmegaFrameI);
                                            if(dFi>0)
                                            
                                                dt=dFi/norm(omegas(:,i)');
                                            
                                                
                                            end
                                        end
                                        if(any(shift))
                                            
                                            if(norm(positionsCenter(i,:)-positions(j,:)+shift)<2.001*particleRadius+norm(posOmegaFrameI))
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
                                            
                                            if(norm(positionsCenter(i,:)-positions(j,:)+[0,shift(2),0])<2.001*particleRadius+norm(posOmegaFrameI))
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
                                            
                                            if(norm(positionsCenter(i,:)-positions(j,:)+[shift(1),0,0])<2.001*particleRadius+norm(posOmegaFrameI))
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
                                    if(~deposited(keys{q}))
                                        touchingInfo{keys{q}}(i)=touchingInfo{keys{q}}(i)+perShift;
                                    end
                                    
                                end
                                if(norm(positions(keys{q},:)-positions(i,:)-touchingInfo{i}(keys{q}))>(2+1e-5)*particleRadius)
                                    touchingInfo{i}.remove(keys{q});
                                    if(~deposited(keys{q}))
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