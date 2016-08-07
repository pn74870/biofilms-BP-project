[x1,y1,z1]=sphere;

[x,y,z]=textread('results.txt','%f %f %f');
figure;
axis('equal');
hold on;
for i=1:length(x)
    surf(x1+x(i),y1+y(i),z1+z(i));
end
%hold off;