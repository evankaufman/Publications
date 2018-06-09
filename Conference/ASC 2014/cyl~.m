function [faces vertices] = cyl(r,h1,h2,n,EAngles)
theta = linspace(0,2*pi,n);
x = r*cos(theta);
y = r*sin(theta);
top = transpose(EAngles321(EAngles)*[x',y',h2*ones(n,1)]');
bottom = transpose(EAngles321(EAngles)*[x',y',h1*ones(n,1)]');
vertices = [bottom;top;0 0 h1; 0 0 h2];
faces = [(1:n)',(n+1:2*n)',((n+1:2*n)+1)'];
faces(end,end) = 1;
faces = [faces;(1:n)',(2:n+1)',((n+1:2*n)+1)'];
faces(end,2) = 1;
faces(end,3) = n+1;
faces=[faces; (1:n-1)' (2:n)' (2*n+1)*ones(n-1,1)];
faces=[faces; (n+1:2*n-1)' (n+2:2*n)' (2*n+2)*ones(n-1,1)];

end