function f=sourceAnalytical(X)

x = X(:,1); y = X(:,2);

%f=[ones(length(x),1);zeros(length(y),1)];

f = [2*y.*sin(x)+cos(x);-cos(x).*(y.^2)+2*cos(x)];
