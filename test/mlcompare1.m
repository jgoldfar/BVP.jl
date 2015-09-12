function mlcompare1()

function dydx = twoode(x,y)
   dydx = [ y(2); -abs(y(1)) ];
end

function res = twobc(ya,yb)
  res = [ ya(1); yb(1) + 2 ];
end
solinit = bvpinit(linspace(0,4,5),[1 0]);
sol = bvp4c(@twoode,@twobc,solinit);
x = linspace(0,4);
y = deval(sol,x);

h=fopen('mldata1.jl','w');
fprintf(h,'mldata=[\n\t[\n\t\t');
fprintf(h,'%e ', y(1,1:end-1));
fprintf(h,'%e\n\t]\t, \t[\n\t\t', y(1,end));
fprintf(h,'%e ',y(2,1:end-1));
fprintf(h,'%e\n\t]\n]\n', y(2,end));
fclose(h);
end