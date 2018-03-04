function matrix = test()

matrix = zeros(2,6);

[nod,pnk] = Lagrange(0.1);
matrix(1,1) = nod;

[nod,pnk] = Newton(0.1);
matrix(1,2) = nod;


[nod,pnk] = SplineLiniar(0.1);
matrix(1,3) = nod;

[nod,pnk] = SplineNatural(0.1);
matrix(1,4) = nod;


[nod,pnk] = SplineTensionat(0.1);
matrix(1,5) = nod;


[nod,pnk] = Trigonometric(0.1);
matrix(1,6) = nod;


[Nk,p] = eval_interpolator_d(1,0.1);
matrix(2,1) = Nk;

[Nk,p] = eval_interpolator_d(2,0.1);
matrix(2,2) = Nk;

[Nk,p] = eval_interpolator_d(3,0.1);
matrix(2,3) = Nk;

[Nk,p] = eval_interpolator_d(4,0.1);
matrix(2,4) = Nk;

[Nk,p] = eval_interpolator_d(5,0.1);
matrix(2,5) = Nk;

[Nk,p] = eval_interpolator_d(6,0.1);
matrix(2,6) = Nk;


matrix 

endfunction


function [Nk,p] = eval_interpolator_d(tip,epss)

x = linspace(1,300,300);

k = 2;
E = 0;
Eprev = 1000;
file = fopen("sunspot.dat","r");
fseek(file,0,SEEK_SET);

for i=1:300
	informatie = fscanf(file,"%s",1);
	puncte(i) = fscanf(file,"%lf",1);
endfor
fclose(file);
nrIter = 1;
while norm(E-Eprev) >= epss && (E<Eprev || nrIter<= 5 ) && k<=9
if k < 9
	interval = floor(linspace(1,300,2^k+1)); 
else
	interval = linspace(1,300,300);
endif

if k < 9
for i= 1:2^k+1
	y(i) = puncte(interval(i));
endfor
else
	for i= 1:300
	y(i) = puncte(interval(i));
endfor
endif



if tip == 1
	for i=1:300
		p(i) = Lagrange_d(interval, y , x(i) );
	endfor
endif

if tip == 2
	for i=1:300
		p(i) = Newton_d(interval, y , x(i) );
	endfor
endif


if tip == 3
	for i=1:300
		p(i) = SplineLiniar_d(interval, y , x(i) );
	endfor
endif

if tip == 4
	for i=1:300
		p(i) = SplineC2natural_d(interval, y , x(i) );
	endfor
endif

if tip == 5
	for i=1:300
		p(i) = SplineC2tensionat_d(interval, y , x(i) );
	endfor
endif

if tip == 6
	for i=1:300
		p(i) = Trigonometric_d(interval, y , x(i) );
	endfor
endif


if nrIter >= 2
	Eprev = E;
endif
E = 0;

if nrIter >= 2
	for i=1:300
		E = E + (abs(puncte(i) - p(i)))^2;
	endfor
endif


E = sqrt(2*pi/301*E);

k = k+1;

nrIter++;

endwhile


k--;
if k < 9
	Nk = 2^k;
else
	Nk = 300;
endif



if E >= Eprev
	Nk = inf;
endif




endfunction



function y =  f(x)
	y = exp(3*cos(x))/(2*pi*besseli(0,3));
endfunction


function E = error(xi,pnk,N) 
	S = 0;
	for i=1:N+1
		S = S + abs(f(xi(i)) - pnk(i))^2;
	endfor
	
	h = (2*pi)/(N+1);
	E = sqrt(S*h);

endfunction

function y =  fder1(x)
	y = (exp(3*cos(x))*(-sin(x)))/(2*pi*besseli(0,3));
endfunction






function [x] = Thomas(a,b,c,d)

n = length(b);

c0(1) = c(1) / b(1);

for i=2:n-1
	c0(i) = c(i) / (b(i) - c0(i-1) * a(i-1));
endfor

d0(1) = d(1)/b(1);

for i=2:n
	d0(i) = (d(i) - d0(i-1) * a(i-1)) / (b(i) - c0(i-1) * a(i-1));
endfor

x(n) = d0(n);

for i=n-1:-1:1
	x(i) = d0(i) - c0(i) * x(i+1);
endfor



endfunction


function [nod,pnk] = Lagrange(epss)


	N = 1000;
	xi = linspace(-pi,pi,N+1);
	k = 2;
	
  	
	Epr = 102;
	 
  E = 0;
  	Eprev = 1000;

	 while norm(E - Eprev) >= epss && E < Eprev


	Nk = power(2,k);
	xk = linspace(-pi,pi,Nk); 
	for i=1:N+1	
		pnk(i)=0;
		for j=1:Nk
			prod=1;
			for l=1:Nk
				if l~=j
					prod=prod*(xi(i)-xk(l))/(xk(j)-xk(l));
				endif	
			endfor
			pnk(i)=pnk(i)+f(xk(j))*prod;
		endfor
	endfor

S = 0;
	for w=1:N+1
		S = S + abs(f(xi(w)) - pnk(w))^2;
	endfor


if k > 2
      Eprev = E;
    endif

	Epr = E;
	h = (2*pi)/(N+1);
	E = sqrt(S*h);

	k = k + 1;
  	endwhile

kk = k;
k = 2^k;
nod = k;

if E >= Epr
	nod = inf;
endif


 for i = 1 : 1001
    a(i) = f(xi(i));
  endfor

 



endfunction


function [nod,pnk] = Newton(epss)

N = 1000;
	xi = linspace(-pi,pi,N+1);
k = 2;
  	

        E = 0;
  	Eprev = 1000;

	 while norm(E - Eprev) >= epss && E < Eprev

	Nk = power(2,k);
	xk = linspace(-pi,pi,Nk); 

	for g=1:N+1

		for w=1:Nk
	d(w) = f(xk(w));
		endfor


for i=1:Nk
	for w=Nk:-1:i+1
		d(w) = (d(w)-d(w-1))/(xk(w)-xk(w-i));
	endfor
endfor

for i=1:Nk
	yk(i)=d(i);
endfor


S = yk(1);
P = 1;

for i=2:Nk
	P = P*(xi(g)-xk(i-1));
	S = S+P*yk(i);
endfor
pnk(g) = S;
	endfor


 

if k > 2
      Eprev = E;
    endif

Epr = E;
E = error(xi,pnk,N);


	k = k + 1;
  	endwhile

nod = 2^k;
if E >= Epr
	nod = inf;
endif



endfunction


function [nod,pnk] = SplineLiniar(epss)

N = 1000;
	xi = linspace(-pi,pi,N+1);
	k = 2;
	E = 0;
  	Eprev = 1000;

	 while norm(E - Eprev) >= epss && E < Eprev
	Nk = power(2,k);
	xk = linspace(-pi,pi,Nk); %punctele pentru nodul Nk

	for g=1:N+1
	if xi(g) < xk(1) || xi(g) > xk(Nk)
	disp("Valorea nu poate fi calculata in afara intervalului");
	else
	i=1;
	while xk(i) < xi(g)
		i++;
	endwhile

	if xk(i) == xi(g)
		rez = f(xk(i));

	else 	
		rez = f(xk(i-1)) + (xi(g) - xk(i-1)) / (xk(i) - xk(i-1)) * (f(xk(i)) - f(xk(i-1)));

	endif

endif

pnk(g) = rez;	
	endfor

    if k > 2
      Eprev = E;
    endif
   
Epr = E;
E = error(xi,pnk,N);

	k = k + 1;
  	endwhile

kk = k;
k = 2^k;
nod = k;
if E >= Epr
	nod = inf;
endif




endfunction



function [nod,pnk] = SplineNatural(epss)

N = 1000;
	xi = linspace(-pi,pi,N+1);
k = 2;
	E = 0;
  	Eprev = 1000;

	 while norm(E - Eprev) >= epss && E < Eprev
	Nk = power(2,k);
	xk = linspace(-pi,pi,Nk); %punctele pentru nodul Nk

	for g=1:N+1

n = length(xk);
a = zeros(1,n);
for i=1:Nk
	y(i) = f(xk(i));
endfor
a = y;
h = zeros(1,n-1);
for i = 1:n-1
	h(i) = xk(i+1) - xk(i);
endfor

bb = zeros(1,n);
dd = zeros(1,n);
aa = zeros(1,n-1);
cc = zeros(1,n-1);
bb(1) = 1;
bb(n) = 1;
d(1) = 0;
d(n) = 0;
for i=2:n-1
	bb(i) = 2*(h(i-1)+h(i));
	dd(i) = 3*(a(i+1)-a(i))/h(i) - 3*(a(i)-a(i-1))/h(i-1);	

endfor

aa = h(1:n-2);
cc = h(2:n-1);
aa(n-1) = 0;
cc = [0 cc];
c = Thomas(aa',bb',cc',dd);

for i=1:n-1
	d(i) = (c(i+1) - c(i))/3*h(i);
	b(i) = (a(i+1) - a(i))/h(i) - h(i)*(2*c(i)+c(i+1))/3;

endfor


for i=1:n-1
	if xi(g) >= xk(i) && xi(g) <= xk(i+1)
		rez = a(i) + b(i)*(xi(g) - xk(i)) + c(i) * (xi(g) - xk(i))^2 + d(i) * (xi(g) - xk(i))^3;
	endif	
endfor	


pnk(g) = rez;	
	endfor

    if k > 2
      Eprev = E;
    endif
  
Epr = E;
E = error(xi,pnk,N);

	k = k + 1;
  	endwhile

kk = k;
k = 2^k;
nod = k;
if E >= Epr
	nod = inf;
endif

	 



endfunction


function [nod,pnk] = SplineTensionat(epss)

N = 1000;
	xi = linspace(-pi,pi,N+1);
k = 2;
	E = 0;
  	Eprev = 1000;

	 while norm(E - Eprev) >= epss && E < Eprev
	Nk = power(2,k);
	xk = linspace(-pi,pi,Nk); %punctele pentru nodul Nk

	for g=1:N+1

n = length(xk);
a = zeros(1,n);
for i=1:Nk
	y(i) = f(xk(i));
endfor
a = y;
h = zeros(1,n-1);
for i = 1:n-1
	h(i) = xk(i+1) - xk(i);
endfor

bb = zeros(1,n);
dd = zeros(1,n);
aa = zeros(1,n-1);
cc = zeros(1,n-1);
bb(1) = 2*h(1);
bb(n) = 2*h(n-1);
d(1) = (f(xk(2))-f(xk(1)))/(xk(2)-xk(1));
d(n) = (f(xk(n))-f(xk(n-1)))/(xk(n)-xk(n-1));
for i=2:n-1
	bb(i) = 2*(h(i-1)+h(i));
	dd(i) = 3*(a(i+1)-a(i))/h(i) - 3*(a(i)-a(i-1))/h(i-1);	

endfor

aa = h;
cc = h;
c = Thomas(aa',bb',cc',dd);

for i=1:n-1
	d(i) = (c(i+1) - c(i))/3*h(i);
	b(i) = (a(i+1) - a(i))/h(i) - h(i)*(2*c(i)+c(i+1))/3;

endfor


for i=1:n-1
	if xi(g) >= xk(i) && xi(g) <= xk(i+1)
		rez = a(i) + b(i)*(xi(g) - xk(i)) + c(i) * (xi(g) - xk(i))^2 + d(i) * (xi(g) - xk(i))^3;
	endif	
endfor	


pnk(g) = rez;	
	endfor

    if k > 2
      Eprev = E;
    endif
  
Epr = E;
E = error(xi,pnk,N);

	k = k + 1;
  	endwhile

kk = k;
k = 2^k;
nod = k;
if E >= Epr
	nod = inf;
endif

	

endfunction




function [nod,pnk] = Trigonometric(epss)

N = 1000;
	xi = linspace(-pi,pi,N+1);
k = 2;
  	
        E = 0;
  	Eprev = 1000;

	 while norm(E - Eprev) >= epss && E < Eprev



	 
	Nk = power(2,k);
	xk = linspace(-pi,pi,Nk+1); 

	for g=1:N+1
	m = Nk/2;
	a=zeros(1,m);
	b=zeros(1,m);
	for j=1:2*m
		for w=1:m+1
			suma=0;
			for i=1:2*m
				suma=suma+f(xk(i))*cos((w-1)*xk(i));
			endfor
			a(w)=1/m*suma;
		endfor

		for w=2:m
			sum2=0;
			for i=1:2*m
				sum2=sum2+f(xk(i))*sin((w-1)*xk(i));
			endfor
			b(w)=1/m*sum2;
		endfor
	endfor

	sum3=0;
	for w=2:m
		sum3=sum3+(a(w)*cos((w-1)*xi(g))+b(w)*sin((w-1)*xi(g)));
	endfor

	S=(a(1)+a(m)*cos(m*xi(g)))/2+sum3;

	pnk(g) = S;
	endfor


if k > 2
      Eprev = E;
    endif

Epr = E;
E = error(xi,pnk,N);


	k = k + 1;
  	endwhile

nod = 2^k;
if E >= Epr
	nod = inf;
endif


	


endfunction


function P=Lagrange_d(xi,yi,x)

suma=0;

n = length(xi);
n = n-1;

for k=1:n+1
	prod=1;
	for i=1:n+1
		if i~=k
			prod=prod*(x-xi(i))/(xi(k)-xi(i));
		endif	
	endfor

		
	suma=suma+yi(k)*prod;
endfor

P=suma;

endfunction


function F=Newton_d(xi,yi,x)

n = length(xi);
n = n-1;

for k=1:n+1
	d(k) = yi(k);
endfor


for i=1:n
	for k=n+1:-1:i+1
		d(k) = (d(k)-d(k-1))/(xi(k)-xi(k-i));
	endfor
endfor

for i=1:n+1
	yi(i)=d(i);
endfor


S = yi(1);
P = 1;

for i=2:n+1
	P = P*(x-xi(i-1));
	S = S+P*yi(i);
endfor

P;
S;

F=S;

endfunction


function [rez] = SplineLiniar_d(x,y,val)

n = length(x);

if val < x(1) || val > x(n)
	disp("Valorea nu poate fi calculata in afara intervalului");
else 
	i=1;
	while x(i) < val
		i++;
	endwhile

	if x(i) == val
		rez = y(i);
		return;
	endif

	rez = y(i-1) + (val - x(i-1)) / (x(i) - x(i-1)) * (y(i) - y(i-1));
endif


endfunction


function [x] = Thomas(a,b,c,d)

n = length(b);

c0(1) = c(1) / b(1);

for i=2:n-1
	c0(i) = c(i) / (b(i) - c0(i-1) * a(i-1));
endfor

d0(1) = d(1)/b(1);

for i=2:n
	d0(i) = (d(i) - d0(i-1) * a(i-1)) / (b(i) - c0(i-1) * a(i-1));
endfor

x(n) = d0(n);

for i=n-1:-1:1
	x(i) = d0(i) - c0(i) * x(i+1);
endfor



endfunction



function yi = SplineC2natural_d(x,y,xi)

n = length(x);
a = zeros(1,n);
a = y;
h = zeros(1,n-1);
for i = 1:n-1
	h(i) = x(i+1) - x(i);
endfor

bb = zeros(1,n);
dd = zeros(1,n);
aa = zeros(1,n-1);
cc = zeros(1,n-1);
bb(1) = 1;
bb(n) = 1;
d(1) = 0;
d(n) = 0;
for i=2:n-1
	bb(i) = 2*(h(i-1)+h(i));
	dd(i) = 3*(a(i+1)-a(i))/h(i) - 3*(a(i)-a(i-1))/h(i-1);	

endfor

aa = h(1:n-2);
cc = h(2:n-1);
aa(n-1) = 0;
cc = [0 cc];
c = Thomas(aa',bb',cc',dd);

for i=1:n-1
	d(i) = (c(i+1) - c(i))/3*h(i);
	b(i) = (a(i+1) - a(i))/h(i) - h(i)*(2*c(i)+c(i+1))/3;
endfor


for i=1:n-1
	if xi >= x(i) && xi <= x(i+1)
		yi = a(i) + b(i)*(xi - x(i)) + c(i) * (xi - x(i))^2 + d(i) * (xi - x(i))^3;
	endif	
endfor




endfunction


function yi = SplineC2tensionat_d(x,y,xi)

n = length(x);
a = zeros(1,n);
a = y;
h = zeros(1,n-1);
for i = 1:n-1
	h(i) = x(i+1) - x(i);
endfor

bb = zeros(1,n);
dd = zeros(1,n);
aa = zeros(1,n-1);
cc = zeros(1,n-1);
bb(1) = 2*h(1);
bb(n) = 2*h(n-1);
d(1) = (y(2)-y(1))/(x(2)-x(1));
d(n) = (y(n)-y(n-1))/(x(n)-x(n-1));
for i=2:n-1
	bb(i) = 2*(h(i-1)+h(i));
	dd(i) = 3*(a(i+1)-a(i))/h(i) - 3*(a(i)-a(i-1))/h(i-1);	

endfor

aa = h;
cc = h;
c = Thomas(aa',bb',cc',dd);

for i=1:n-1
	d(i) = (c(i+1) - c(i))/3*h(i);
	b(i) = (a(i+1) - a(i))/h(i) - h(i)*(2*c(i)+c(i+1))/3;

endfor


for i=1:n-1
	if xi >= x(i) && xi <= x(i+1)
		rez = a(i) + b(i)*(xi - x(i)) + c(i) * (xi - x(i))^2 + d(i) * (xi - x(i))^3;
	endif	
endfor	

yi = rez;


endfunction


function y =  f(x)
	y = exp(3*cos(x))/(2*pi*besseli(0,3));
endfunction

function y = Trigonometric_d(xj,yj,x)

n = length(xj);
m = (n+1)/2;
a = zeros(1,m+1);
b = zeros(1,m);

for k=1:m+1
	suma = 0;
	for j=1:2*m-1
	
	suma += yj(j)*cos((k-1)*xj(j));

	endfor	

	a(k) = 1/m*suma;
	
endfor

b(1) = 0;
for k=2:m

	suma=0;
	for j=1:2*m-1
		suma += yj(j)*sin((k-1)*xj(j));
	endfor
		
	b(k)= 1/m*suma;
endfor


suma = 0;
for k=1:m
	suma += a(k)*cos((k-1)*x) + b(k)*sin((k-1)*x);
endfor
y = (a(1)+a(m+1)*cos((m+1)*x))/2+suma;


endfunction

