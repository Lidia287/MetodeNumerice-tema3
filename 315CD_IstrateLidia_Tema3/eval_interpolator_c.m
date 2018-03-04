function nod = eval_interpolator_c(tip,epss)

N = 1000;
xi = linspace(-pi,pi,N+1); %punctele in care se calculeaza suma din eroare

if tip == 1

	[nod,pnk] = Lagrange(epss);
	 

endif


if tip == 2

	[nod,pnk] = Newton(epss);
		
endif



if tip == 3

	[nod,pnk] = SplineLiniar(epss);

endif


if tip == 4

	[nod,pnk] = SplineNatural(epss);
endif



if tip == 5

	[nod,pnk] = SplineTensionat(epss);
	
endif


if tip == 6

	[nod,pnk] = Trigonometric(epss);
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
%2^k
E;
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

E;

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
	disp("Valoarea nu poate fi calculata in afara intervalului");
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
   % E = 0;
Epr = E;
E = error(xi,pnk,N);

E;
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

E;
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

E;
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

E;

	k = k + 1;
  	endwhile

nod = 2^k;
if E >= Epr
	nod = inf;
endif


	


endfunction
