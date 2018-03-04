Fisier README pentru tema 3 

Cazul continuu - eps = 0.1 

1.Lagrange
	E =  0.47513
	E =  0.18757
	E =  0.21549
	nod = Inf

2.Newton
	E =  0.47513
	E =  0.18757
	E =  0.21549
	nod = Inf


3.Spline Liniar
	E =  0.48658
	E =  0.12988
	E =  0.025654
	E =  0.0061438	
	nod =  64


4.Spline Cubic Natural
	E =  0.47177
	E =  0.088267
	E =  0.018813
	nod =  32
	

5.Spline Cubic Tensionat
	E =  0.47625
	E =  0.087916
	E =  0.018814
	nod =  32



6.Trigonometric
	E =  0.18013
	E =  0.027732
	E =    1.7354e-04
	nod =  32


Cel mai exact interpolant este Spline Liniar ,
eroarea fiind cea mai apropiata de 0 pentru eps = 0.1 .

Cel mai repede converg Spline Cubic Natural ,
Spline Cubic Tensionat , Trigonometric deoarece
sunt necesare numai 3 noduri .

Cazul discret - eps = 0.1


1.Lagrange
	E =  156.18
	E =  5649.7
	E =    5.9312e+07
	E =    3.5489e+16
	nod = Inf

2.Newton
	E =  156.18
	E =  5649.7
	E =    5.9312e+07
	E =    3.5476e+16
	nod = Inf


3.Spline Liniar
	E =  115.10
	E =  118.10
	E =  122.07
	E =  57.708
	E =  28.546
	E =  11.021
	E = 0
	nod =  300

4.Spline Cubic Natural
	E =    7.9991e+04
	E =    1.3614e+04
	E =  1863.9
	E =  1811.2
	E =  127.15
	E =  19.067
	E = 0
	nod =  300


5.Spline Cubic Tensionat
	E =    9.6077e+04
	E =    1.4662e+04
	E =  1879.9
	E =  1905.5
	nod = Inf



6.Trigonometric
	E =  282.38
	E =  324.33
	E =  257.53
	E =  286.57
	nod = Inf



Am implementat o functie separata pentru calculul erorii

function E = error(xi,pnk,N) 
	S = 0;
	for i=1:N+1
		S = S + abs(f(xi(i)) - pnk(i))^2;
	endfor
	
	h = (2*pi)/(N+1);
	E = sqrt(S*h);

endfunction


Atat scriptul de test.m si test_grafic.m au 
un timp de rulare de aproximativ 10 minute
