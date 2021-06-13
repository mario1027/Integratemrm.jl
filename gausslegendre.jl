function gausslegendre(f::Function,n::Int64,a,b)
	w = zeros(n)
	x = zeros(n)
	m  = i = j = t = t1 = pp = p1 = p2 = p3 = 0
	eps = 3E-14   # ajuste apropiado
	m= floor(Int,((n + 1)/2))
	quadra=0
	for i in 1:m
		t = cos(pi*(Float64(i) - 0.25)/(Float64(n) + 0.5) )
		t1=1
		while (abs(t - t1) ) >= eps
			p1 = 1
			p2 = 0
			for j in 1:n
				p3 = p2
				p2 = p1 
                p1 = ((2.0*Float64(j)-1.0)*t*p2 - (Float64(j)-1.0)*p3)/(Float64(j))
			end
			pp = n*(t*p1 - p2)/(t*t - 1.0) 
            t1 = t
			t = t1-p1/pp
		end
		x[i] = - t
		x[n - i+1] = t 
        w[i] = 2.0/( (1.0 - t*t)*pp*pp) 
        w[n - i+1] = w[i]
	end
	for i in 1:n
		x[i] = x[i]*(b - a)/2. + (b + a)/2. 
        w[i] = w[i]*(b - a)/2. 
	end
	return sum(f(x[i])*w[i] for i in 1:n)
end
