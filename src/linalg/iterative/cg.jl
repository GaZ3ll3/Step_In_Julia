export cg

function cg(A::Union(AbstractMatrix, Function),b::Vector; 
	tol::Real=1e-6,maxIter::Int64=1000,
	M::Union(AbstractMatrix, Function)=x->x ,x0::Vector=[],out::Int64=0)

	n = length(b)
	nb = norm(b)
	
	Ap = zeros(eltype(b),n) 
	Af =  isa(A,Function) ? A : x->A_mul_B!(Ap, A, x) # save space
	Mf =  isa(M,Function) ? M : x->M\x
	
	if isempty(x0)
		x0 = zeros(eltype(b),n)
		r = copy(b)
	else
		r = b - Af(x0)
	end	
	z = Mf(r)
	p = copy(z)	
			
	resvec = zeros(maxIter)

	iter   = 1 
	flag   = 1
	
	for iter=1:maxIter
		Ap = Af(p)
		gamma = dot(r,z)
		alpha = gamma/dot(p,Ap)
		
		# if A is ill-conditioned
		if alpha==Inf || alpha<0
			flag = 2; break
		end
		
		BLAS.axpy!(n,alpha,p,1,x0,1) # x = alpha*p+x	
		BLAS.axpy!(n,-alpha,Ap,1,r,1) # r -= alpha*Ap 

		resvec[iter] = norm(r)/nb
		# converge
		if  resvec[iter] <= tol
			flag = 0; break
		end
		
		z    = Mf(r)
		beta = dot(z,r)/gamma

		# p = z + beta*p
		p = BLAS.scal!(n,beta,p,1)
		p = BLAS.axpy!(n,1.0,z,1,p,1)
	end
	# control output
	if out>=0
		if flag == 1
			println("Does not converge within maximum iterations.")
		elseif flag == 2
			println("Matrix A in cg has to be positive definite.")
		end
	end
	return x0,flag, resvec[iter],  resvec, iter
end