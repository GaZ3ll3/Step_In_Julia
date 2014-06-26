export cg

function cg(A::Union(AbstractMatrix, Function),b::Vector; 
	tol::Real=length(b)*eps(),maxIter::Int64=length(b),
	M::Union(AbstractMatrix, Function)=x->x ,x::Vector=[],out::Int=0)

	n = length(b)
	nb = norm(b)
	
	Ap = zeros(eltype(b),n) 
	Af =  isa(A,Function) ? A : x->A*x
	Mf =  isa(M,Function) ? M : x->M\x
	
	if isempty(x)
		x = zeros(eltype(b),n)
		r = copy(b)
	else
		r = b - Af(x)
	end	
	z = Mf(r)
	p = copy(z)	
	
	iter   = 1 
	err    = norm(r)/nb
	flag   = -1
	
	for iter=1:maxIter
		Ap = Af(p)
		gamma = dot(r,z)
		alpha = gamma/dot(p,Ap)
		
		# postive definite required
		if alpha==Inf || alpha<0
			flag = -2; break
		end
		
		BLAS.axpy!(n,alpha,p,1,x,1) # x = alpha*p+x	
		BLAS.axpy!(n,-alpha,Ap,1,r,1) # r -= alpha*Ap 
	
		err = norm(r)/nb

		if err <= tol
			flag = 0; break
		end
		
		z    = Mf(r)
		beta = dot(z,r)/gamma

		# p = z + beta*p
		p = BLAS.scal!(n,beta,p,1)
		p = BLAS.axpy!(n,1.0,z,1,p,1)
	end
	
	if out>=0
		if flag==-1
			println("Does not converge within maximum iterations.")
		elseif flag==-2
			println("Matrix A in cg has to be positive definite.")
		end
	end
	return x,flag, err, iter
end