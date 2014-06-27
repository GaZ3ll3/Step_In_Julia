export gmres

function gmres(A::Union(AbstractMatrix, Function), b::Vector; restart::Int64=20,
	tol::Real=1e-6,M::Union(AbstractMatrix,Function)=x->x, 
	maxIter::Int64 = length(b), x::Vector=[], out::Int64=0)

	Af =  isa(A,Function) ? A : x->A*x
	Mf =  isa(M,Function) ? M : x->M\x
	
	# initialization
	n  = length(b)

	if isempty(x)
		x = zeros(n)
		r = Mf(b)
	else
		r = Mf(b-Af(x))
	end
	
	bnrm2 = norm(b)
	if bnrm2 == 0.0; bnrm2 = 1.0; end
	
	err = norm( r ) / bnrm2
	if err < tol; return x, err; end
	
	V     = zeros(n,restart + 1)
	H     = zeros(restart + 1,restart)  
	cs    = zeros(restart)  
	sn    = zeros(restart)
	s     = zeros(restart + 1) # allocate once
	w     = zeros(n)
	
	if iseltype(b,Complex) || iseltype(r,Complex) || iseltype(A,Complex)
		s  = complex(s)
		b  = complex(b)
		r  = complex(r)
		x  = complex(x)
		V  = complex(V)
		H  = complex(H)
		cs = complex(cs)
		sn = complex(sn)
	end

	resvec = zeros((restart + 1)*maxIter)
	
	iter = 0
	flag = -1
	cnt  = 1

	@inbounds for iter = 1:maxIter
		s[1]      = norm( r );  
		V[:,1]    = r / s[1];
		
		@inbounds for i = 1:restart
			w = Af(V[:,i])
			# A_mul_B!(w, A, V[:,i]) # save space if A is matrix
			w = Mf(w)
			
			@inbounds for k = 1:i # basis using Gram-Schmidt
				H[k,i] = dot(w,V[:,k])
				BLAS.axpy!(n,-H[k,i],V[:,k],1,w,1)
			end
			H[i+1,i] = norm( w )
			V[:,i+1] = w / H[i+1,i]
			
			@inbounds for k = 1:i-1 # apply Givens rotation
				temp     =  cs[k]*H[k,i] + sn[k]*H[k+1,i]
				H[k+1,i] = -sn[k]*H[k,i] + cs[k]*H[k+1,i]
				H[k,i]   = temp
			end
			
			# Approximate residual norm
			cs[i],sn[i] = rotmat( H[i,i], H[i+1,i] )
			s[i+1] = -sn[i]*s[i]
			s[i]   = cs[i]*s[i]
			H[i,i] = cs[i]*H[i,i] + sn[i]*H[i+1,i]
			H[i+1,i] = 0.0
			err  = abs(s[i+1]) / bnrm2
	
			resvec[cnt] = err
			
			if err <= tol
				y  = H[1:i,1:i] \ s[1:i]
				x += V[:,1:i]*y
				flag = 0; break
			end
			cnt = cnt+1
		end
		if  err <= tol
			flag = 0
			break
		end
		y  = H[1:restart,1:restart]\s[1:restart]
		x += V[:,1:restart]*y
		
		r = b - Af(x)
		r = Mf(r)
		
		s[restart+1] = norm(r)
		resvec[cnt] = abs(s[restart+1]) / bnrm2
		
		if err <= tol
		   flag = 0; break
		end
	end
	
	if flag==-1
		println(@sprintf("gmres iterated maxIter (=%d) times without achieving the desired tolerance.",maxIter))
	elseif out>=1
		println(@sprintf("gmres achieved desired tolerance at iteration %d. Residual norm is %1.2e.",iter,resvec[cnt]))
	end
	return x,flag,resvec[cnt],iter,resvec[1:cnt]
end

function rotmat( a, b )
# c,s = rotmat(a,b)
#  Givens rotation matrix parameters for a and b.

	if b==0.0
		c = 1.0
		s = 0.0
	elseif abs(b) > abs(a)
		temp = a / b
		s    = 1.0/sqrt( 1.0 + temp^2 )
		c    = temp * s
	else
		temp = b/a
		c    = 1.0/sqrt( 1.0 + temp^2 )
		s    = temp * c
	end
	return c, s
end