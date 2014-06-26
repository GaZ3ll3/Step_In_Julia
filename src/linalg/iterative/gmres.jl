export gmres

function gmres(A::Union(AbstractMatrix, Function), b::Vector; 
	tol::Real=length(b)*eps(),M::Union(AbstractMatrix,Function)=x->x, 
	maxIter::Int64 = length(b), out::Int64=0)


	n = length(b)
	nb = norm(b)