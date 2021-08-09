[U,S,V] = svd(A) ;

[m,n] = size(A) ;

r = rank(A) ;

[Q,R,P] = qr([V(1:r,1:r) V(1:r,r+1:n) ;
	    V(r+1:n,1:r) V(r+1:n,r+1:n)]) ;


