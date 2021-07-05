function [x,xc]=cubetest(file, np, dx, nc, dxc)

  if ( nargin < 5 ) dxc = dx ; endif
  if ( nargin < 4 ) nc = 0 ; endif 
  
  x = (rand(np, 3)-0.5) ;
  x(:,1) = x(:,1)*(dx(2)-dx(1)) + 0.5*(dx(2)+dx(1)) ;
  x(:,2) = x(:,2)*(dx(4)-dx(3)) + 0.5*(dx(4)+dx(3)) ;
  x(:,3) = x(:,3)*(dx(6)-dx(5)) + 0.5*(dx(6)+dx(5)) ;

  r2 = sum(x.^2, 2) ;

  f = exp(-r2) ;
  
  fid = fopen(file, "w") ;

  fprintf(fid, "%d %d %d %d\n", np, 3, 1, 0) ;

  dat = [x f]' ;

  fprintf(fid, "%1.16e %1.16e %1.16e %1.16e\n", dat) ;

  xc = [linspace(dxc(1), dxc(2), nc)' ...
	linspace(dxc(3), dxc(4), nc)' ...
	linspace(dxc(5), dxc(6), nc)'] ;

  r2 = sum(xc.^2, 2) ;

  f = exp(-r2) ;

  fprintf(fid, "%d\n", nc) ;

  dat = [xc f]' ;

  fprintf(fid, "%1.16e %1.16e %1.16e %1.16e\n", dat) ;

  fclose(fid) ;
