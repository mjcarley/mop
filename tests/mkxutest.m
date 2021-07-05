## generate input data for Example 4.3 of Xu, 2004

x = [-1 -1; 0 -1; 1 -1; -1 0; 0 0; 1 0; -1 1; -1 2] ;

[np, dim] = size(x) ;

fid = fopen("xu43test.dat", "w") ;

fprintf(fid, "%d %d %d %d\n", np, dim, 0, 0) ;

fprintf(fid, "%f %f\n", x') ;

fclose(fid) ;
