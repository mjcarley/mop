function [x,f,w,y,g]=ldtfile(file)

fid = fopen(file, 'r') ;

dat = fscanf(fid, '%d', 4) ;

np = dat(1) ;
nc = dat(2) ;
nf = dat(3) ;
nw = dat(4) ;

dat = fscanf(fid, '%f', (nc+nw+nf)*np) ;
dat = reshape(dat, nc+nw+nf, np)' ;

x = dat(:,1:nc) ;
w = dat(:,nc+1:nc+nw) ;
f = dat(:,nc+nw+1:end) ;

fclose(fid) ;
