#ifndef _MOPBLOCK_H_INCLUDED_
#define _MOPBLOCK_H_INCLUDED_

#include <glib.h>

gint block_derivatives(gdouble *x, gint n,
		       gint dim, gint pmax, gint *d, 
		       gdouble *xp) ;
gint block_powers(gdouble *x, gint n, gint dim,
		  gint pmax, gdouble *xp) ;
gdouble block_multipower(gdouble *xp, gint i, gint n, gint dim,
			 gint *p) ;

#endif /*_MOPBLOCK_H_INCLUDED_*/
