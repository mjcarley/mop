/* mop-logging.c
 * 
 * Copyright (C) 2006 Michael Carley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/**
 * @defgroup logging Logging functions
 * @{
 * 
 */

#include <stdio.h>

#include <glib.h>

#include "mop.h"
#include "mop-logging.h"

/* void mop_logging_func(const gchar *log_domain, */
/* 		      GLogLevelFlags log_level, */
/* 		      const gchar *message, */
/* 		      gpointer data[]) ; */
/* const gchar *mop_logging_string(GLogLevelFlags level) ; */

static const gchar *mop_logging_string(GLogLevelFlags level)

{
  const gchar *strings[] = {"RECURSION", 
			    "FATAL",
			    "ERROR",
			    "CRITICAL",
			    "WARNING",
			    "MESSAGE",
			    "INFO",
			    "DEBUG"} ;

  switch (level) {
  default: g_assert_not_reached() ; break ;
  case G_LOG_FLAG_RECURSION: return strings[0] ; break ;
  case G_LOG_FLAG_FATAL: return strings[1] ; break ;
  case G_LOG_LEVEL_ERROR: return strings[2] ; break ;
  case G_LOG_LEVEL_CRITICAL: return strings[3] ; break ;
  case G_LOG_LEVEL_WARNING: return strings[4] ; break ;
  case G_LOG_LEVEL_MESSAGE: return strings[5] ; break ;
  case G_LOG_LEVEL_INFO: return strings[6] ; break ;
  case G_LOG_LEVEL_DEBUG: return strings[7] ; break ;
  }

  return NULL ;
}

static void mop_logging_func(const gchar *log_domain,
		      GLogLevelFlags log_level,
		      const gchar *message,
		      gpointer data[])

{
  FILE *f = (FILE *)data[MOP_LOGGING_DATA_FID] ;
  gchar *p = (gchar *)data[MOP_LOGGING_DATA_PREFIX] ;
  GLogLevelFlags level = *(GLogLevelFlags *)data[MOP_LOGGING_DATA_LEVEL] ;
  gint (*exit_func)(void) = data[MOP_LOGGING_DATA_EXIT_FUNC] ;

  if ( log_level > level ) return ;

  fprintf(f, "%s%s-%s: %s\n", p, 
	  G_LOG_DOMAIN, mop_logging_string(log_level & G_LOG_LEVEL_MASK),
	  message) ;

  if ( log_level <= G_LOG_LEVEL_ERROR ) {
    if ( exit_func != NULL ) exit_func() ;
  }

  return ;
}

/** 
 * Initialize MOP logging
 * 
 * @param f file stream for messages
 * @param p string to prepend to messages
 * @param log_level maximum logging level to handle (see g_log)
 * @param exit_func function to call if exiting on an error
 * 
 * @return 0 on success
 */

gint mop_logging_init(FILE *f, gchar *p, 
		      GLogLevelFlags log_level,
		      gpointer exit_func)

{
  static gpointer data[MOP_LOGGING_DATA_WIDTH] ;
  static GLogLevelFlags level ;

  if ( f != NULL ) 
    data[MOP_LOGGING_DATA_FID] = f ;
  else
    data[MOP_LOGGING_DATA_FID] = stderr ;    
  if ( p != NULL ) 
    data[MOP_LOGGING_DATA_PREFIX] = g_strdup(p) ;
  else
    data[MOP_LOGGING_DATA_PREFIX] = g_strdup("") ;

  level = log_level ;
  data[MOP_LOGGING_DATA_LEVEL] = &level ;    
    
  g_log_set_handler (G_LOG_DOMAIN, 
		     G_LOG_FLAG_RECURSION |
		     G_LOG_FLAG_FATAL |   
		     G_LOG_LEVEL_ERROR |
		     G_LOG_LEVEL_CRITICAL |
		     G_LOG_LEVEL_WARNING |
		     G_LOG_LEVEL_MESSAGE |
		     G_LOG_LEVEL_INFO |
		     G_LOG_LEVEL_DEBUG,
		     (GLogFunc)mop_logging_func, data);

  return 0 ;
}

/**
 * @}
 * 
 */
