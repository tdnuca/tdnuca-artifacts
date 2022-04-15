#ifndef __INTHIST_EXTRA_H__
#define __INTHIST_EXTRA_H__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _NANOS_H_
#pragma message "Compiling sequential version"
#ifdef MANUAL_SOCKET
#undef MANUAL_SOCKET
#endif
#define MANUAL_SOCKET 0
#endif

#ifndef MANUAL_SOCKET
#define MANUAL_SOCKET 1
#endif

#if !MANUAL_SOCKET
#define count_sockets(x) ;
#define current_socket(x) ;

#else
#define count_sockets(x) nanos_get_num_sockets((x));
#define current_socket(x) nanos_current_socket((x));
#endif

#ifndef USE_OMPSS
#define USE_OMPSS 1
#endif

#ifndef SEQUENCE_MODE
#define SEQUENCE_MODE USE_OMPSS
#endif

#ifndef INTHIST_CHECK
#define INTHIST_CHECK 0
#endif

#ifndef USE_RANDOMIMAGE
#define USE_RANDOMIMAGE 1
#endif


#endif // __INTHIST_EXTRA_H__
