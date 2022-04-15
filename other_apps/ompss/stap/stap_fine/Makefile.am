bin_PROGRAMS = stap-fine

stap_fine_SOURCES = stap_bench_trans.c
stap_fine_LDADD = $(NANOX_LDLIBS)
stap_fine_LINK = $(MCC) --ompss $(AM_CFLAGS) $(CFLAGS) -o $@ $(NANOX_LDFLAGS)
AM_CFLAGS += -Wno-unused-result
CPPFLAGS += $(NANOX_CPPFLAGS)

include $(top_srcdir)/mk/rules.ompss.mk
