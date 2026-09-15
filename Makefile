#
#	Top-level Makefile for native Monte Carlo builds.
#
#	Build the native Monte Carlo executable with:
#
#		make local-build
#
#	This builds the demo against the host toolchain and the UxHw
#	compatibility shim in `submodules/compat`, rather than against the
#	Signaloid cloud compiler.
#

#
#	Include `config.mk` first, so that the variables it sets cannot
#	overwrite the ones this Makefile defines below. `config.mk` is the
#	source list shared with the Signaloid cloud build.
#
include src/config.mk

CC		= gcc

#
#	Use `gnu11` rather than `c11`: the UxHw compatibility shim calls the
#	POSIX `random()` and `srandom()`, which a strict ISO C dialect does
#	not declare. Under `c11` they are implicitly declared as returning
#	`int`, while they in fact return `long`.
#
CFLAGS		= -std=gnu11 -O2 -Isrc/
LDFLAGS		=
LIBS		= -lgsl -lgslcblas -lm

#
#	Platform-specific paths (macOS with MacPorts or Homebrew).
#
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    ifneq ($(wildcard /opt/local/include/gsl),)
        CFLAGS  += -I/opt/local/include
        LDFLAGS += -L/opt/local/lib
    else ifneq ($(wildcard /opt/homebrew/include/gsl),)
        CFLAGS  += -I/opt/homebrew/include
        LDFLAGS += -L/opt/homebrew/lib
    else ifneq ($(wildcard /usr/local/include/gsl),)
        CFLAGS  += -I/usr/local/include
        LDFLAGS += -L/usr/local/lib
    endif
endif

BINARY		= demo-native-mc
SRC_DIR		= src

#
#	`config.mk` omits `uxhw.c`, because the Signaloid cloud compiler
#	provides the UxHw API natively. The native Monte Carlo build has to
#	compile the compatibility shim in explicitly.
#
PROJECT_C	= $(addprefix $(SRC_DIR)/,$(SOURCES)) $(SRC_DIR)/uxhw.c

PROJECT_C_OBJ	= $(PROJECT_C:.c=.o)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(BINARY): $(PROJECT_C_OBJ)
	$(CC) $(LDFLAGS) -o $@ $(PROJECT_C_OBJ) $(LIBS)

.PHONY: local-build clean

local-build: $(BINARY)

clean:
	rm -f $(PROJECT_C_OBJ) $(BINARY)
