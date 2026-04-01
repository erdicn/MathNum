CC      = gcc
TARGET  = exe.exe
INCDIR  = ../Include
PROJDIR = $(CURDIR)
BUILDDIR = build

INCLUDES = -I. -I$(INCDIR)
CFLAGS   = -O3 -Wall -Wextra -Wshadow -Wcast-align -Wunused-const-variable=2 \
           -Wfloat-equal -Wundef -Wcast-qual -Wunreachable-code \
           $(INCLUDES)
LDFLAGS  = -lm -lgsl -lgslcblas 

SRCS    = $(wildcard *.c) $(shell find $(INCDIR) -name '*.c') #$(wildcard $(INCDIR)/*.c)

INC_SUBDIRS = $(shell find $(INCDIR) -type d)

# This converts "main.c" to "build/main.o" and "../includes/util.c" to "build/util.o"
# We use 'notdir' to keep the build directory flat and clean
OBJS    = $(addprefix $(BUILDDIR)/, $(notdir $(SRCS:.c=.o)))

# VPATH tells Make where to look for source files if they aren't in the current dir
vpath %.c $(INC_SUBDIRS) . #$(INCDIR) .

par: LDFLAGS += -fopenmp
par: CFLAGS += -fopenmp
par: $(BUILDDIR) $(TARGET)

all: $(BUILDDIR) $(TARGET)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILDDIR)/%.o: %.c
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

debug: CFLAGS = -g -O0 -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer
debug: LDFLAGS += -fsanitize=address
debug: clean all

.PHONY: clean test debug

run: all
	./$(TARGET)

clean:
	rm -rf $(BUILDDIR) $(TARGET)

