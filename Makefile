PROJECT = tap

SRCDIR = src
OBJDIR = obj
BINDIR = bin
DEPDIR = .depend
INCLUDEDIR = include

$(shell mkdir -p $(DEPDIR) > /dev/null)

INCLUDEFLAG = -I $(INCLUDEDIR)
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

SOURCES := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(INCLUDEDIR)/*.h)
OBJECTS := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
RM = rm -f

CC = gcc
CFLAGS = -std=c99 -Wall $(INCLUDEFLAG) $(DEPFLAGS)
CFLAGS += -Wextra -Wwrite-strings -Wno-parentheses
CFLAGS += -Wpedantic -Warray-bounds 
CFLAGS += -fmax-errors=10
DEBUGFLAGS = -g -O0
RELEASEFLAGS = -O3
PROFILEFLAGS = -pg $(DEBUGFLAGS)
LINKER = gcc
LFLAGS = -Wall -lm $(INCLUDEFLAG)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td


# ------- all target: build the main project ------

.PHONY: all
all: $(BINDIR)/$(PROJECT)

$(BINDIR)/$(PROJECT): $(OBJECTS)
	$(LINKER) $^ $(LFLAGS) -o $@ 
	
# ---------- release target: extra optimization ----

.PHONY: release
release: CFLAGS += $(RELEASEFLAGS)
release: $(BINDIR)/$(PROJECT)

# ---------- debug target---------------------------

.PHONY: debug
debug: CFLAGS += $(DEBUGFLAGS)
debug: $(BINDIR)/$(PROJECT)

# ---------- profile target-------------------------

.PHONY: debug
profile: CFLAGS += $(PROFILEFLAGS)
profile: $(BINDIR)/$(PROJECT)

# ---------- compile objects

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(DEPDIR)/%.d
	$(CC) $(CFLAGS) -c $< $(INCLUDEFLAG) -o $@
	$(POSTCOMPILE)

$(OBJDIR)/%.o: $(TESTDIR)/%.c
	$(CC) $(CFLAGS_TEST) -c $< $(INCLUDEFLAG_TEST) -o $@
	$(POSTCOMPILE)

# ---------- clean/clear
.PHONY: clean clear
clean clear:
	@$(RM) $(OBJECTS)

.PHONY: remove
remove: clean
	@$(RM) $(BINDIR)/$(PROJECT)

.DELETE_ON_ERROR:

# -------- dependency tracking
$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS))))


