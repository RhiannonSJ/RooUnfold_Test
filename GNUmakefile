#===============================================================================
# File and Version Information:
#      $Id: GNUmakefile,v 1.20 2010-07-17 01:00:53 adye Exp $
#
# Description:
#      Makefile for the RooUnfold package
#
# Instructions:
#      o Make sure the ROOTSYS environment variable is set and points
#        to your ROOT release, and that $ROOTSYS/bin is in your PATH.
#
#      o run 'make <target>'
#        - Default target makes shared library (libRooUnfold.so) for use
#          in ROOT.
#        - Add NOROOFIT=1 to build test programs (RooUnfoldTest*)
#          without RooFit (this is default if RooFit is not available).
#        - Add SHARED=1 to link test executables with shared library
#          (libRooUnfold.so). Otherwise links with static library
#          (libRooUnfold.a).
#        - Add ROOTBUILD=debug for debug version.
#        - Add VERBOSE=1 to show main commands as they are executed
#
# Build targets:
#      shlib   - make libRooUnfold.so (default target)
#      include - make dependency files (*.d)
#      lib     - make libRooUnfold.a
#      bin     - make lib and example programs
#      commands- show commands to make each type of target
#      html    - make documentation in htmldoc subdirectory
#      cleanbin- delete test binaries and objects
#      clean   - delete all intermediate and final build objects
#
# Author: Tim Adye <T.J.Adye@rl.ac.uk>
#
#===============================================================================

# === ROOT setup ===============================================================
-include $(ROOTSYS)/test/Makefile.arch
ifeq ($(ROOTCONFIG),)
ROOTCONFIG    = $(ROOTSYS)/bin/root-config
endif
ifeq ($(ARCH),)
# === This section is just in case ROOT's test/Makefile.arch didn't work =======
out := $(shell echo "$(ROOTSYS)/test/Makefile.arch not found - trying a basic Linux config" >&2)
ARCH          =   $(shell $(ROOTCONFIG) --arch)
ROOTLIBS      =   $(shell $(ROOTCONFIG) --libs)
ROOTINCLUDES  = -I$(shell $(ROOTCONFIG) --incdir)
CXXFLAGS      =   $(shell $(ROOTCONFIG) --cflags)
CXX           = g++
CXXFLAGS     += -Wall -fPIC
LD            = g++
LDFLAGS       =
SOFLAGS       = -shared
ObjSuf        = o
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"
ifneq ($(findstring debug,$(ROOTBUILD)),)
CXXFLAGS     += -g
LDFLAGS      += -g
endif
else
ROOTINCLUDES  = -I$(shell $(ROOTCONFIG) --incdir)
endif
ifeq ($(VERBOSE),1)
_             =
else
_             = @
endif

# === RooUnfold directories and options ========================================

ifneq ($(findstring g++,$(CXX)),)
MFLAGS        = -MM
endif
SRCDIR        = $(CURDIR)/src/
WORKDIR       = $(CURDIR)/tmp/$(ARCH)/
LIBDIR        = $(CURDIR)/
SHLIBDIR      = $(CURDIR)/
EXEDIR        = $(CURDIR)/
EXESRC        = $(CURDIR)/examples/
INCLUDES      = -I$(SRCDIR)
HTMLDOC       = htmldoc

# === Internal configuration ===================================================

PACKAGE       = RooUnfold
OBJDIR        = $(WORKDIR)obj/
DEPDIR        = $(WORKDIR)dep/

ifeq ($(NOROOFIT),)
ifneq ($(shell $(ROOTCONFIG) --has-roofit),yes)
out := $(shell echo "This version of ROOT does not support RooFit. We will build the test programs without it." >&2)
NOROOFIT      = 1
endif
endif

ifneq ($(NOROOFIT),)
CPPFLAGS     += -DNOROOFIT
else
ifneq ($(wildcard $(ROOTSYS)/lib/libRooFitCore.$(DllSuf)),)
ROOFITLIBS   += -lRooFit -lRooFitCore -lThread -lMinuit -lHtml
else
ROOFITLIBS   += -lRooFit -lMinuit -lHtml
endif
endif

MAIN          = $(notdir $(wildcard $(EXESRC)*.cxx))
MAINEXE       = $(addprefix $(EXEDIR),$(patsubst %.cxx,%$(ExeSuf),$(MAIN)))
MAINDEP       = $(addprefix $(DEPDIR),$(patsubst %.cxx,%.d,$(MAIN)))
ROOTSYS      ?= ERROR_RootSysIsNotDefined
HLIST         = $(filter-out $(SRCDIR)$(PACKAGE)_LinkDef.h,$(wildcard $(SRCDIR)*.h)) $(SRCDIR)$(PACKAGE)_LinkDef.h
CINTFILE      = $(WORKDIR)$(PACKAGE)Cint.cxx
CINTOBJ       = $(OBJDIR)$(PACKAGE)Cint.o
LIBFILE       = $(LIBDIR)lib$(PACKAGE).a
SHLIBFILE     = $(SHLIBDIR)lib$(PACKAGE).$(DllSuf)

ifneq ($(SHARED),)
LIBS          = -L$(SHLIBDIR)
LINKLIB       = $(SHLIBFILE)
LINKLIBOPT    = -l$(PACKAGE)
else
LIBS          = -L$(LIBDIR)
LINKLIB       = $(LIBFILE)
LINKLIBOPT    = -Wl,-static -l$(PACKAGE) -Wl,-Bdynamic
endif

# List of all object files to build
OLIST         = $(addprefix $(OBJDIR),$(patsubst %.cxx,%.o,$(notdir $(wildcard $(SRCDIR)*.cxx))))

ifeq ($(MFLAGS),)

# Can't make dependency files, so make every compilation dependent on all headers.
HDEP          = $(HLIST)

else

# List of all dependency file to make
DLIST         = $(addprefix $(DEPDIR),$(patsubst %.cxx,%.d,$(notdir $(wildcard $(SRCDIR)*.cxx $(EXESRC)*.cxx))))

endif

ifeq ($(NOROOFIT),)
ifeq ($(DLIST),)
# Since we can't check the dependencies, include RooFit on all binaries
ROOFITCLIENTS = $(patsubst %.cxx,$(OBJDIR)%.o,$(MAIN))
else
# List of programs that use RooFit. Should only be those in $(EXESRC)
ROOFITCLIENTS = $(patsubst $(DEPDIR)%.d,$(OBJDIR)%.o,$(shell fgrep -l '/RooAbsArg.h ' $(DLIST) 2>/dev/null))
endif
endif

MAINCINTOBJ = $(addprefix $(OBJDIR),$(patsubst %.d,%Cint.o,$(notdir $(shell grep -l '^ *$(EXESRC).*\.h ' $(MAINDEP) 2>/dev/null))))

# === Implicit rules ===========================================================

# Implicit rule making all dependency Makefiles included at the end of this makefile
$(DEPDIR)%.d : $(SRCDIR)%.cxx
	@echo "Making $@"
	@mkdir -p $(DEPDIR)
	@rm -f $@
	$(_)set -e; \
	 $(CXX) $(MFLAGS) $(CPPFLAGS) $(INCLUDES) $(ROOTINCLUDES) $< \
	 | sed 's,\($(notdir $*)\.o\) *:,$(OBJDIR)\1 $@ :,g' > $@; \
	 [ -s $@ ] || rm -f $@

$(DEPDIR)%.d : $(EXESRC)%.cxx
	@echo "Making $@"
	@mkdir -p $(DEPDIR)
	@rm -f $@
	$(_)set -e; \
	 $(CXX) $(MFLAGS) $(CPPFLAGS) $(INCLUDES) $(ROOTINCLUDES) $< \
	 | sed 's,\($(notdir $*)\.o\) *:,$(OBJDIR)\1 $@ :,g' > $@; \
	 [ -s $@ ] || rm -f $@

# Implicit rule to compile all classes
$(OBJDIR)%.o : $(SRCDIR)%.cxx $(HDEP)
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	$(_)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $(OBJDIR)$(notdir $@) $(INCLUDES)

# Implicit rule to compile main program
$(OBJDIR)%.o : $(EXESRC)%.cxx $(HDEP)
	@echo "Compiling main program $<"
	@mkdir -p $(OBJDIR)
	$(_)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $(OBJDIR)$(notdir $@) $(INCLUDES)

# === Explicit rules ===========================================================

default : shlib

# Rule to make ROOTCINT output file
$(CINTOBJ) : $(HLIST)
	@mkdir -p $(WORKDIR)
	@mkdir -p $(OBJDIR)
	@echo "Running rootcint for $(SRCDIR)$(PACKAGE)_LinkDef.h"
	$(_)cd $(SRC) ; $(ROOTSYS)/bin/rootcint -f $(CINTFILE) -c -p $(INCLUDES) $(HLIST)
	@echo "Compiling $(CINTFILE)"
	$(_)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(CINTFILE) -o $(CINTOBJ) $(INCLUDES)

# Rule to make test ROOTCINT output file
$(MAINCINTOBJ) : $(wildcard $(EXESRC)*.h)
	@mkdir -p $(WORKDIR)
	@mkdir -p $(OBJDIR)
	@echo "Running rootcint for $(addprefix $(EXESRC),$(patsubst %Cint.o,%.cxx,$(notdir $@)))"
	$(_)cd $(SRC) ; $(ROOTSYS)/bin/rootcint -f $(addprefix $(WORKDIR),$(patsubst %.o,%.cxx,$(notdir $@))) -c -p $(INCLUDES) $(shell sed -n 's=^ *\($(EXESRC).*\.h\) .*=\1=p' $(addprefix $(DEPDIR),$(patsubst %Cint.o,%.d,$(notdir $@))))
	@echo "Compiling $(addprefix $(WORKDIR),$(patsubst %.o,%.cxx,$(notdir $@)))"
	$(_)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(addprefix $(WORKDIR),$(patsubst %.o,%.cxx,$(notdir $@))) -o $@ $(INCLUDES) -DNOINLINE

# Rule to combine objects into a library
$(LIBFILE) : $(OLIST) $(CINTOBJ)
	@echo "Making $(LIBFILE)"
	@mkdir -p $(LIBDIR)
	@rm -f $(LIBFILE)
	$(_)ar q $(LIBFILE) $(OLIST) $(CINTOBJ)
	@ranlib $(LIBFILE)

# Rule to combine objects into a shared library
$(SHLIBFILE) : $(OLIST) $(CINTOBJ)
	@echo "Making $(SHLIBFILE)"
	@mkdir -p $(SHLIBDIR)
	@rm -f $(SHLIBFILE)
	$(_)$(LD) $(SOFLAGS) $(LDFLAGS) $(OLIST) $(CINTOBJ) $(OutPutOpt)$(SHLIBFILE) $(ROOTLIBS)

$(MAINEXE) : $(EXEDIR)%$(ExeSuf) : $(OBJDIR)%.o $(LINKLIB) $(filter $(OBJDIR)%Cint.o,$(MAINCINTOBJ))
	@echo "Making executable $@"
	@mkdir -p $(EXEDIR)
	$(_)$(LD) $(LDFLAGS) $< $(OutPutOpt)$@ $(LIBS) $(filter $(OBJDIR)$(patsubst %$(ExeSuf),%Cint.o,$(notdir $@)),$(MAINCINTOBJ)) $(LINKLIBOPT) $(ROOTLIBS) $(if $(findstring $<,$(ROOFITCLIENTS)),$(ROOFITLIBS))

# Useful build targets
include: $(DLIST)
lib: $(LIBFILE)
shlib: $(SHLIBFILE)
bin: shlib $(MAINEXE)

commands :
	@echo "Make $(DEPDIR)%.d:	$(CXX) $(MFLAGS) $(CPPFLAGS) $(INCLUDES) $(ROOTINCLUDES) $(SRCDIR)%.cxx | sed 's,\(%\.o\) *:,$(OBJDIR)\1 $(DEPDIR)%.d :,g' > $(DEPDIR)%.d"
	@echo
	@echo "Compile $(SRCDIR)%.cxx:	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(SRCDIR)%.cxx -o $(OBJDIR)%.o $(INCLUDES)"
	@echo
	@echo "Make $(SHLIBFILE):	$(LD) $(SOFLAGS) $(LDFLAGS) *.o $(OutPutOpt)$(SHLIBFILE) $(ROOTLIBS)"
	@echo
	@echo "Make executable $(EXEDIR)RooUnfoldTest$(ExeSuf):	$(LD) $(LDFLAGS) $(OBJDIR)RooUnfoldTest.o $(OutPutOpt)$(EXEDIR)RooUnfoldTest$(ExeSuf) $(LIBS) $(LINKLIBOPT) $(ROOTLIBS) $(ROOFITLIBS)"

clean : cleanbin
	rm -f $(DLIST)
	rm -f $(CINTFILE) $(basename $(CINTFILE)).h
	rm -f $(OLIST) $(CINTOBJ)
	rm -f $(LIBFILE)
	rm -f $(SHLIBFILE)
	rm -f $(MAINCINTOBJ) $(addprefix $(WORKDIR),$(patsubst %.o,%.cxx,$(notdir $(MAINCINTOBJ))) $(patsubst %.o,%.h,$(notdir $(MAINCINTOBJ))))

cleanbin :
	rm -f $(addprefix $(OBJDIR),$(patsubst %.cxx,%.o,$(MAIN)))
	rm -f $(MAINEXE)

$(HTMLDOC)/index.html : $(SHLIBFILE)
	@echo "Making HTML documentation in $(HTMLDOC)"
	@( echo 'gSystem->Load("libRooUnfold");'; \
	   echo 'THtml h; h.SetOutputDir("$(HTMLDOC)");'; \
	   echo 'h.MakeAll();';\
	   echo '.q') \
	 | root -l -b

html : $(HTMLDOC)/index.html

.PHONY : include shlib lib bin default clean cleanbin html

ifneq ($(DLIST),)
-include $(DLIST)
endif
