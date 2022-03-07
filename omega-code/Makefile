CXX=g++
CXXFLAGS= -O3 -Istonne/include/ -Istonne/external/ -std=c++1y #-ltcmalloc  -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free 
#DEBUGFLAGS=-D DEBUG_MEM_OUTPUT -D DEBUG_MSWITCH_FUNC
#-D DEBUG_MEM_OUTPUT -D DEBUG_MEM_INPUT -D DEBUG_ASWITCH_CONFIG -D DEBUG_ASWITCH_FUNC -D DEBUG_MSWITCH_CONFIG -D DEBUG_MSWITCH_FUNC
BIN_OMEGA=omega
BIN_STONNE=stonne_single
SOURCE=$(wildcard stonne/src/*.cpp)
OBJSDIR=stonne/objs
OBJS_OMEGA=$(patsubst stonne/src/%,$(OBJSDIR)/%,$(patsubst %.cpp,%.o,$(SOURCE)))
OBJS_OMEGA := $(filter-out stonne/objs/main.o, $(OBJS_OMEGA))
OBJS_STONNE=$(patsubst stonne/src/%,$(OBJSDIR)/%,$(patsubst %.cpp,%.o,$(SOURCE)))
OBJS_STONNE := $(filter-out stonne/objs/omega.o, $(OBJS_STONNE))
#OBJS=$(patsubst stonne/src/%,$(OBJSDIR)/%,$(patsubst %.cpp,%.o,$(SOURCE)))
INCLUDES=$(wildcard stonne/include/*.h)
$(warning   OMEGA is $(OBJS_OMEGA))
$(warning   STONNE is $(OBJS_STONNE))

all: $(BIN_OMEGA) $(BIN_STONNE)

$(BIN_OMEGA): $(OBJSDIR) $(OBJS_OMEGA)
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS)  -o $@ $(OBJS_OMEGA)  #-pthread -ltcmalloc


$(BIN_STONNE): $(OBJSDIR) $(OBJS_STONNE)
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS)  -o $@ $(OBJS_STONNE)  #-pthread -ltcmalloc


$(OBJSDIR):
	mkdir -p $@

$(OBJSDIR)/%.o: stonne/src/%.cpp $(INCLUDES)
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) -c $< -o $@  #-ltcmalloc



clean:
	rm -rf $(OBJSDIR) 


