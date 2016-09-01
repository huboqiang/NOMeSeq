OBJ_DIR = ./
BIN_DIR = ./

SOURCES= main.cpp gzstream.cpp
OBJECTS= main.o   gzstream.o

BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= depth_stat
BIN_OBJECT   = $(patsubst %,$(BIN_DIR)/%,$(PROGRAM))

INCLUDES= -I 
LIBRARY = -L 
LDFLAGS = -lz -lpthread

$(BIN_OBJECT): $(BUILT_OBJECTS)
	@$(CXX) -o $(BIN_OBJECT) $(BUILT_OBJECTS)  $(LDFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(LIBRARY)
$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX)      -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(LIBRARY)
	@echo $(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(LIBRARY)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/main.o  $(OBJ_DIR)/gzstream.o

.PHONY: clean
