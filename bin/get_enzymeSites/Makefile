UTILITIES_DIR = /data/Analysis/huboqiang/software/bedtools-2.17.0/src/utils/
OBJ_DIR = ./
BIN_DIR = ./

SOURCES= main.cpp Fasta.cpp split.cpp
OBJECTS= main.o   Fasta.o   split.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= get_enzymeSite
BIN_OBJECT   = $(patsubst %,$(BIN_DIR)/%,$(PROGRAM))

FASTA = Fasta
INCLUDES= -I $(patsubst %,$(UTILITIES_DIR)/%,$(FASTA))
LIBRARY = -L $(patsubst %,$(UTILITIES_DIR)/%,$(FASTA))
LDFLAGS =

$(BIN_OBJECT): $(BUILT_OBJECTS)
	@$(CXX) -o $(BIN_OBJECT) $(BUILT_OBJECTS)  $(LDFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(LIBRARY)
$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX)      -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(LIBRARY)
	@echo $(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(LIBRARY)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/main.o  $(OBJ_DIR)/Fasta.o   $(OBJ_DIR)/split.o.o

.PHONY: clean