# Directories
BUILD_DIR := build
INCLUDE_DIR := include
SRC_DIR := src

# Library
LIBRARY := data_handler
LIB_DIR := lib/DataHandler
LIB_BUILD := $(LIB_DIR)/lib
LIB_INCLUDE := $(LIB_DIR)/include

# Compiler
CXX := g++

# Flags
WFLAGS := -Wall -Wextra -Werror -Wshadow 
MFLAGS := -ffloat-store -fno-fast-math
CFLAGS := $(WFLAGS) $(MFLAGS) -I$(INCLUDE_DIR) -I$(LIB_INCLUDE)

# Files
PROJECT := filter
TARGET := $(BUILD_DIR)/$(PROJECT)
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))

# Static C++ Code Analyser
CPPCHECK := cppcheck

.PHONEY: all clean run cppcheck

# Linking
$(TARGET): $(OBJECTS)
	@mkdir -p $(dir $@)
	$(CXX) $^ -L$(LIB_BUILD) -l$(LIBRARY) -o $@ 

# Compiling 
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	# Make the DataHandler Library
	$(MAKE) -C $(LIB_DIR)
	# Compile the project
	@mkdir -p $(dir $@)
	$(CXX) -g -I $(INCLUDE_DIR) $(CFLAGS) -c $^ -o $@ 

all: $(TARGET)
	
clean:
	$(MAKE) -C $(LIB_DIR) clean
	rm -rf $(BUILD_DIR)

run: $(TARGET)
	$(TARGET)

# Static C++ Code Analyser
cppcheck: 
	$(CPPCHECK) --quiet --enable=all --suppress=missingIncludeSystem --error-exitcode=1 -I $(INCLUDE_DIR) $(SRC_DIR)
