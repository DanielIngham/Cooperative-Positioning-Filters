PROJECT := filter
PROJECT_DIR =$(CURDIR)

# Directories
BUILD_DIR :=$(PROJECT_DIR)/build
INCLUDE_DIR :=$(PROJECT_DIR)/include
SRC_DIR :=$(PROJECT_DIR)/src

# Library
LIB_DIR :=$(PROJECT_DIR)/lib
LIBRARIES :=$(LIB_DIR)/DataHandler/lib/libdata_handler.a

# Compiler
CXX := g++

# Flags
WFLAGS := -Wall -Wextra -Werror -Wshadow 
MFLAGS := -ffloat-store -fno-fast-math
CFLAGS := $(WFLAGS) $(MFLAGS) -I$(INCLUDE_DIR) -I$(LIB_DIR)

# Files
TARGET := $(BUILD_DIR)/$(PROJECT)
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))

# Static C++ Code Analyser
CPPCHECK := cppcheck

.PHONY: all clean run cppcheck

# Linking
$(TARGET): $(OBJECTS)
	# Make the DataHandler Library
	$(MAKE) -C $(LIB_DIR)/DataHandler
	# Linking
	@mkdir -p $(dir $@)
	$(CXX) $^  $(LIBRARIES) -o $@ 

# Compiling 
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	# Compile the project
	@mkdir -p $(dir $@)
	$(CXX) -g $(CFLAGS) -c $^ -o $@ 

all: $(TARGET)
	
clean:
	$(MAKE) -C $(LIB_DIR)/DataHandler clean
	rm -rf $(BUILD_DIR)


run: $(TARGET)
	# Set an environment variable called PROJECT_DIR 
	PROJECT_DIR=$(PROJECT_DIR) $(TARGET)

# Static C++ Code Analyser
cppcheck: 
	$(CPPCHECK) --quiet --enable=all --suppress=missingIncludeSystem --error-exitcode=1 -I $(INCLUDE_DIR) $(SRC_DIR)
