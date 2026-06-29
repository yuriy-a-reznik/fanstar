# Makefile for FANSTAR - Faster An* Lattice Quantization
# 
# Usage:
#   make              - build the test executable
#   make clean        - remove build artifacts
#   make rebuild      - clean and build
#   make help         - show this help message
#
# Variables:
#   CC         - C compiler (default: gcc)
#   CFLAGS     - compiler flags (default: -O2 -Wall -Wextra)
#   DEBUG      - set to 1 for debug build with -g -O0

# ============================================================================
# Configuration
# ============================================================================

# Directories
SRC_DIR     := ../src
INC_DIR     := ../include
BUILD_DIR   := .
BIN_DIR     := $(BUILD_DIR)
OBJ_DIR     := $(BUILD_DIR)/obj

# Source files
SOURCES     := $(SRC_DIR)/fanstar_test.c \
              $(SRC_DIR)/fanstar_CS.c \
              $(SRC_DIR)/fanstar_VC.c \
              $(SRC_DIR)/fanstar_MCQ.c \
              $(SRC_DIR)/fanstar_MCSQ.c \
              $(SRC_DIR)/fanstar_MCQ_opt.c \
              $(SRC_DIR)/fanstar_MCSQ_opt.c \
              $(SRC_DIR)/fanstar_faster.c

# Object files
OBJECTS     := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SOURCES))

# Executable
TARGET      := $(BIN_DIR)/fanstar_test

# Compiler and flags
CC          ?= gcc
CFLAGS_BASE := -Wall -Wextra -std=c99
CFLAGS_OPT  := -O2 -march=native
CFLAGS_DBG  := -g -O0 -DDEBUG

# Determine optimization flags based on DEBUG variable
ifeq ($(DEBUG),1)
    CFLAGS := $(CFLAGS_BASE) $(CFLAGS_DBG)
else
    CFLAGS := $(CFLAGS_BASE) $(CFLAGS_OPT)
endif

# Include directories
CFLAGS      += -I$(INC_DIR)

# Linker flags
LDFLAGS     := -lm

# ============================================================================
# Build Rules
# ============================================================================

# Default target
.PHONY: all
all: $(TARGET)

# Create object directory
$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

# Compile source files to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Link object files to create executable
$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Build complete: $(TARGET)"

# ============================================================================
# Maintenance Targets
# ============================================================================

.PHONY: clean
clean:
	@rm -rf $(OBJ_DIR)
	@rm -f $(TARGET)
	@echo "Clean complete"

.PHONY: rebuild
rebuild: clean all

.PHONY: help
help:
	@echo "FANSTAR Build System"
	@echo "===================="
	@echo ""
	@echo "Targets:"
	@echo "  make              - Build the test executable (default)"
	@echo "  make clean        - Remove all build artifacts"
	@echo "  make rebuild      - Clean and build"
	@echo "  make help         - Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  CC=<compiler>     - Specify C compiler (default: gcc)"
	@echo "  CFLAGS=<flags>    - Specify additional compiler flags"
	@echo "  DEBUG=1           - Build with debug symbols and no optimization"
	@echo ""
	@echo "Examples:"
	@echo "  make                          # Standard optimized build"
	@echo "  make DEBUG=1                  # Debug build"
	@echo "  make CC=clang                 # Use clang instead of gcc"
	@echo "  make rebuild                  # Clean and rebuild"
	@echo ""
	@echo "Output:"
	@echo "  Executable: $(TARGET)"
	@echo "  Objects:    $(OBJ_DIR)/"

.PHONY: info
info:
	@echo "Build Configuration:"
	@echo "  Compiler:    $(CC)"
	@echo "  CFLAGS:      $(CFLAGS)"
	@echo "  LDFLAGS:     $(LDFLAGS)"
	@echo "  Source dir:  $(SRC_DIR)"
	@echo "  Include dir: $(INC_DIR)"
	@echo "  Build dir:   $(BUILD_DIR)"
	@echo "  Target:      $(TARGET)"

# ============================================================================
# Dependency tracking (optional, for advanced builds)
# ============================================================================

# Uncomment the lines below to enable dependency tracking
# This allows make to rebuild when header files change

# DEPS := $(patsubst $(OBJ_DIR)/%.o, $(OBJ_DIR)/%.d, $(OBJECTS))
# 
# $(OBJ_DIR)/%.d: $(SRC_DIR)/%.c | $(OBJ_DIR)
# 	@$(CC) -MM $(CFLAGS) $< > $@
# 
# -include $(DEPS)

# ============================================================================
# Phony targets (not files)
# ============================================================================

.PHONY: all clean rebuild help info
