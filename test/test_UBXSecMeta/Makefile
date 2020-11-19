.PHONY: clean

CC = g++
CFLAGS = -Wall -Wextra -Werror -Wno-sign-compare -Wshadow -Wdouble-promotion -std=c++17

INCLUDES = -I$(UBSMEAR_DIR)/inc/ -I$(UBSMEAR_DIR)/test/
TARGET   = test

$(TARGET): $(TARGET).cc
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(TARGET).cc

clean:
	rm -f $(TARGET)
