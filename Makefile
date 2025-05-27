# Компилятор и флаги
CXX = g++
CXXFLAGS = -std=c++11 -Wall
PARALLEL_FLAGS = -fopenmp

# Исходные файлы (.cpp)
CPP_SOURCES = A0_START.cpp Timer.cpp

# Заголовочные файлы (.h)
HEADERS = $(wildcard *.h)

# Имя исполняемого файла
TARGET = ./x64/program
PARALLEL_TARGET = ./x64/parallel

# Правила компиляции
all: $(TARGET)

$(TARGET): $(CPP_SOURCES) $(HEADERS)
	mkdir -p x64
	$(CXX) $(CXXFLAGS) $(CPP_SOURCES) -o $(TARGET)

parallel: $(PARALLEL_TARGET)

$(PARALLEL_TARGET): $(CPP_SOURCES) $(HEADERS)
	mkdir -p x64
	$(CXX) $(CXXFLAGS) $(PARALLEL_FLAGS) $(CPP_SOURCES) -o $(PARALLEL_TARGET)

clean:
	rm -f $(TARGET) $(PARALLEL_TARGET)
