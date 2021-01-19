# Project: SMLC++ -- Sparse Matrix Library in C++ Using Old Yale Format.

CPP      = g++.exe
OBJ      = main.o sparse_matrix.o sparse_linearsys.o
LINKOBJ  = main.o sparse_matrix.o sparse_linearsys.o
LIBS     = -L"C:/MinGW64/lib32" -L"C:/dev_cpp/MinGW64/x86_64-w64-mingw32/lib32" -static-libgcc -g -m32 -g3
INCS     = -I"C:/MinGW64/include" -I"C:/dev_cpp/MinGW64/x86_64-w64-mingw32/include" -I"C:/MinGW64/lib/gcc/x86_64-w64-mingw32/4.8.1/include"
CXXINCS  = -I"C:/dev_cpp/MinGW64/include" -I"C:/dev_cpp/MinGW64/x86_64-w64-mingw32/include" -I"C:/dev_cpp/MinGW64/lib/gcc/x86_64-w64-mingw32/4.8.1/include" -I"C:/dev_cpp/MinGW64/lib/gcc/x86_64-w64-mingw32/4.8.1/include/c++"
BIN      = splib_pp.exe
CXXFLAGS = $(CXXINCS) -m32 -g3 -std=c++11 -g
RM       = rm.exe -f

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BIN) all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o $(BIN) $(LIBS)

main.o: main.cpp
	$(CPP) -c main.cpp -o main.o $(CXXFLAGS)

sparse_matrix.o: sparse_matrix.cpp
	$(CPP) -c sparse_matrix.cpp -o sparse_matrix.o $(CXXFLAGS)

sparse_linearsys.o: sparse_linearsys.cpp
	$(CPP) -c sparse_linearsys.cpp -o sparse_linearsys.o $(CXXFLAGS)
