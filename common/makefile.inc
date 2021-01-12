ifneq (,$(wildcard /usr/local/opt/llvm/bin/clang++))
  CXX = /usr/local/opt/llvm/bin/clang++
  CXXFLAGS += -I/usr/local/opt/llvm/include
  LDFLAGS  += -L/usr/local/opt/llvm/lib
else
  CXX ?= g++
  CXXFLAGS += -fno-loop-interchange -fno-loop-unroll-and-jam
endif

CXXFLAGS += -std=c++11 -O3 -fopenmp \
            -march=native $(LDFLAGS)

clang-tidy: $(SRCS)
	clang-tidy $(filter %.cpp, $^) -- -Wall

tar: submission.tar.gz

submission.tar.gz: $(REQUIRED_FILES)
	tar czvf submission.tar.gz $^

check:
	@tar tf submission.tar.gz $(REQUIRED_FILES)
