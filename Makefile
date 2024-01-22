# For NVIDIA gpu
CC := clang++ 
CFLAGS := -fsycl -fsycl-targets=nvptx64-nvidia-cuda 

# For intel cpu
#CC := icpx
#CFLAGS := -fsycl 


#LDFLAGS=-L/usr/local/lib/vmd/plugins/LINUXAMD64/molfile/
#LDLIBS=-lmylib

SRCS := main.cpp InputData.cpp force_lj.cpp 
OBJS := $(SRCS:.cpp=.o)
EXE := test

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXE) $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXE)
