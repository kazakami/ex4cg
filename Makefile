objs =  Vector.o Matrix4x4d.o Quaternion.o image.o vrml.o imageLoader.o pmxLoader.o vmdLoader.o
CXXFLAGS = -Wall -std=c++11 -O3 -fopenmp# -pg
LDLIBS	= -lm -lglut -lGLU -lGL -lpng -lXext -lX11 -fopenmp

.SUFFIXES: .cpp .o
.PHONY: all last clean


all: kadai1 kadai2 test_gl pmx_test last
last: last_pmx last_vrml


kadai1: kadai1.o $(objs)
	$(CXX) -o $@ $< $(objs) $(LDFLAGS) $(LDLIBS)
kadai2: kadai2.o $(objs)
	$(CXX) -o $@ $< $(objs) $(LDFLAGS) $(LDLIBS)
test_gl: test_opengl.o $(objs)
	$(CXX) -o $@ $< $(objs) $(LDFLAGS) $(LDLIBS)
pmx_test: pmx_test.o $(objs)
	$(CXX) -o $@ $< $(objs) $(LDFLAGS) $(LDLIBS)

last_vrml: last_vrml.o $(objs)
	$(CXX) -o $@ $< $(objs) $(LDFLAGS) $(LDLIBS)
last_pmx: last_pmx.o $(objs)
	$(CXX) -o $@ $< $(objs) $(LDFLAGS) $(LDLIBS)



.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<


clean: 
	$(RM) *.o kadai1 kadai2 test_gl pmx_test
