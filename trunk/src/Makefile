include ../Makefile.inc

CXXFLAGS+= -I../include

LIB = liblshkit.a
SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,%.o,$(SRCS))
DEPS = $(patsubst %.o,%.d,$(OBJS))

.PHONY:	all clean

all:	$(LIB)

$(LIB):	$(OBJS)
	$(AR) $(ARFLAGS) $@ $^

.cpp.o:
	$(CXX) -c -MMD -MP $(CXXFLAGS) $< -o $@
	@sed -i -e '1s,\($*\)\.o[ :]*,\1.o $*.d: ,' $*.d

clean:
	#-$(MAKE) -C 3rd-party clean
	rm -f $(OBJS)
	rm -f $(DEPS)
	rm -f $(LIB)

-include $(DEPS)

