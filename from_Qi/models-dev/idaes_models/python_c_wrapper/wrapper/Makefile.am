# Header files
include_HEADERS = wrapper.h

# Use libtool to create a library in the lib directory
lib_LTLIBRARIES = libwrapper.la

# Specify the library sources
libwrapper_la_SOURCES =  wrapper.c
    
# Add extra libraries to this. Using -x c++ to force all files to be compile as c++
libwrapper_la_CFLAGS= $(CPPFLAGS) -c -Wall -Wextra -fPIC
    
# Make this a shared library
libwrapper_la_LDFLAGS= $(ADOLC_LIBS) -ladolc -shared -lm $(PYTHON_LIBS) -lpython2.7
