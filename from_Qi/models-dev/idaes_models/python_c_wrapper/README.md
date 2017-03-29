This directory includes the implementation of the Python-C wrapper.

**Directory contents:**

- wrapper: This directory contains the implementation of the wrapper and defines the macro FUNCGEN. It includes a test script that runs a dummy method and spits out some values.
- sampledir: This directory includes 2 python modules that define methods to be wrapped by the wrap script. These are used in the unit tests. 
- pyadolc_properties : This directory includes a python module that attempts to implement 2 property methods using pyadolc (work in progress).
- python_properties : This directory includes a python module that includes implementation of 2 property methods using substitution with direct values instead of AD. 
- prop_lib : This is a directory where the tests output a wrapped properties library. 
- autogenerate : This directory includes implementation of the wrap.py script and supporting wrapper.py. These provide basic functionality to wrap a directory of python modules into a C library to be used by pyomo.
- autogenerate/wrap_decorator : This directory defines the (dummy) @wrapthis decorator. 
- autogenerate/example : This directory contains a minimal pyomo model that uses the V_vap ExternalFunction. Use as template to test a generated .so library. 
- build_wrapper.sh (obsolete) : This script builds both the wrapper library and the minimal physical properties library and then runs the solver code.

**Basics:**

1. Any python method may be wrapped via the wrapper library and then used in a properties library (written in C). To start, you need to add the method you're attempting to wrap to a Python module. Say this is "phys_prop.py". 

2. The wrapper library is dependent on:

   * Python: add the path to python.h to your includes, link against lpython-2.7.
   * idaes_models : Needs to be setup already on your machine via setup.py.
   * ADOL-C / ASL / ipopt: Add the relevant include directories, link against adol-c. 

3. The wrapper library implements the FUNCGEN macro. This macro inputs a Python method's name and module and implements all the C/Python hooks necessary to call this method from C code. 

**Restrictions on Python methods to be wrapped**

1. All Python methods must input the following 3 arguments in the shown order:

   * n - Integer corresponding to arglist->n.
   * at - List corresponding to arglist->at.
   * ra - List corresponding to arglist-ra.
   
   Optionally, methods may also input the 2 following arguments in the shown order:
   
   * derivs - List corresponding to arglist->dervis. 
   * hes - List corresponding to arglist->hes.

   Example of an acceptable Python method prototype and return value:

   ```
   def V_vap_py(n, at, ra, derivs, hes):
       V = -1
       return [V,derivs,hes]
   ```

2. For return values, methods must return a list of the following format:

   ```
   [return_value,(optional) derivs, (optional) hes]
   ```

**Autogeneration**

To wrap a directory of Python module, clone the models repo then navigate to src/autogenerate and do the following:

1. In your Python module, import the @wrapthis decorator and use it as shown below:

   ```
   from idaes_models.python_c_wrapper.autogenerate.wrap_decorator.wrap_decorator import wrapthis

   @wrapthis
   def V_vap_py_first_five_args(n, at, ra, derivs, hes):
      ...
   ```

2. Add the directory of your Python module to the Python path so that the autogenerate script can see it, e.g:

   ```
   export PYTHONPATH=$PYTHONPATH:/home/jovyan/models/models-fork/models/src/sampledir/ 
   ```

3. Set 3 environment variables as follows:
    
   ```
   export MODELS_HOME=<models repo home directory>
   export CPPFLAGS="-I<ASL include directory> -I<Python include directory>"
   export ADOLC_LIBS=<ADOL-C lib directory>
   ```
   Example values: 

   ```
   export MODELS_HOME=~/models/models-fork/models
   export CPPFLAGS="-I/usr/local/src/Ipopt-3.12.5/ThirdParty/ASL -I/opt/conda/pkgs/python-2.7.12-0/include/python2.7"
   export ADOLC_LIBS=/usr/local/lib64
   ```

   You may also pass these values directly later on (step 5).

4. If on NERSC, load the needed modules by running the following:

   ```
   source $MODELS_HOME/src/autogenerate/setup_nersc_env.sh
   ```

5. Wrap:

   ```
   wrap --path ~/models/models-fork/models/src/sampledir --dest ~/models/models-fork/models/src/prop_lib
   ```

   This configures and compiles the wrapper library, autogenerates C code for a properties library based on your Python methods,
   and compiles that library placing it in the location pointed to by "dest". 

   Sample output:

   ```
   Properties library has been generated @ /home/jovyan/models/models-fork/models/src/prop_lib/phys_prop_C_auto.so
   Set environment variables using 'source /home/jovyan/models/models-fork/models/src/prop_lib/python_c_wrapper_helper.sh'
   ```

6. Run the helper script to set up the needed environment variables: 

   ```
   source /home/jovyan/models/models-fork/models/src/prop_lib/python_c_wrapper_helper.sh
   ```
   ...then point pyomo to the library location:

   ```
   model.V_func = ExternalFunction(library="/home/jovyan/models/models-fork/models/src/prop_lib/phys_prop_C_auto.so", function="V_vap_py")
   ```

**Manual Wrapping**

*Wrapper library compilation:*

3. Compile the wrapper library in the current Docker image:

   a. Browse to wrapper directory

      ```
      cd idaes_models/python_c_wrapper/wrapper

      ```

   b. Locate the ASL include directory (e.g: /usr/local/src/Ipopt-3.12.5/ThirdParty/ASL) and the Python include directory (e.g: /opt/conda/pkgs/python-2.7.12-0/include/python2.7). Also, locate the ADOL-C library location (e.g: /usr/local/lib64). Run the configure command passing the 2 include paths and the library path as shown:

      ```
      ./configure CPPFLAGS="-I/usr/local/src/Ipopt-3.12.5/ThirdParty/ASL \
      -I/opt/conda/pkgs/python-2.7.12-0/include/python2.7" --with-adolc=/usr/local/lib64

      ```
   c. Modify LD_LIBRARY_PATH path to include the paths to the ADOL-C libraries, for e.g in the Docker image:

      ```
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib64:$MODELS_HOME/idaes_models/python_c_wrapper/wrapper/.libs
      
      ```
   d. Compile the wrapper library using the existing Makefile.

      ```
      make
      ```
      A file libwrapper.so is generated under src/wrapper/.libs.

*Using the wrapper library in your C code (no autogeneration):*

4. Suppose your C code is located in a file called phys_prop.c. To define a method V_vap located in your phys_prop.py module in your wrapper library you will need to add the following lines:

   ```
   // Include the wrapper library header file and Python header file.
   #include <Python.h>
   #include "wrapper.h"
   ...
   ...
   // Wrap your python method V_vap using the wrapper library:
   FUNCGEN(phys_prop,V_vap);
   ```

6. In your funcadd function, you will define V_vap in the usual way. For e.g:

   ```
   void funcadd(AmplExports *ae)
   {
       int t = 0;
       addfunc("V_vap", (rfunc)V_vap_py, t, 2, NULL);

   }
   ```

8. Finally, when building your library, add the relevant include paths of the wrapper library and link against the wrapper library and Python, for e.g your INCLUDE and LD_FLAGS variables in your Makefile could look like this:

   ```
   INCLUDE = -I/home/jovyan/ADOL-C-2.6.2/ADOL-C/include \
   -I/usr/local/src/Ipopt-3.12.5/ThirdParty/ASL/solvers -I/usr/local/src/Ipopt-3.12.5/ThirdParty/ASL \
   -I/opt/conda/pkgs/python-2.7.12-0/include/python2.7 \
   -I$(MODELS_HOME)/idaes_models/python_c_wrapper/wrapper
   LDFLAGS = -L/home/jovyan/adolc_base/lib64 -ladolc -L$(MODELS_HOME)/src/wrapper/.libs -lwrapper -shared -lm -lpython2.7

   ```