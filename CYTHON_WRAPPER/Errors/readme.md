### 1. Clang:

we tried to compile cython wrapper with Clang, and it seems Clang cannot handle overwriting the same methods which were extracted from different c++ header files; therefore, it produces this error. 

<img src="Clang.png" width="900" height="450">



### 2. Constructors :

the most problemtic part of the code belongs to the second contructors in general when the initialization of the memory is not nullary. At this point we got a bit confused whether to use c++ types in the arguments of constructors infront of new function or not. New function returns a pointer; therefore, everything infront of it should be set as a c++ object rather than a python one. Therefore we left the constructor as it is (without any type ) infront of the new function, in order to revoke the C++ constructor which was declared in the header file.

On the other hand, the __cinit__ in the Binstar class accepts every type both python and c++, but we have to be careful whether we want to put python objects or nor, since python objects affect the performance of the programme.



<img src="con.png" width="900" height="200">
