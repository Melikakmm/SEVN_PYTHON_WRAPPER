
# Extern block:

extern doesn't wrap the functions or variables automatically , it just calls them from the source or the header file. it is equivalent to #include <header.h>. 

Macros should be defined in extern blocks like either real functions or variables. 

semicolons must be ommited.

class declarations: because we don't want to use these classes like IO as an opaque type, we don't declare the body of the class, we just declare the name.


# extension types :

The conventional way to wrap a C++ class in Cython is with an extension type. We name it RNG to avoid clashing with the MT_RNG name, although there are ways to allow them to have the same name (see Chapter 6). Typically, a wrapper extension type has a pointer to a heap-allocated instance of the C++ class it is wrapping.


Storing a pointer to a heap-allocated C++ object in an extension type works in all instances. If the C++ class provides a nullary (no- argument) constructor, we can store a stack-allocated object directly —that is, no pointer indirection required. This removes the need to allocate and delete the instance, and there are efficiency gains as well.



An instance of a C++ class is an object that is created based on the definition of a C++ class. In C++, a class is a blueprint or template for creating objects, and an instance of a class is an actual object that is created in memory. The instance of a C++ class contains data members and member functions defined in the class. When an instance of a C++ class is created, memory is allocated for the data members of the object, and the constructor of the class is called to initialize the object's state. The instance can then be used to call member functions and access data members. When the instance is no longer needed, its destructor is called to release the memory used by the object.

