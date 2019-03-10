# Fractal Flames

This is an implementation of the fractal flames algorithm by Scott Draves. 
I originally had the repo at https://bitbucket.org/alrikai/fractalflames/wiki/Home but have since migrated it here. 

It spawns its own threads and operates asynchronously, generating images and writing the output frames to a shared thread safe buffer that the caller registers (see main.cpp for an example of this). 
