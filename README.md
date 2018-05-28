# scicomp-project
Just a small 'final project' for python scientific computation course.

This consists of two things: 

Scripts generate_simple_geo.py and generate_spiral_geo.py can be used to generate 2D finite element domains. These are used in my research at the moment. New thing to me was the use of pygmsh. Gmsh also needs to be installed for full usage experience.

Scripts test_pod.py and pod.py do some model order reduction. I was planning to use scikit-learn surrogate models to find maximums of the error function but didn't get them to work. So I used scipy optimize routines instead.
