# depth-perception
(This README was added yearly 2 years after the project ended. Sorry for the lack of detail. I don't remember that much.)
## Dependencies
- STB library for opening images (https://github.com/nothings/stb/blob/master/stb_image.h)
- OpenCL for parallel processing
- Win32 for displaying renderings
## About Repo Contents
- The "Navigation Specification" PDF describes the purpose of this project, the systems at play in the "Navigation" folder's code, 
and some other features I didn't get to finish/start implementing. In the Contents page, I embolden the features that I did implement.
- The "Navigation" folder has the primary code written for this project, including the code (OpenCL kernel code and C++) for
performing depth perception and then projecting the measured points into a sparse 3D array of surfaces (stored as planes with colors).
- The "Visual Testing" folder has the code I used for visually testing my primary code. The main thing it does is render the plane map,
where the camera can be moved around to view the virtual world from different angles.
