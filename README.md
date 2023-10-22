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
## Demo
The visual test executable (compiled for windows) can be found in the "demo1" release, and screenshots from the test can be found in
the "Screenshots" folder.
In this test, two images (MCImg1 and MCImg2) represent pictures of the real world taken by two adjacent and parallel cameras 
(as described in the specification). Next, the depths of various points throughtout the perceived region are calculated, and the points
are projected into a virtual world where the origin marks the position of one of the two parallel cameras. This virtual world is
percieved by a virtual camera, which the user can move around forward/left/backward/right/up/down using the w/a/s/d/c/x keys, 
respectively. By dragging the mouse, the camera can be rotated, and the t key can toggle between a view of the plane map and a view of
the triangle mesh generated immediately after depth perception is performed. Many pairs of the screenshots are taken from the same 
position/angle, but they look different because they were rendered in different modes using the t key. The new view, after the camera
is moved or rotated, isn't rendered until the space key is pressed.
