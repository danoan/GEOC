## GEOC

Geometric Estimators on digital Curves (GEOC) library computes some
multigrid convergent estimators along digital curves. Currently 
implemented are: Local length, curvature, tangent. This library is part
of the project [BTools](https://github.com/danoan/BTools), developped to
support the publication [1].

 [1] Antunes, D., Lachaud, J.O., Talbot H.: Digital 
 curvature evolution model for image segmentation. In:
 Couprie, M., Cousty, J., Kenmochi, Y., Mustafa, N. (eds.) 
 Discrete Geometry for Computer Imagery. pp 15-26. Springer
 International Publishing, Cham (2019).
 
## Dependencies

1. [libboost1.66.0-dev](https://www.boost.org/users/history/version_1_66_0.html)
2. [opencv-3.3.0](https://opencv.org/releases.html)
3. [eigen-3.36](http://eigen.tuxfamily.org/index.php?title=Main_Page)
4. [DGtal0.9](https://dgtal.org/download/)

# Installation
```
mkdir build
cd build
cmake ..
(ccmake ..) - environment variables configuration
make install
```