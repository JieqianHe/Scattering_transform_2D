# 2d-signal-synthesis

2D signal synthesis using wavelet scattering coefficients. Written in Python.

Authors: Jieqian He

### functions

The first four functions build up gabor wavelet family both in space and frequnecy.

The next three functions compute wavelet transforms in space and frequnecy.

The next two functions deal with 2D signal synthesis. Notice in function *synthesis*, we minimize a random conbination of wavelet coefficients iteratively until all wavelets are included. We use python build in function *least squares* to minimize difference while they simply use gradient descent with finite difference to find the local minimum.

### Results and funture plan
The synthesized image has wavelet coefficients close to the target but the image itself does not look like the target. The next thing is to add more information, such as the covariance of wavelet coefficients. We also need to specify the gradients of target function to save time.


