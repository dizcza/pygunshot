# C++ implementation of pygunshot [1]

Muzzle blast and shock wave components are implemented.

## Compile and run

```
g++ gunshot.cpp usage.cpp -o gunshot.o
./gunshot.o
```

## Visualize the results in Python

```python
import matplotlib.pyplot as plt
import numpy as np
pmb = np.fromfile("output.float32", dtype=np.float32)
plt.plot(pmb)
plt.show()
```

## Known bugs

(Same for Python implementation.)

* Many pygunshot bugs have been fixed (sorted by severity):
  - [ainf - the speed of sound - is not a part of the Friedlander amplitude](https://github.com/metu-sparg/pygunshot/blob/76004698fc5b5f7c34012ff0983fd6da094d272a/pygunshot/muzzleblast.py#L27) according to [1,3]
  - [do not divide Pb by 100](https://github.com/metu-sparg/pygunshot/blob/76004698fc5b5f7c34012ff0983fd6da094d272a/pygunshot/muzzleblast.py#L158) according to [3]
  - [do not multiply lp by 10](https://github.com/metu-sparg/pygunshot/blob/76004698fc5b5f7c34012ff0983fd6da094d272a/pygunshot/muzzleblast.py#L113) according to [3]
  - [`coneAngle(M)`, not `coneAngle(uexit)`](https://github.com/metu-sparg/pygunshot/blob/76004698fc5b5f7c34012ff0983fd6da094d272a/pygunshot/process.py#L45) according to [3]
  - [wrong coefficients in the momentum index](https://github.com/metu-sparg/pygunshot/blob/76004698fc5b5f7c34012ff0983fd6da094d272a/pygunshot/muzzleblast.py#L133) according to [3]
  
  With these bugs fixed, however, the simulated muzzle blast at the mic is over-squeezed (the positive phase duration is very small). A dirty fix is to multiply the scaling length `l` by 10. A better fix would be a fit to real data.
  
* The default parameters of a Berlage wave claimed in [2] do not produce the nice picture shown in the paper.

## References

1. Hacıhabiboğlu, H. (2017). Procedural Synthesis of Gunshot Sounds Based on Physically Motivated Models. In Game Dynamics (pp. 47-69). Springer, Cham.
2. Dobrynin, Y., Maksymov, M., & Boltenkov, V. (2020). Development of a Method for Determining the Wear of Artillery Barrels by Acoustic Fields of Shots. Eastern-European Journal of Enterprise Technologies, 3(5), 105.
3. Fansler, K. S., Thompson, W. P., Carnahan, J. S., & Patton, B. J. (1993). A parametric investigation of muzzle blast. ARMY RESEARCH LAB ABERDEEN PROVING GROUND MD.
