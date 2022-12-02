# AEPH Input Guide

This input guide is a companion of MATLAB program **AEPH** from *Anisotropic elasticity with Matlab*, Chyanbin Hwu, Springer, 2021. Additional information can be found in Sec. 3.1.4 and App. E in the book.

> All the labels used in AEPH are **numbers**, i.e. all data in the input files must be numeric.

> Parameters enclosed by brackets `[ ]` are optional.

## Table of Contents
- [input_control.txt](#input_controltxt)
- [input_elastic.txt](#input_elastic1txt)
- [input_loadstr.txt](#input_loadstrtxt)
  - [Ltype=4](#ltype4-bem-for-elliptical-hole)
  - [Ltype=9](#ltype9-bem-for-polygonal-hole)
  - [Ltype=411](#ltype411-uniform-load)
  - [Ltype=611](#ltype611-uniform-load-elliptical-hole)
  - [Ltype=614](#ltype614-point-load-elliptical-hole)
  - [Ltype=622](#ltype622-uniform-load-polygon-like-hole)
- [input_variable.txt](#input_variabletxt)
- [input_xn.txt](#input_xntxt)
- [input_node1.txt](#input_node1txt)
- [input_bc.txt](#input_bctxt)

## input_control.txt

### Format

```tex
Nmat Dtype Ltype Otype Icheck E0 h0 eps
Etype Ttype Ptype Vtype
[Etype Ttype Ptype Vtype]
[...]
```

Those in the 1st line are global control parameters.
The material property parameters start from the 2nd line. **Each line is a different material.**

> Parameters enclosed by brackets `[ ]` are optional.

E.g.

```matlab
1 2 611 2 0 1e9 1e-3 1e-6
1 0 0 0
```

### Parameters

#### Nmat
Total number of materials.
  
#### Dtype
Valid problem dimension labels are:
  - **1, 11, 111, 112, 12, 121, 122**  
    Generalized plane strain.
  - **2, 21, 211, 212, 22, 221, 222**  
    Generalized plane stress.
  - **3, 31, 311, 312, 32, 321, 322**  
    Coupled stretching-bending.
  - **4, 14, 104**  
    Three-dimensional.

#### Ltype
Solution label, often corresponds to its section number in the book.
  
#### Otype
Output points label. Used with `input_variable.txt`, see [that section](#input_variabletxt) for more information. Valid labels are:
- **1, 11, 12, 13**: Curve  
  Additional output for analytical solutions (**Ltype > 100**):
  - **11**: $\sigma_{ss}$, hoop stress (hole)
  - **12**: $\sigma_{ss}, \sigma_{nn}, \sigma_{sn}$, internal surface traction (inclusion/interface?)
  - **13**: $\sigma_{nn}$, normal stress (contact?)
- **2**: Area.
- **3**: Discrete points.
  
#### Icheck
  Internal check label. If and only if `Icheck=1`, values of **A**, **B**, **mu**, **N1**, **N2**, **N3**, **N** are verified through alternative approach and some identities. In some cases, also verifies the solution.
 
#### E0
  Reference Young's modulus for the nondimensionalization of **N2**, **N3**, and **N**. Suggested value=<mark>1.0e9</mark>.
  
#### h0
  Reference thickness used for the nondimensionalization of Aij,Bij and Dij. Suggested value=<mark>1.0e-3</mark>.
  
#### eps
  Perturbation ratio and error tolerance, suggested value=<mark>1.0e-6</mark>.
  
  > The suggested values of E0, h0, and eps are taken from comments in `Main.m`.
  
#### Etype
  Used with `input_elastic#.txt`, see [that section](#input_elastic1txt) for more informations. Valid elastic property labels:
  
  - **0**  
    No elastic properties.
    
    > Even though it is listed here, it does not mean you have to use it. Think again before doing anything stupid.
    
  - **1**  
    Isotropic.
    
  - **2**  
    Orthotropic.
    
  - **3**  
    Anisotropic, values in `input_elastic#.txt` are Cij.
    
  - **4**  
    Anisotropic, values in `input_elastic#.txt` are Sij.
    
  - **5**  
    Unidirectional fiber-reinforced composite.
    
  - **6**  
    Composite laminate.
    
#### Ttype
Used with `input_thermal#.txt`, see that section for more information. (Unfortunately, because I have never run any thermal problem with AEPH, "that section" **does not** exist yet. Same for Ptype and Vtype.) Valid thermal property labels:
  - **0**  
    No thermal properties.
  - **1**  
    Isotropic.
  - **2**  
    Orthotropic.
  - **3**  
    Anisotropic, values in `input_thermal#.txt` are kij and betaij.
  - **4**  
    Anisotropic, values in `input_thermal#.txt` are kij and alphaij.
  - **5**  
    Unidirectional fiber-reinforced composites.
  - **6**  
    Composite laminates.
    
#### Ptype
  Piezoelectric properties input label. Used with `input_piezo#.txt`, see Ch. 11 of the book for more information. Valid labels are:
  - **0**  
    No piezoelectric properties.
  - **1, 2, 3, 4**
  - **5**  
    Electro-elastic laminates.
  - **11, 12, 13, 14, 15, 16, 17, 18, 19**  
    Magneto-electro-elastic materials.
    
#### Vtype
  Viscoelastic properties input label. Valid labels are:
  - **0**  
    No viscoelastic properties.
  - **1, 2**  
    Isotropic.
  - **3**  
    Standard linear viscoelastic solids.
  - **4**  
    Prony series.

## input_elastic1.txt

Here goes the elastic properties.

> Each material has a input file; there are as many `input_elastic#.txt` files as the materials.
> The number **1** in the file name is its material number (order of appearance in `input_control.txt`).

Its format depends on [Etype](#etype) in `input_control.txt`.

### Etype=1, Isotropic

```
E v
```

- **E**: Young's modulus.
- **v**: Poisson's ratio.

### Etype=2, Orthotropic

```
E1 E2 E3 G23 G31 G12 v23 v13 v12 angle
```

- **E1, E2, E3**: Young's moduli.
- **G23, G31, G12**: shear moduli.
- **v23, v13, v12**: Poisson's ratios.
- **angle**: orientation in x1-x2 plane in degrees, directed counterclockwisely from positive x1-axis to the principal material direction.

### Etype=3, General anisotropic, Cij

```
C11 C21 C31 C41 C51 C61
C12 C22 C32 C42 C52 C62
C13 C23 C33 C43 C53 C63
C14 C24 C34 C44 C54 C64
C15 C25 C35 C45 C55 C65
C16 C26 C36 C46 C56 C66
```

- **Cij**: elements of stiffness matrix.

### Etype=4, General anisotropic, Sij

```
S11 S21 S31 S41 S51 S61
S12 S22 S32 S42 S52 S62
S13 S23 S33 S43 S53 S63
S14 S24 S34 S44 S54 S64
S15 S25 S35 S45 S55 S65
S16 S26 S36 S46 S56 S66
```

- **Sij**: elements of compliance matrix.

### Etype=5, Unidirectional fiber-reinforced composite

```
E1 E2 G12 v12 angle
```

- **E1, E2,**: Young's moduli.
- **G12**: shear modulus.
- **v12**: Poisson's ratio.
- **angle**: fiber orientation in degrees, directed counterclockwisely from x1-axis to the principal material direction.

### Etype=6, Composite laminates

The first line defines the number of layers and materials, followed by elastic properties of each material, and the rest are layer arrangements. The first layer has the smallest `z`. The reference plane is the mid-plane.

```tex
nLayer nMat
E1 E2 G12 v12
[E1 E2 G12 v12]
[...]
mat angle thk
[mat angle thk]
[...]
```

E.g.

```
4 2
138e9 9e9 6.9e9 0.3
9e9 9e9 6e9 0.25
1 45 0.001
2 0 0.001
2 45 0.001
1 -45 0.001
```

In this example, there are 4 layers and 2 materials.
Material 1 is defined by `138e9 9e9 6.9e9 0.3`, and material 2 by `9e9 9e9 6e9 0.25`.
Layer 1 is defined by `1 45 0.001`; it is made of material **1** with fiber angle **45°** and thickness **0.001**, and the range of its `z` coordinate is [-0.002, -0.001].

- **nLayer**: total number of layers.
- **nMat**: total number of materials.
- **E1, E2**: Young's moduli.
- **G12**: shear modulus.
- **v12**: Poisson's ratio.
- **mat**: material number of this layer. e.g. `mat=1` is the 1st material defined in this file.
- **angle**: fiber orientation in degrees, directed counterclockwisely from x1-axis to the principal material direction.
- **thk**: layer thickness.

## input_loadstr.txt

`input_loadstr.txt` provides the parameters used in the solutions, e.g. value of stress at infinity, location and value of point force, and parameters of elliptical hole/inclusion; therefore, its format varies for each `Ltype`. `loadstr` stands for "load and structure".

> Only those I use frequently are listed here, for a full list, see appendix E.4 of *Anisotropic elastic plate with MATLAB*.

Its format depends on [Ltype](#ltype) in `input_control.txt`.

### Ltype=4, BEM for elliptical hole

```
elemType GausPts x0 y0 angle a b
```

- **elemType**: defines the element type. Valid labels are:
  
  - **1**: linear element.
    
  - **2**: quadratic element.
    
  - **3**: linear element with cubic deflection.
    
- **GausPts**: number of Gaussain points used for the line integral.
  
- **x0**, **y0**: center of ellipse.
  
- **angle**: angle from positive x1-axis to the major axis of the ellipse, in degrees.
  
- **a**, **b**: half major and minor axes of the ellipse.
  

### Ltype=9, BEM for polygonal hole

```
elem GausPts x0 y0 angle a b epsilon k
```

- **elem**: element type. Valid labels are:
  
  - **1**: linear element.
    
  - **2**: quadratic element.
    
  - **3**: linear element with cubic deflection.
    
- **GausPts**: number of Gaussain points used for the line integral.
  
- **x0**, **y0**: center of ellipse.
  
- **angle**: angle from positive x1-axis to the major axis of the ellipse, in degrees.
  
- **a, c, epsilon, k**: parameters defining hole contour by
  
  $$
  x_1=a(\cos\psi+\epsilon\cos k\psi)\,, \qquad
  x_2=a(c\sin\psi-\epsilon\sin k \psi)\,.
  $$

### Ltype=411, Uniform load

```
loadLabel val1 [val2 val3 ...]
```

Avaliable `loadLabel` and their corresponding `val` are:

- **Unidirectional tension with an angle, loadLabel=1**
  
  ```tex
  1 sigma angle
  ```
  
  - **sigma**: value of unidirectional tension.
  - **angle**: angle of the load, in degrees, counterclockwise from positive x1-axis.
- **Biaxial tension, loadLabel=2**
  
  ```tex
  2 sigma1 sigma2
  ```
  
  - **sigma1, sigma2**: value of tension loads.
- **Inplane shear, loadLabel=3**
  
  ```tex
  3 sigma
  ```
  
  - **sigma**: value of inplane shear load.
- **Antiplane shear, loadLabel=4**
  
  ```tex
  4 sigma13 sigma23
  ```
  
  - **sigma13, sigma23**: value of anti-plane shear load $\tau_{13}$ and $\tau_{23}$.
- **Stress components, loadLabel=5**
  
  ```tex
  5 sigma11 sigma22 sigma23 sigma13 sigma12
  ```
  
  - **sigma11, sigma22, sigma23, sigma13, sigma12**: stress at infinity $\sigma_{11}, \sigma_{22}, \sigma_{23}, \sigma_{13}, \sigma_{12}$.

### Ltype=611, Uniform load, elliptical hole

```
a b loadLabel val1 [val2 val3 ...]
```

- **a, b**: semi-major and minor axes of the ellipse
  
- **loadLabel**: valid labels are: **1, 2, 3, 4, 5**. See [Ltype 411](#ltype411-uniform-load) for loadLabels and their corresponding `val1 val2 ...`.
  

### Ltype=614, Point load, elliptical hole

```
a b p1 p2 p3 [p4] x1 x2
```

- **a, b**: semi-major and minor axes of the ellipse
  
- **p1, p2, p3, p4**: components of point force, `p4` for piezoeletric material.
  
- **x1, x2**: (x1, x2) is the location of the point force.
  

### Ltype=622, Uniform load, polygon-like hole

```
a c epsilon k loadLabel val1 [val2 val3 ...]
```

- **a, c, epsilon, k**: parameters defining hole contour by
  
  $$
  x_1=a(\cos\psi+\epsilon\cos k\psi)\,, \qquad
  x_2=a(c\sin\psi-\epsilon\sin k \psi)\,.
  $$
- **loadLabel**: valid labels are: **1, 2, 3, 4, 5**. See [Ltype 411](#ltype411-uniform-load) for loadLabels and their corresponding `val1 val2 ...`.
  

## input_variable.txt

Despite its mysterious name, `input_variable.txt` defines the points at which the results (e.g. displacement, stress, strain) are calculated and output to `result.txt` and figures.

Its format depends on [Otype](#otype) in `input_control.txt`.  
This is only a portion of all the available options. Consult *Anisotropic elastic plate with MATLAB* for the full list.

> The first number in `input_variable.txt` controls the type of curve/surface. For example, this number is `1` for piecewise line segments, and `2` for circular area. Don't forget to include it in this file.

### Otype=1, 11, 12, 13, Curve

- **Piecewise line segments**
  
  ```tex
  1 nLines x1 y1 x2 y2 nPts [x1 y1 x2 y2 nPts ...]
  ```
  
  - **nLines**: total number of line segments.
  - **x1, y1, x2, y2**: the line segment starts from (x1, y1) and ends at (x2, y2).
  - **nPts**: total number of points on this line segment, including its end points.
- **Arc**
  
  ```tex
  2 x0 y0 r start_angle end_angle nPts
  ```
  
  - **x0, y0**: center of the arc.
  - **r**: radius.
  - **start_angle, end_angle**: starting and ending angles of the arc, in degrees.
  - **nPts**: total number of points on the arc, end points included.
- **Slanted elliptical curve**
  
  ```tex
  3 x0 y0 a b start_angle end_angle nPts slant_angle
  ```
  
  - **x0, y0**: center of the ellipse
  - **a, b**: lengths of half major and minor axes.
  - **start_angle, end_angle**: range of the ellipse parameter $\psi$, in degrees.
  - **slant_angle**: rotation angle of the curve, counterclockwise from positive x1-axis.
- **Slanted polygon-like curve**
  
  ```tex
  4 x0 y0 a c epsilon k psi_0 psi_1 nPts slant_angle
  ```
  
  - **x0, y0**: center of the ellipse
    
  - **a, c, epsilon, k**: parameters of the curve
    
    $x_1 = a(\cos\psi+\epsilon\cos k\psi),\ x_2=a(c\sin\psi-\epsilon\sin k\psi).$
    
  - **psi_0, psi_1**: range of the ellipse parameter $\psi$, in degrees.
    
  - **slant_angle**: rotation angle of the curve, counterclockwise from positive x1-axis.
    

### Otype=2, Area

- **Slanted rectangle**
  
  ```tex
  1 x1 y1 x2 y2 nXPts nYPts angle
  ```
  
  - **x1, y1, x2, y2**: two corners, (x1, y1) and (x2, y2), defining the rectangle.
    
  - **nXPts, nYPts**: number of points divided in x and y directions; `nXPts*nYPts` points in total.
    
  - **angle**: slant angle, in degrees, counterclockwise. I guess it rotates the rectangle, defined by x1, y1, x2, and y2, with respect to its center.
    
- **Circular**
  Disk, ring, or sector.
  
  ```tex
  2 x0 y0 r1 theta1 r2 theta2 nRPts nThPts
  ```
  
  - **x0, y0**: center of the disk/ring.
  - **r1, r2**: inner and outer radii.
  - **theta1, theta2**: starting and ending angle, in degrees.
  - **nRPts, nThPts**: number of points divided in r and theta directions.

### Otype=3, User defined points

```tex
nPts
x y
[x y]
[...]
```

- **nPts**: total number of points.
  
- **x, y**: coordinate of the points.
  

## input_xn.txt

Coordinates of nodes. Arranged **clockwise**. Each row contains the coordinates of one node. The nodes are numbered automatically from `1` to `total number of nodes`.

```
x1 x2 [x3]
x1 x2 [x3]
[...]
```

## input_node1.txt

Though it is called `input_node1.txt`, this file defines the **elements** by nodal connectivity. Each row defines one element.

> The number **1** in the file name is its region number.

```
n1 n2 [n3 n4 n5 n6 n7 n8]
n1 n2 [n3 n4 n5 n6 n7 n8]
[...]
```

- **n1**: index of the first node of the element.
  
- **n2**, **n3**, **n4**, **n5**, **n6**, **n7**, **n8**: index of the second (n2) to eighth node (n8) of the element, if any.
  

## input_bc.txt

Defines the boundary conditions of nodes. Must be defined for **all** the nodes.

```
label1 label2 label3 [label4 ...] value1 value2 value3 [value4 ...]
```

- **label1**, **label2**, **label3**, **...** : type of boundary condition of the 1st to 8th degree of freedom, if any. Valid labels are:
  
  - **0**: traction-prescribed.
    
  - **1**: displacement-prescribed.
    
  - **2**: contact.
    
- **value1**, **value2**, **value3**, **...** : prescribed values of the 1st to 8th degree of freedom, if any.
