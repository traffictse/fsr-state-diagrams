# Usage
First, you need to import the target FSR in either `SuccCalQuick.m` (for feedback functions in a recursion form) or `SuccCalFF.m` (for feedback functions in a polynomial form on finite fields). For example, for a degree-4 FSR `f(x) = x_1 + x_2 + x_4`:
* In `SuccCalQuick.m`, it should be `matrix(:, degree + 1) = mod(matrix(:, 1) + matrix(:, 2) + matrix(:, degree), 2);`.
* In `SuccCalFF.m`, there is an old bug. Not recommended currently.

Second, you could simply use `SubgCalQuick(degree,option_plot,option_cal)` to compute and plot the state diagram and its properties.
* option_plot  = 1, 2, 3, 4
  * 1: plot the state diagram
  * 2: plot the conjugated matrix and the associated matrix as two graphs
  * 3: plot the state diagram and two matrices as three graphs
  * 4: no plotting
* option_cal = 1, 2
  * 1: target FSR in `SuccCalQuick.m`
  * 2: target FSR in `SuccCalFF.m`

All parameters would be exported in a file `datas.txt`, including detailed structure of each subgraph.

For the option `layout` of the `plot` function, `layered` and `force` make no big difference. Still, `force` makes it more like a brocoli. 
