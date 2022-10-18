# O-tile structure construction when $m\equiv 2 \pmod{4}$
This program will try to help you construct a set of O-tile structures with $(\frac{m^2}{4}+m)$-tiles in $\mathbb{C}^m \otimes \mathbb{C}^m$ for $6\leq m\leq 38$ and $m\equiv 2 \pmod 4$.

And each constructed O-tile structure has $\frac{m^2-4m}{4}$ pieces of $2\times 2$ tiles, $m$ pieces of $1\times 2$ tiles and $m$ pieces of $2\times 1$ tiles.

- (1)	$\frac{m^2-4m}{4}$ pieces of $2\times 2$ tiles are divided into two cases: 

- case 1: $\frac{m}{2}$ pieces of $2\times 2$ tiles $\\{ j, j+\frac{m}{2}\\} \times \\{ pj, pj+\frac{m}{2}\\}, 0\leq j\leq \frac{m}{2}-1$.

- case 2: $\frac{m^2-6m}{4}$ pieces of $2\times 2$ tiles $M_{i,j}\times N_{i,j}$, where $M_{i,j} = \\{ j, j+l_{i}\\} , N_{i,j} = \\{ pj+a_{i}, pj+a_{i}+l_{i}\\} , 0\leq i\leq \frac{m-10}{4}, 0\leq j\leq m-1$.
- (2) $m$ pieces of $1\times 2$ tiles $\\{ j\\} \times N_{j}$, where $N_{j} = \\{ pj+x, pj+y\\} , 0\leq j\leq m-1$.
- (3) $m$ pieces of $2\times 1$ tiles $M_{j}\times \\{ j\\} $, where $M_{j} =\\{ p^{-1}(j-z), p^{-1}(j-w)\\} , 0\leq j\leq m-1$.

And $l_{i}, a_{i}, x, y, z$ and $w$ correspond to $len\[i\], arr\[i\], combo\[i\]\[0\], combo\[i\]\[1\], combo\[5-i\]\[0\]$ and $combo\[5-i\]\[1\]$ in this program respectively.


## Required:
- CodeBlocks, Visual Studio 2022 or Dev-C++


## To compile:
* Install gcc (g++) that supports -std=c++11 (only Dev-C++)  
Tools---Compiler Options---General---Add the following commands when calling the compiler: -std=c++11  


## Others:
- We present some running results of the program in folder 'Some results'. 
