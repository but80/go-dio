package matlab

// Round calculates rounding.
func Round(x float64) int {
	if 0 < x {
		return int(x + .5)
	}
	return int(x - .5)
}

// histc counts the number of values in vector x that fall between the
// elements in the edges vector (which must contain monotonically
// nondecreasing values). n is a length(edges) vector containing these counts.
// No elements of x can be complex.
// http://www.mathworks.co.jp/help/techdoc/ref/histc.html
//
// Input:
//   x              : Input vector
//   edges          : Input matrix (1-dimension)
//
// Output:
//   index          : Result counted in vector x
// Caution:
//   Lengths of index and edges must be the same.
func histc(x []float64, edges []float64, index []int) {
	count := 1

	i := 0
	for ; i < len(edges); i++ {
		index[i] = 1
		if edges[i] >= x[0] {
			break
		}
	}
	for ; i < len(edges); i++ {
		if edges[i] < x[count] {
			index[i] = count
		} else {
			index[i] = count
			i--
			count++
		}
		if count == len(x) {
			break
		}
	}
	count--
	for i++; i < len(edges); i++ {
		index[i] = count
	}
}

// Interp1 interpolates to find yi, the values of the underlying function Y
// at the points in the vector or array xi. x must be a vector.
// http://www.mathworks.co.jp/help/techdoc/ref/Interp1.html
//
// Input:
//   x          : Input vector (Time axis)
//   y          : Values at x[n]
//   xi         : Required vector
//
// Output:
//   yi         : Interpolated vector
func Interp1(x, y []float64, xi []float64, yi []float64) {
	h := make([]float64, len(x)-1)
	s := make([]float64, len(xi))
	k := make([]int, len(xi))

	for i := 0; i < len(x)-1; i++ {
		h[i] = x[i+1] - x[i]
	}

	histc(x, xi, k)

	for i := range xi {
		s[i] = (xi[i] - x[k[i]-1]) / h[k[i]-1]
	}

	for i := range xi {
		yi[i] = y[k[i]-1] + s[i]*(y[k[i]]-y[k[i]-1])
	}
}
