package world

func matlabRound(x float64) int {
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

// interp1 interpolates to find yi, the values of the underlying function Y
// at the points in the vector or array xi. x must be a vector.
// http://www.mathworks.co.jp/help/techdoc/ref/interp1.html
//
// Input:
//   x          : Input vector (Time axis)
//   y          : Values at x[n]
//   xi         : Required vector
//
// Output:
//   yi         : Interpolated vector
func interp1(x, y []float64, xi []float64, yi []float64) {
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

// decimate carries out down sampling by both IIR and FIR filters.
// Filter coeffiencts are based on filterForDecimate().
//
// Input:
//   x          : Input signal
//   r          : Coefficient used for down sampling
//                (fs after down sampling is fs/r)
// Output:
//   y          : Output signal
func decimate(x []float64, r int, y []float64) {
	xLength := len(x)
	const nFact = 9
	tmp1 := make([]float64, xLength+nFact*2)
	tmp2 := make([]float64, xLength+nFact*2)

	for i := 0; i < nFact; i++ {
		tmp1[i] = 2*x[0] - x[nFact-i]
	}
	for i := nFact; i < nFact+xLength; i++ {
		tmp1[i] = x[i-nFact]
	}
	for i := nFact + xLength; i < 2*nFact+xLength; i++ {
		tmp1[i] = 2*x[xLength-1] - x[xLength-2-(i-(nFact+xLength))]
	}

	filterForDecimate(tmp1[:2*nFact+xLength], r, tmp2)
	for i := 0; i < 2*nFact+xLength; i++ {
		tmp1[i] = tmp2[2*nFact+xLength-i-1]
	}
	filterForDecimate(tmp1[:2*nFact+xLength], r, tmp2)
	for i := 0; i < 2*nFact+xLength; i++ {
		tmp1[i] = tmp2[2*nFact+xLength-i-1]
	}

	nout := (xLength-1)/r + 1
	nbeg := r - r*nout + xLength

	count := 0
	for i := nbeg; i < xLength+nFact; i += r {
		y[count] = tmp1[i+nFact-1]
		count++
	}
}

// filterForDecimate calculates the coefficients of low-pass filter and
// carries out the filtering. This function is only used for decimate().
func filterForDecimate(x []float64, r int, y []float64) {
	xLength := len(x)

	// filter Coefficients
	var a [3]float64
	var b [2]float64

	switch r {
	case 11: // fs : 44100 (default)
		a[0] = 2.450743295230728
		a[1] = -2.06794904601978
		a[2] = 0.59574774438332101
		b[0] = 0.0026822508007163792
		b[1] = 0.0080467524021491377
	case 12: // fs : 48000
		a[0] = 2.4981398605924205
		a[1] = -2.1368928194784025
		a[2] = 0.62187513816221485
		b[0] = 0.0021097275904709001
		b[1] = 0.0063291827714127002
	case 10:
		a[0] = 2.3936475118069387
		a[1] = -1.9873904075111861
		a[2] = 0.5658879979027055
		b[0] = 0.0034818622251927556
		b[1] = 0.010445586675578267
	case 9:
		a[0] = 2.3236003491759578
		a[1] = -1.8921545617463598
		a[2] = 0.53148928133729068
		b[0] = 0.0046331164041389372
		b[1] = 0.013899349212416812
	case 8: // fs : 32000
		a[0] = 2.2357462340187593
		a[1] = -1.7780899984041358
		a[2] = 0.49152555365968692
		b[0] = 0.0063522763407111993
		b[1] = 0.019056829022133598
	case 7:
		a[0] = 2.1225239019534703
		a[1] = -1.6395144861046302
		a[2] = 0.44469707800587366
		b[0] = 0.0090366882681608418
		b[1] = 0.027110064804482525
	case 6: // fs : 24000 and 22050
		a[0] = 1.9715352749512141
		a[1] = -1.4686795689225347
		a[2] = 0.3893908434965701
		b[0] = 0.013469181309343825
		b[1] = 0.040407543928031475
	case 5:
		a[0] = 1.7610939654280557
		a[1] = -1.2554914843859768
		a[2] = 0.3237186507788215
		b[0] = 0.021334858522387423
		b[1] = 0.06400457556716227
	case 4: // fs : 16000
		a[0] = 1.4499664446880227
		a[1] = -0.98943497080950582
		a[2] = 0.24578252340690215
		b[0] = 0.036710750339322612
		b[1] = 0.11013225101796784
	case 3:
		a[0] = 0.95039378983237421
		a[1] = -0.67429146741526791
		a[2] = 0.15412211621346475
		b[0] = 0.071221945171178636
		b[1] = 0.21366583551353591
	case 2: // fs : 8000
		a[0] = 0.041156734567757189
		a[1] = -0.42599112459189636
		a[2] = 0.041037215479961225
		b[0] = 0.16797464681802227
		b[1] = 0.50392394045406674
	default:
		a[0] = 0.0
		a[1] = 0.0
		a[2] = 0.0
		b[0] = 0.0
		b[1] = 0.0
	}

	// Filtering on time domain.
	var w [3]float64
	var wt float64
	for i := 0; i < xLength; i++ {
		wt = x[i] + a[0]*w[0] + a[1]*w[1] + a[2]*w[2]
		y[i] = b[0]*wt + b[1]*w[0] + b[1]*w[1] + b[0]*w[2]
		w[2] = w[1]
		w[1] = w[0]
		w[0] = wt
	}
}
