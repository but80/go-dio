package world

import (
	"gonum.org/v1/gonum/fourier"
)

// Commands for FFT (This is the same as FFTW)
const (
	fftForward = 1 + iota
	fftBackward
	fftEstimate
)

// Struct used for FFT
type fftPlan struct {
	fft  *fourier.FFT
	n    int
	sign int
	// flags uint
	cIn  []complex128
	in   []float64
	cOut []complex128
	out  []float64
	// input []float64
	// ip    []int
	// w     []float64
}

func backwardFFT(p fftPlan) {
	p.fft.Sequence(p.out, p.cIn[:p.n/2+1])
}

func forwardFFT(p fftPlan) {
	p.fft.Coefficients(p.cOut[:p.n/2+1], p.in)
}

// func originalBackwardFFT(p fftPlan) {
// 	if p.cOut == nil {
// 		// c2r
// 		p.input[0] = real(p.cIn[0])
// 		p.input[1] = real(p.cIn[p.n/2])
// 		for i := 1; i < p.n/2; i++ {
// 			p.input[i*2] = real(p.cIn[i])
// 			p.input[i*2+1] = -imag(p.cIn[i])
// 		}
// 		rdft(p.n, -1, p.input, p.ip, p.w)
// 		for i := 0; i < p.n; i++ {
// 			p.out[i] = p.input[i] * 2.0
// 		}
// 	} else {
// 		// c2c
// 		for i := 0; i < p.n; i++ {
// 			p.input[i*2] = real(p.cIn[i])
// 			p.input[i*2+1] = imag(p.cIn[i])
// 		}
// 		cdft(p.n*2, -1, p.input, p.ip, p.w)
// 		for i := 0; i < p.n; i++ {
// 			p.cOut[i] = complex(
// 				p.input[i*2],
// 				-p.input[i*2+1],
// 			)
// 		}
// 	}
// }

// func originalForwardFFT(p fftPlan) {
// 	if p.cIn == nil {
// 		// r2c
// 		for i := 0; i < p.n; i++ {
// 			p.input[i] = p.in[i]
// 		}
// 		rdft(p.n, 1, p.input, p.ip, p.w)
// 		p.cOut[0] = complex(p.input[0], 0.0)
// 		for i := 1; i < p.n/2; i++ {
// 			p.cOut[i] = complex(
// 				p.input[i*2],
// 				-p.input[i*2+1],
// 			)
// 		}
// 		p.cOut[p.n/2] = complex(p.input[1], 0.0)
// 	} else {
// 		// c2c
// 		for i := 0; i < p.n; i++ {
// 			p.input[i*2] = real(p.cIn[i])
// 			p.input[i*2+1] = imag(p.cIn[i])
// 		}
// 		cdft(p.n*2, 1, p.input, p.ip, p.w)
// 		for i := 0; i < p.n; i++ {
// 			p.cOut[i] = complex(
// 				p.input[i*2],
// 				-p.input[i*2+1],
// 			)
// 		}
// 	}
// }

func fftPlanDftC2r1d(n int, in []complex128, out []float64, flags uint) fftPlan {
	var output fftPlan
	output.fft = fourier.NewFFT(n)
	output.n = n
	output.in = nil
	output.cIn = in
	output.out = out
	output.cOut = nil
	output.sign = fftBackward
	// output.flags = flags
	// output.input = make([]float64, n)
	// output.ip = make([]int, n)
	// output.w = make([]float64, n*5/4)

	// output.ip[0] = 0
	// makewt(output.n>>2, output.ip, output.w)
	// makect(output.n>>2, output.ip, output.w[output.n>>2:])
	return output
}

func fftPlanDftR2c1d(n int, in []float64, out []complex128, flags uint) fftPlan {
	var output fftPlan
	output.fft = fourier.NewFFT(n)
	output.n = n
	output.in = in
	output.cIn = nil
	output.out = nil
	output.cOut = out
	output.sign = fftForward
	// output.flags = flags
	// output.input = make([]float64, n)
	// output.ip = make([]int, n)
	// output.w = make([]float64, n*5/4)

	// output.ip[0] = 0
	// makewt(output.n>>2, output.ip, output.w)
	// makect(output.n>>2, output.ip, output.w[output.n>>2:])
	return output
}

func fftExecute(p fftPlan) {
	if p.sign == fftForward {
		forwardFFT(p)
	} else {
		// ifft
		backwardFFT(p)
	}
}

// // The following functions are reffered by
// // http://www.kurims.kyoto-u.ac.jp/~ooura/index.html

// func cdft(n, isgn int, a []float64, ip []int, w []float64) {
// 	nw := ip[0]
// 	if isgn >= 0 {
// 		cftfsub(n, a, ip, nw, w)
// 	} else {
// 		cftbsub(n, a, ip, nw, w)
// 	}
// }

// func rdft(n, isgn int, a []float64, ip []int, w []float64) {
// 	nw := ip[0]
// 	nc := ip[1]

// 	if isgn >= 0 {
// 		if n > 4 {
// 			cftfsub(n, a, ip, nw, w)
// 			rftfsub(n, a, nc, w[nw:])
// 		} else if n == 4 {
// 			cftfsub(n, a, ip, nw, w)
// 		}
// 		xi := a[0] - a[1]
// 		a[0] += a[1]
// 		a[1] = xi
// 	} else {
// 		a[1] = 0.5 * (a[0] - a[1])
// 		a[0] -= a[1]
// 		if n > 4 {
// 			rftbsub(n, a, nc, w[nw:])
// 			cftbsub(n, a, ip, nw, w)
// 		} else if n == 4 {
// 			cftbsub(n, a, ip, nw, w)
// 		}
// 	}
// }

// func makewt(nw int, ip []int, w []float64) {
// 	ip[0] = nw
// 	ip[1] = 1
// 	if nw > 2 {
// 		nwh := nw >> 1
// 		delta := math.Atan(1.0) / float64(nwh)
// 		wn4r := math.Cos(delta * float64(nwh))
// 		w[0] = 1
// 		w[1] = wn4r
// 		if nwh == 4 {
// 			w[2] = math.Cos(delta * 2)
// 			w[3] = math.Sin(delta * 2)
// 		} else if nwh > 4 {
// 			makeipt(nw, ip)
// 			w[2] = 0.5 / math.Cos(delta*2)
// 			w[3] = 0.5 / math.Cos(delta*6)
// 			for j := 4; j < nwh; j += 4 {
// 				w[j] = math.Cos(delta * float64(j))
// 				w[j+1] = math.Sin(delta * float64(j))
// 				w[j+2] = math.Cos(3 * delta * float64(j))
// 				w[j+3] = -math.Sin(3 * delta * float64(j))
// 			}
// 		}
// 		nw0 := 0
// 		for nwh > 2 {
// 			nw1 := nw0 + nwh
// 			nwh >>= 1
// 			w[nw1] = 1
// 			w[nw1+1] = wn4r
// 			if nwh == 4 {
// 				wk1r := w[nw0+4]
// 				wk1i := w[nw0+5]
// 				w[nw1+2] = wk1r
// 				w[nw1+3] = wk1i
// 			} else if nwh > 4 {
// 				wk1r := w[nw0+4]
// 				wk3r := w[nw0+6]
// 				w[nw1+2] = 0.5 / wk1r
// 				w[nw1+3] = 0.5 / wk3r
// 				for j := 4; j < nwh; j += 4 {
// 					wk1r := w[nw0+2*j]
// 					wk1i := w[nw0+2*j+1]
// 					wk3r := w[nw0+2*j+2]
// 					wk3i := w[nw0+2*j+3]
// 					w[nw1+j] = wk1r
// 					w[nw1+j+1] = wk1i
// 					w[nw1+j+2] = wk3r
// 					w[nw1+j+3] = wk3i
// 				}
// 			}
// 			nw0 = nw1
// 		}
// 	}
// }

// func makeipt(nw int, ip []int) {
// 	ip[2] = 0
// 	ip[3] = 16
// 	m := 2
// 	for l := nw; l > 32; l >>= 2 {
// 		m2 := m << 1
// 		q := m2 << 3
// 		for j := m; j < m2; j++ {
// 			p := ip[j] << 2
// 			ip[m+j] = p
// 			ip[m2+j] = p + q
// 		}
// 		m = m2
// 	}
// }

// func makect(nc int, ip []int, c []float64) {
// 	ip[1] = nc
// 	if nc > 1 {
// 		nch := nc >> 1
// 		delta := math.Atan(1.0) / float64(nch)
// 		c[0] = math.Cos(delta * float64(nch))
// 		c[nch] = 0.5 * c[0]
// 		for j := 1; j < nch; j++ {
// 			c[j] = 0.5 * math.Cos(delta*float64(j))
// 			c[nc-j] = 0.5 * math.Sin(delta*float64(j))
// 		}
// 	}
// }
