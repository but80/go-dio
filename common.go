package world

import (
	"math"
)

func myMaxInt(x, y int) int {
	if x < y {
		return y
	}
	return x
}

func myMinInt(x, y int) int {
	if x < y {
		return x
	}
	return y
}

func getSuitableFFTSize(sample int) int {
	e := uint(math.Log2(float64(sample)))
	return 1 << (e + 1)
}

func nuttallWindow(y []float64) {
	for i := range y {
		tmp := float64(i) / float64(len(y)-1)
		y[i] = 0.355768 - 0.487396*math.Cos(2.0*math.Pi*tmp) +
			0.144232*math.Cos(4.0*math.Pi*tmp) -
			0.012604*math.Cos(6.0*math.Pi*tmp)
	}
}
