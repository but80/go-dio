package world

import (
	"fmt"
)

func ExampleFFTPlan() {
	const fftSize = 8
	y := []float64{1, -1, .5, -.5, .25, -.25, .125, 0}
	ySpectrum := make([]complex128, fftSize)
	forwardFFT := fftPlanDftR2c1d(fftSize, y, ySpectrum, 0)
	fftExecute(forwardFFT)
	for _, v := range ySpectrum {
		fmt.Printf("%+4.3f\n", v)
	}
	fmt.Println()

	y2 := make([]float64, fftSize)
	backwardFFT := fftPlanDftC2r1d(fftSize, ySpectrum, y2, 0)
	fftExecute(backwardFFT)
	for _, v := range y2 {
		fmt.Printf("%+4.3f\n", v)
	}
	fmt.Println()

	// Output:
	// (+0.125+0.000i)
	// (+0.573+0.509i)
	// (+0.625+0.750i)
	// (+0.927+1.259i)
	// (+3.625+0.000i)
	// (+0.000+0.000i)
	// (+0.000+0.000i)
	// (+0.000+0.000i)
	//
	// +8.000
	// -8.000
	// +4.000
	// -4.000
	// +2.000
	// -2.000
	// +1.000
	// +0.000
}
