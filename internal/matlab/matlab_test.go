package matlab

import "fmt"

func ExampleInterp1() {
	x := []float64{0, 2, 6}
	y := []float64{10, 20, 0}
	xi := []float64{0, 1, 2, 3, 4, 5, 6}
	yi := []float64{0, 0, 0, 0, 0, 0, 0}

	Interp1(x, y, xi, yi)
	for _, v := range yi {
		fmt.Printf("%4.1f\n", v)
	}

	// Output:
	// 10.0
	// 15.0
	// 20.0
	// 15.0
	// 10.0
	//  5.0
	//  0.0
}
