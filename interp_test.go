package dio

import "fmt"

func ExampleInterp1() {
	xy := []interval{
		{interval: 10, location: 0},
		{interval: 20, location: 2},
		{interval: 0, location: 6},
	}
	xi := []float64{0, 1, 2, 3, 4, 5, 6}
	yi := []float64{0, 0, 0, 0, 0, 0, 0}

	interp1(xy, xi, yi)
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
