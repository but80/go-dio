package dio

// struct for getFourZeroCrossingIntervals()
// "negative" means "zero-crossing point going from positive to negative"
// "positive" means "zero-crossing point going from negative to positive"
type zeroCrossings struct {
	negativeIntervalLocations []float64
	negativeIntervals         []float64
	numberOfNegatives         int
	positiveIntervalLocations []float64
	positiveIntervals         []float64
	numberOfPositives         int
	peakIntervalLocations     []float64
	peakIntervals             []float64
	numberOfPeaks             int
	dipIntervalLocations      []float64
	dipIntervals              []float64
	numberOfDips              int
}

func newZeroCrossings(n int) *zeroCrossings {
	return &zeroCrossings{
		negativeIntervalLocations: make([]float64, n),
		positiveIntervalLocations: make([]float64, n),
		peakIntervalLocations:     make([]float64, n),
		dipIntervalLocations:      make([]float64, n),
		negativeIntervals:         make([]float64, n),
		positiveIntervals:         make([]float64, n),
		peakIntervals:             make([]float64, n),
		dipIntervals:              make([]float64, n),
	}
}

// zeroCrossingEngine calculates the zero crossing points from positive to
// negative. Thanks to Custom.Maid http://custom-made.seesaa.net/ (2012/8/19)
func zeroCrossingEngine(filteredSignal []float64, yLength int,
	fs float64, intervalLocations, intervals []float64) int {
	negativeGoingPoints := make([]int, yLength)

	for i := 0; i < yLength-1; i++ {
		if 0.0 < filteredSignal[i] && filteredSignal[i+1] <= 0.0 {
			negativeGoingPoints[i] = i + 1
		}
	}
	negativeGoingPoints[yLength-1] = 0

	edges := make([]int, yLength)
	count := 0
	for i := 0; i < yLength; i++ {
		if negativeGoingPoints[i] > 0 {
			edges[count] = negativeGoingPoints[i]
			count++
		}
	}

	if count < 2 {
		return 0
	}

	fineEdges := make([]float64, count)
	for i := 0; i < count; i++ {
		fineEdges[i] =
			float64(edges[i]) - filteredSignal[edges[i]-1]/
				(filteredSignal[edges[i]]-filteredSignal[edges[i]-1])
	}

	for i := 0; i < count-1; i++ {
		intervals[i] = fs / (fineEdges[i+1] - fineEdges[i])
		intervalLocations[i] = (fineEdges[i] + fineEdges[i+1]) / 2.0 / fs
	}

	return count - 1
}
