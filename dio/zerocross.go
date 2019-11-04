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

// getFourZeroCrossingIntervals() calculates four zero-crossing intervals.
// (1) Zero-crossing going from negative to positive.
// (2) Zero-crossing going from positive to negative.
// (3) Peak, and (4) dip. (3) and (4) are calculated from the zero-crossings of
// the differential of waveform.
func (s *Session) getFourZeroCrossingIntervals(filteredSignal []float64, zeroCrossings *zeroCrossings) {
	// xLength / 4 (old version) is fixed at 2013/07/14
	maximumNumber := s.yLength
	zeroCrossings.negativeIntervalLocations = make([]float64, maximumNumber)
	zeroCrossings.positiveIntervalLocations = make([]float64, maximumNumber)
	zeroCrossings.peakIntervalLocations = make([]float64, maximumNumber)
	zeroCrossings.dipIntervalLocations = make([]float64, maximumNumber)
	zeroCrossings.negativeIntervals = make([]float64, maximumNumber)
	zeroCrossings.positiveIntervals = make([]float64, maximumNumber)
	zeroCrossings.peakIntervals = make([]float64, maximumNumber)
	zeroCrossings.dipIntervals = make([]float64, maximumNumber)

	zeroCrossings.numberOfNegatives = zeroCrossingEngine(filteredSignal,
		s.yLength, s.actualFS, zeroCrossings.negativeIntervalLocations,
		zeroCrossings.negativeIntervals)

	for i := 0; i < s.yLength; i++ {
		filteredSignal[i] = -filteredSignal[i]
	}
	zeroCrossings.numberOfPositives = zeroCrossingEngine(filteredSignal,
		s.yLength, s.actualFS, zeroCrossings.positiveIntervalLocations,
		zeroCrossings.positiveIntervals)

	for i := 0; i < s.yLength-1; i++ {
		filteredSignal[i] = filteredSignal[i] - filteredSignal[i+1]
	}
	zeroCrossings.numberOfPeaks = zeroCrossingEngine(filteredSignal,
		s.yLength-1, s.actualFS, zeroCrossings.peakIntervalLocations,
		zeroCrossings.peakIntervals)

	for i := 0; i < s.yLength-1; i++ {
		filteredSignal[i] = -filteredSignal[i]
	}
	zeroCrossings.numberOfDips = zeroCrossingEngine(filteredSignal,
		s.yLength-1, s.actualFS, zeroCrossings.dipIntervalLocations,
		zeroCrossings.dipIntervals)
}
