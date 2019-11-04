package world

import "math"

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

// designLowCutFilter calculates the coefficients the filter.
func designLowCutFilter(n, fftSize int, lowCutFilter []float64) {
	for i := 1; i <= n; i++ {
		lowCutFilter[i-1] = 0.5 - 0.5*math.Cos(float64(i)*2.0*math.Pi/float64(n+1))
	}
	for i := n; i < fftSize; i++ {
		lowCutFilter[i] = 0.0
	}
	sumOfAmplitude := 0.0
	for i := 0; i < n; i++ {
		sumOfAmplitude += lowCutFilter[i]
	}
	for i := 0; i < n; i++ {
		lowCutFilter[i] = -lowCutFilter[i] / sumOfAmplitude
	}
	for i := 0; i < (n-1)/2; i++ {
		lowCutFilter[fftSize-(n-1)/2+i] = lowCutFilter[i]
	}
	for i := 0; i < n; i++ {
		lowCutFilter[i] = lowCutFilter[i+(n-1)/2]
	}
	lowCutFilter[0] += 1.0
}

// getSamplesForDIO calculates the number of samples required for Dio().
//
// Input:
//   fs             : Sampling frequency [Hz]
//   xLength        : Length of the input signal [Sample].
//   framePeriod    : Frame shift [msec]
//
// Output:
//   The number of samples required to store the results of Dio()
func getSamplesForDIO(fs, xLength int, framePeriod float64) int {
	return int(1000.0*float64(xLength)/float64(fs)/framePeriod) + 1
}

// getSpectrumForEstimation calculates the spectrum for estimation.
// This function carries out downsampling to speed up the estimation process
// and calculates the spectrum of the downsampled signal.
func getSpectrumForEstimation(x []float64, yLength int, actualFS float64, fftSize, decimationRatio int, ySpectrum []complex128) {
	xLength := len(x)
	y := make([]float64, fftSize)

	// Initialization
	for i := 0; i < fftSize; i++ {
		y[i] = 0.0
	}

	// Downsampling
	if decimationRatio != 1 {
		decimate(x, decimationRatio, y)
	} else {
		copy(y[:xLength], x)
	}

	// Removal of the DC component (y = y - mean value of y)
	meanY := 0.0
	for i := 0; i < yLength; i++ {
		meanY += y[i]
	}
	meanY /= float64(yLength)
	for i := 0; i < yLength; i++ {
		y[i] -= meanY
	}
	for i := yLength; i < fftSize; i++ {
		y[i] = 0.0
	}

	forwardFFT := fftPlanDftR2c1d(fftSize, y, ySpectrum, fftEstimate)
	fftExecute(forwardFFT)

	// Low cut filtering (from 0.1.4). Cut off frequency is 50.0 Hz.
	cutoffInSample := matlabRound(actualFS / cutOff)
	designLowCutFilter(cutoffInSample*2+1, fftSize, y)

	filterSpectrum := make([]complex128, fftSize)
	forwardFFT.cOut = filterSpectrum
	fftExecute(forwardFFT)

	for i := 0; i <= fftSize/2; i++ {
		ySpectrum[i] *= filterSpectrum[i]
	}
}

// getF0CandidateFromRawEvent() calculates F0 candidate contour in 1-ch signal
func getF0CandidateFromRawEvent(boundaryF0, fs float64,
	ySpectrum []complex128, yLength, fftSize int, f0Floor, f0Ceil float64,
	temporalPositions []float64, f0Length int,
	f0Score, f0Candidate []float64) {
	filteredSignal := make([]float64, fftSize)
	getFilteredSignal(matlabRound(fs/boundaryF0/2.0), fftSize, ySpectrum,
		yLength, filteredSignal)

	var zeroCrossings zeroCrossings
	getFourZeroCrossingIntervals(filteredSignal, yLength, fs,
		&zeroCrossings)

	getF0CandidateContour(&zeroCrossings, boundaryF0, f0Floor, f0Ceil,
		temporalPositions, f0Length, f0Candidate, f0Score)
}

// getF0CandidatesAndScores calculates all f0 candidates and their scores.
func getF0CandidatesAndScores(boundaryF0List []float64,
	numberOfBands int, actualFS float64, yLength int,
	temporalPositions []float64, f0Length int,
	ySpectrum []complex128, fftSize int, f0Floor, f0Ceil float64,
	rawF0Candidates, rawF0Scores [][]float64) {
	f0Candidate := make([]float64, f0Length)
	f0Score := make([]float64, f0Length)

	// Calculation of the acoustics events (zero-crossing)
	for i := 0; i < numberOfBands; i++ {
		getF0CandidateFromRawEvent(boundaryF0List[i], actualFS, ySpectrum,
			yLength, fftSize, f0Floor, f0Ceil, temporalPositions, f0Length,
			f0Score, f0Candidate)
		for j := 0; j < f0Length; j++ {
			// A way to avoid zero division
			rawF0Scores[i][j] = f0Score[j] /
				(f0Candidate[j] + mySafeGuardMinimum)
			rawF0Candidates[i][j] = f0Candidate[j]
		}
	}
}

// getBestF0Contour calculates the best f0 contour based on scores of
// all candidates. The F0 with highest score is selected.
func getBestF0Contour(f0Length int, f0Candidates, f0Scores [][]float64,
	numberOfBands int, bestF0Contour []float64) {
	for i := 0; i < f0Length; i++ {
		tmp := f0Scores[0][i]
		bestF0Contour[i] = f0Candidates[0][i]
		for j := 1; j < numberOfBands; j++ {
			if tmp > f0Scores[j][i] {
				tmp = f0Scores[j][i]
				bestF0Contour[i] = f0Candidates[j][i]
			}
		}
	}
}

// fixStep1 is the 1st step of the postprocessing.
// This function eliminates the unnatural change of f0 based on allowedRange.
func fixStep1(bestF0Contour []float64, f0Length,
	voiceRangeMinimum int, allowedRange float64, f0Step1 []float64) {
	f0Base := make([]float64, f0Length)
	// Initialization
	for i := 0; i < voiceRangeMinimum; i++ {
		f0Base[i] = 0.0
	}
	for i := voiceRangeMinimum; i < f0Length-voiceRangeMinimum; i++ {
		f0Base[i] = bestF0Contour[i]
	}
	for i := f0Length - voiceRangeMinimum; i < f0Length; i++ {
		f0Base[i] = 0.0
	}

	// Processing to prevent the jumping of f0
	for i := 0; i < voiceRangeMinimum; i++ {
		f0Step1[i] = 0.0
	}
	for i := voiceRangeMinimum; i < f0Length; i++ {
		if math.Abs((f0Base[i]-f0Base[i-1])/
			(mySafeGuardMinimum+f0Base[i])) <
			allowedRange {
			f0Step1[i] = f0Base[i]
		} else {
			f0Step1[i] = 0.0

		}
	}
}

// fixStep2 is the 2nd step of the postprocessing.
// This function eliminates the suspected f0 in the anlaut and auslaut.
func fixStep2(f0Step1 []float64, f0Length,
	voiceRangeMinimum int, f0Step2 []float64) {
	for i := 0; i < f0Length; i++ {
		f0Step2[i] = f0Step1[i]
	}

	center := (voiceRangeMinimum - 1) / 2
	for i := center; i < f0Length-center; i++ {
		for j := -center; j <= center; j++ {
			if f0Step1[i+j] == 0 {
				f0Step2[i] = 0.0
				break
			}
		}
	}
}

// getNumberOfVoicedSections() counts the number of voiced sections.
func getNumberOfVoicedSections(f0 []float64, f0Length int,
	positiveIndex, negativeIndex []int) (int, int) {
	positiveCount := 0
	negativeCount := 0
	for i := 1; i < f0Length; i++ {
		if f0[i] == 0 && f0[i-1] != 0 {
			negativeIndex[negativeCount] = i - 1
			negativeCount++
		} else {
			if f0[i-1] == 0 && f0[i] != 0 {
				positiveIndex[positiveCount] = i
				positiveCount++
			}
		}
	}
	return positiveCount, negativeCount
}

// selectBestF0() corrects the f0[currentIndex] based on
// f0[currentIndex + sign].
func selectBestF0(currentF0, pastF0 float64,
	f0Candidates [][]float64, numberOfCandidates,
	targetIndex int, allowedRange float64) float64 {
	referenceF0 := (currentF0*3.0 - pastF0) / 2.0

	minimumError := math.Abs(referenceF0 - f0Candidates[0][targetIndex])
	bestF0 := f0Candidates[0][targetIndex]

	var currentError float64
	for i := 1; i < numberOfCandidates; i++ {
		currentError = math.Abs(referenceF0 - f0Candidates[i][targetIndex])
		if currentError < minimumError {
			minimumError = currentError
			bestF0 = f0Candidates[i][targetIndex]
		}
	}
	if math.Abs(1.0-bestF0/referenceF0) > allowedRange {
		return 0.0
	}
	return bestF0
}

// fixStep3() is the 3rd step of the postprocessing.
// This function corrects the f0 candidates from backward to forward.
func fixStep3(f0Step2 []float64, f0Length int,
	f0Candidates [][]float64, numberOfCandidates int,
	allowedRange float64, negativeIndex []int, negativeCount int,
	f0Step3 []float64) {
	for i := 0; i < f0Length; i++ {
		f0Step3[i] = f0Step2[i]
	}

	for i := 0; i < negativeCount; i++ {
		var limit int
		if i == negativeCount-1 {
			limit = f0Length - 1
		} else {
			limit = negativeIndex[i+1]
		}
		for j := negativeIndex[i]; j < limit; j++ {
			f0Step3[j+1] =
				selectBestF0(f0Step3[j], f0Step3[j-1], f0Candidates,
					numberOfCandidates, j+1, allowedRange)
			if f0Step3[j+1] == 0 {
				break
			}
		}
	}
}

// fixStep4() is the 4th step of the postprocessing.
// This function corrects the f0 candidates from forward to backward.
func fixStep4(f0Step3 []float64, f0Length int,
	f0Candidates [][]float64, numberOfCandidates int,
	allowedRange float64, positiveIndex []int, positiveCount int,
	f0Step4 []float64) {
	for i := 0; i < f0Length; i++ {
		f0Step4[i] = f0Step3[i]
	}

	for i := positiveCount - 1; i >= 0; i-- {
		limit := 1
		if i != 0 {
			limit = positiveIndex[i-1]
		}
		for j := positiveIndex[i]; j > limit; j-- {
			f0Step4[j-1] =
				selectBestF0(f0Step4[j], f0Step4[j+1], f0Candidates,
					numberOfCandidates, j-1, allowedRange)
			if f0Step4[j-1] == 0 {
				break
			}
		}
	}
}

// fixF0Contour() calculates the definitive f0 contour based on all f0
// candidates. There are four steps.
func fixF0Contour(framePeriod float64, numberOfCandidates,
	fs int, f0Candidates [][]float64,
	bestF0Contour []float64, f0Length int, f0Floor,
	allowedRange float64, fixedF0Contour []float64) {
	voiceRangeMinimum := int(0.5+1000.0/framePeriod/f0Floor)*2 + 1

	if f0Length <= voiceRangeMinimum {
		return
	}

	f0Tmp1 := make([]float64, f0Length)
	f0Tmp2 := make([]float64, f0Length)

	fixStep1(bestF0Contour, f0Length, voiceRangeMinimum, allowedRange, f0Tmp1)
	fixStep2(f0Tmp1, f0Length, voiceRangeMinimum, f0Tmp2)

	positiveIndex := make([]int, f0Length)
	negativeIndex := make([]int, f0Length)
	positiveCount, negativeCount := getNumberOfVoicedSections(f0Tmp2, f0Length, positiveIndex, negativeIndex)
	fixStep3(f0Tmp2, f0Length, f0Candidates, numberOfCandidates,
		allowedRange, negativeIndex, negativeCount, f0Tmp1)
	fixStep4(f0Tmp1, f0Length, f0Candidates, numberOfCandidates,
		allowedRange, positiveIndex, positiveCount, fixedF0Contour)
}

// getFilteredSignal calculates the signal that is the convolution of the
// input signal and low-pass filter.
// This function is only used in rawEventByDio()
func getFilteredSignal(halfAverageLength, fftSize int,
	ySpectrum []complex128, yLength int, filteredSignal []float64) {
	lpf := make([]float64, fftSize)
	// Nuttall window is used as a low-pass filter.
	// Cutoff frequency depends on the window length.
	nuttallWindow(lpf[:halfAverageLength*4])

	lpfSpectrum := make([]complex128, fftSize)
	forwardFFT := fftPlanDftR2c1d(fftSize, lpf, lpfSpectrum, fftEstimate)
	fftExecute(forwardFFT)

	// Convolution
	lpfSpectrum[0] *= ySpectrum[0]
	for i := 1; i <= fftSize/2; i++ {
		lpfSpectrum[i] *= ySpectrum[i]
		lpfSpectrum[fftSize-i-1] = lpfSpectrum[i]
	}

	inverseFFT := fftPlanDftC2r1d(fftSize, lpfSpectrum, filteredSignal, fftEstimate)
	fftExecute(inverseFFT)

	// Compensation of the delay.
	indexBias := halfAverageLength * 2
	for i := 0; i < yLength; i++ {
		filteredSignal[i] = filteredSignal[i+indexBias]
	}
}

// checkEvent returns 1, provided that the input value is over 1.
// This function is for RawEventByDio().
func checkEvent(x int) int {
	if 0 < x {
		return 1
	}
	return 0
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
func getFourZeroCrossingIntervals(filteredSignal []float64, yLength int,
	actualFS float64, zeroCrossings *zeroCrossings) {
	// xLength / 4 (old version) is fixed at 2013/07/14
	maximumNumber := yLength
	zeroCrossings.negativeIntervalLocations = make([]float64, maximumNumber)
	zeroCrossings.positiveIntervalLocations = make([]float64, maximumNumber)
	zeroCrossings.peakIntervalLocations = make([]float64, maximumNumber)
	zeroCrossings.dipIntervalLocations = make([]float64, maximumNumber)
	zeroCrossings.negativeIntervals = make([]float64, maximumNumber)
	zeroCrossings.positiveIntervals = make([]float64, maximumNumber)
	zeroCrossings.peakIntervals = make([]float64, maximumNumber)
	zeroCrossings.dipIntervals = make([]float64, maximumNumber)

	zeroCrossings.numberOfNegatives = zeroCrossingEngine(filteredSignal,
		yLength, actualFS, zeroCrossings.negativeIntervalLocations,
		zeroCrossings.negativeIntervals)

	for i := 0; i < yLength; i++ {
		filteredSignal[i] = -filteredSignal[i]
	}
	zeroCrossings.numberOfPositives = zeroCrossingEngine(filteredSignal,
		yLength, actualFS, zeroCrossings.positiveIntervalLocations,
		zeroCrossings.positiveIntervals)

	for i := 0; i < yLength-1; i++ {
		filteredSignal[i] = filteredSignal[i] - filteredSignal[i+1]
	}
	zeroCrossings.numberOfPeaks = zeroCrossingEngine(filteredSignal,
		yLength-1, actualFS, zeroCrossings.peakIntervalLocations,
		zeroCrossings.peakIntervals)

	for i := 0; i < yLength-1; i++ {
		filteredSignal[i] = -filteredSignal[i]
	}
	zeroCrossings.numberOfDips = zeroCrossingEngine(filteredSignal,
		yLength-1, actualFS, zeroCrossings.dipIntervalLocations,
		zeroCrossings.dipIntervals)
}

// getF0CandidateContourSub calculates the f0 candidates and deviations.
// This is the sub-function of getF0Candidates() and assumes the calculation.
func getF0CandidateContourSub(
	interpolatedF0Set [4][]float64, f0Length int, f0Floor float64,
	f0Ceil, boundaryF0 float64, f0Candidate, f0Score []float64) {
	for i := 0; i < f0Length; i++ {
		f0Candidate[i] = (interpolatedF0Set[0][i] +
			interpolatedF0Set[1][i] + interpolatedF0Set[2][i] +
			interpolatedF0Set[3][i]) / 4.0

		f0Score[i] = math.Sqrt(((interpolatedF0Set[0][i]-f0Candidate[i])*
			(interpolatedF0Set[0][i]-f0Candidate[i]) +
			(interpolatedF0Set[1][i]-f0Candidate[i])*
				(interpolatedF0Set[1][i]-f0Candidate[i]) +
			(interpolatedF0Set[2][i]-f0Candidate[i])*
				(interpolatedF0Set[2][i]-f0Candidate[i]) +
			(interpolatedF0Set[3][i]-f0Candidate[i])*
				(interpolatedF0Set[3][i]-f0Candidate[i])) / 3.0)

		if f0Candidate[i] > boundaryF0 || f0Candidate[i] < boundaryF0/2.0 ||
			f0Candidate[i] > f0Ceil || f0Candidate[i] < f0Floor {
			f0Candidate[i] = 0.0
			f0Score[i] = maximumValue
		}
	}
}

// getF0CandidateContour() calculates the F0 candidates based on the
// zero-crossings.
func getF0CandidateContour(zeroCrossings *zeroCrossings,
	boundaryF0, f0Floor, f0Ceil float64,
	temporalPositions []float64, f0Length int,
	f0Candidate, f0Score []float64) {
	if 0 == checkEvent(zeroCrossings.numberOfNegatives-2)*
		checkEvent(zeroCrossings.numberOfPositives-2)*
		checkEvent(zeroCrossings.numberOfPeaks-2)*
		checkEvent(zeroCrossings.numberOfDips-2) {
		for i := 0; i < f0Length; i++ {
			f0Score[i] = maximumValue
			f0Candidate[i] = 0.0
		}
		return
	}

	var interpolatedF0Set [4][]float64
	for i := 0; i < 4; i++ {
		interpolatedF0Set[i] = make([]float64, f0Length)
	}

	interp1(zeroCrossings.negativeIntervalLocations,
		zeroCrossings.negativeIntervals,
		zeroCrossings.numberOfNegatives,
		temporalPositions, f0Length, interpolatedF0Set[0])
	interp1(zeroCrossings.positiveIntervalLocations,
		zeroCrossings.positiveIntervals,
		zeroCrossings.numberOfPositives,
		temporalPositions, f0Length, interpolatedF0Set[1])
	interp1(zeroCrossings.peakIntervalLocations,
		zeroCrossings.peakIntervals, zeroCrossings.numberOfPeaks,
		temporalPositions, f0Length, interpolatedF0Set[2])
	interp1(zeroCrossings.dipIntervalLocations,
		zeroCrossings.dipIntervals, zeroCrossings.numberOfDips,
		temporalPositions, f0Length, interpolatedF0Set[3])

	getF0CandidateContourSub(interpolatedF0Set, f0Length, f0Floor,
		f0Ceil, boundaryF0, f0Candidate, f0Score)
}

// dioGeneralBody estimates the F0 based on Distributed Inline-filter
// Operation.
func dioGeneralBody(x []float64, fs int,
	framePeriod, f0Floor, f0Ceil, channelsInOctave float64,
	speed int, allowedRange float64,
	temporalPositions, f0 []float64) {
	xLength := len(x)
	numberOfBands := 1 + int(math.Log2(f0Ceil/f0Floor)*channelsInOctave)
	boundaryF0List := make([]float64, numberOfBands)
	for i := 0; i < numberOfBands; i++ {
		boundaryF0List[i] = f0Floor * math.Pow(2.0, float64(i+1)/channelsInOctave)
	}

	// normalization
	decimationRatio := myMaxInt(myMinInt(speed, 12), 1)
	yLength := 1 + int(xLength/decimationRatio)
	actualFS := float64(fs) / float64(decimationRatio)
	fftSize := getSuitableFFTSize(yLength +
		matlabRound(actualFS/cutOff)*2 + 1 +
		4*int(1.0+actualFS/boundaryF0List[0]/2.0))

	// Calculation of the spectrum used for the f0 estimation
	ySpectrum := make([]complex128, fftSize)
	getSpectrumForEstimation(x, yLength, actualFS, fftSize,
		decimationRatio, ySpectrum)

	f0Candidates := make([][]float64, numberOfBands)
	f0Scores := make([][]float64, numberOfBands)
	f0Length := getSamplesForDIO(fs, xLength, framePeriod)
	for i := 0; i < numberOfBands; i++ {
		f0Candidates[i] = make([]float64, f0Length)
		f0Scores[i] = make([]float64, f0Length)
	}

	for i := 0; i < f0Length; i++ {
		temporalPositions[i] = float64(i) * framePeriod / 1000.0
	}

	getF0CandidatesAndScores(boundaryF0List, numberOfBands,
		actualFS, yLength, temporalPositions, f0Length, ySpectrum,
		fftSize, f0Floor, f0Ceil, f0Candidates, f0Scores)

	// Selection of the best value based on fundamental-ness.
	// This function is related with SortCandidates() in MATLAB.
	bestF0Contour := make([]float64, f0Length)
	getBestF0Contour(f0Length, f0Candidates, f0Scores,
		numberOfBands, bestF0Contour)

	// // Postprocessing to find the best f0-contour.
	fixF0Contour(framePeriod, numberOfBands, fs, f0Candidates,
		bestF0Contour, f0Length, f0Floor, allowedRange, f0)
}

// DioOption is the struct to order the parameter for Dio().
type DioOption struct {
	F0Floor          float64
	F0Ceil           float64
	ChannelsInOctave float64
	FramePeriod      float64 // msec
	Speed            int     // (1, 2, ..., 12)
	AllowedRange     float64 // Threshold used for fixing the F0 contour.
}

// GetSamplesForDIO calculates the number of samples required for Dio().
//
// Input:
//   fs             : Sampling frequency [Hz]
//   xLength        : Length of the input signal [Sample].
//   framePeriod    : Frame shift [msec]
//
// Output:
//   The number of samples required to store the results of Dio()
func GetSamplesForDIO(fs, xLength int, framePeriod float64) int {
	return int(1000.0*float64(xLength)/float64(fs)/framePeriod) + 1
}

// Dio estimates the F0.
//
// Input:
//   x                    : Input signal
//   xLength              : Length of x
//   fs                   : Sampling frequency
//   option               : Struct to order the parameter for DIO
//
// Output:
//   temporal_positions   : Temporal positions.
//   f0                   : F0 contour.
func Dio(x []float64, fs int, option *DioOption, temporalPositions, f0 []float64) {
	dioGeneralBody(x, fs, option.FramePeriod, option.F0Floor,
		option.F0Ceil, option.ChannelsInOctave, option.Speed,
		option.AllowedRange, temporalPositions, f0)
}

// InitializeDioOption allocates the memory to the struct and sets the
// default parameters.
//
// Output:
//   option   : Struct for the optional parameter.
func InitializeDioOption(option *DioOption) {
	// You can change default parameters.
	option.ChannelsInOctave = 2.0
	option.F0Ceil = ceilF0
	option.F0Floor = floorF0
	option.FramePeriod = 5

	// You can use the value from 1 to 12.
	// Default value 11 is for the fs of 44.1 kHz.
	// The lower value you use, the better performance you can obtain.
	option.Speed = 1

	// You can give a positive real number as the threshold.
	// The most strict value is 0, and there is no upper limit.
	// On the other hand, I think that the value from 0.02 to 0.2 is reasonable.
	option.AllowedRange = 0.1
}
