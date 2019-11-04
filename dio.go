package world

import (
	"math"

	"gonum.org/v1/gonum/fourier"
)

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

// getSpectrumForEstimation calculates the spectrum for estimation.
// This function carries out downsampling to speed up the estimation process
// and calculates the spectrum of the downsampled signal.
func (s *DioSession) getSpectrumForEstimation() {
	y := make([]float64, s.fftSize)

	// Downsampling
	if s.decimationRatio != 1 {
		decimate(s.x, s.decimationRatio, y)
	} else {
		copy(y[:len(s.x)], s.x)
	}

	// Removal of the DC component (y = y - mean value of y)
	meanY := 0.0
	for i := 0; i < s.yLength; i++ {
		meanY += y[i]
	}
	meanY /= float64(s.yLength)
	for i := 0; i < s.yLength; i++ {
		y[i] -= meanY
	}
	for i := s.yLength; i < s.fftSize; i++ {
		y[i] = 0.0
	}

	s.fft.Coefficients(s.ySpectrum[:s.fftSize/2+1], y)

	// Low cut filtering (from 0.1.4). Cut off frequency is 50.0 Hz.
	cutoffInSample := matlabRound(s.actualFS / cutOff)
	designLowCutFilter(cutoffInSample*2+1, s.fftSize, y)

	filterSpectrum := make([]complex128, s.fftSize/2+1)
	s.fft.Coefficients(filterSpectrum, y)

	for i := 0; i <= s.fftSize/2; i++ {
		s.ySpectrum[i] *= filterSpectrum[i]
	}
}

// getF0CandidateFromRawEvent() calculates F0 candidate contour in 1-ch signal
func (s *DioSession) getF0CandidateFromRawEvent(boundaryF0 float64) {
	filteredSignal := make([]float64, s.fftSize)
	s.getFilteredSignal(matlabRound(s.actualFS/boundaryF0/2.0), filteredSignal)

	var zeroCrossings zeroCrossings
	s.getFourZeroCrossingIntervals(filteredSignal, &zeroCrossings)

	s.getF0CandidateContour(&zeroCrossings, boundaryF0)
}

// getF0CandidatesAndScores calculates all f0 candidates and their scores.
func (s *DioSession) getF0CandidatesAndScores() {
	// Calculation of the acoustics events (zero-crossing)
	for i := 0; i < s.numberOfBands; i++ {
		s.getF0CandidateFromRawEvent(s.boundaryF0List[i])
		for j := 0; j < s.f0Length; j++ {
			// A way to avoid zero division
			s.f0Scores[i][j] = s.f0Score[j] / (s.f0Candidate[j] + mySafeGuardMinimum)
			s.f0Candidates[i][j] = s.f0Candidate[j]
		}
	}
}

// getBestF0Contour calculates the best f0 contour based on scores of
// all candidates. The F0 with highest score is selected.
func (s *DioSession) getBestF0Contour(bestF0Contour []float64) {
	for i := 0; i < s.f0Length; i++ {
		tmp := s.f0Scores[0][i]
		bestF0Contour[i] = s.f0Candidates[0][i]
		for j := 1; j < s.numberOfBands; j++ {
			if tmp > s.f0Scores[j][i] {
				tmp = s.f0Scores[j][i]
				bestF0Contour[i] = s.f0Candidates[j][i]
			}
		}
	}
}

// fixStep1 is the 1st step of the postprocessing.
// This function eliminates the unnatural change of f0 based on allowedRange.
func (s *DioSession) fixStep1(bestF0Contour []float64, voiceRangeMinimum int, f0Step1 []float64) {
	f0Base := make([]float64, s.f0Length)
	// Initialization
	for i := 0; i < voiceRangeMinimum; i++ {
		f0Base[i] = 0.0
	}
	for i := voiceRangeMinimum; i < s.f0Length-voiceRangeMinimum; i++ {
		f0Base[i] = bestF0Contour[i]
	}
	for i := s.f0Length - voiceRangeMinimum; i < s.f0Length; i++ {
		f0Base[i] = 0.0
	}

	// Processing to prevent the jumping of f0
	for i := 0; i < voiceRangeMinimum; i++ {
		f0Step1[i] = 0.0
	}
	for i := voiceRangeMinimum; i < s.f0Length; i++ {
		if math.Abs((f0Base[i]-f0Base[i-1])/
			(mySafeGuardMinimum+f0Base[i])) <
			s.option.AllowedRange {
			f0Step1[i] = f0Base[i]
		} else {
			f0Step1[i] = 0.0
		}
	}
}

// fixStep2 is the 2nd step of the postprocessing.
// This function eliminates the suspected f0 in the anlaut and auslaut.
func (s *DioSession) fixStep2(f0Step1 []float64, voiceRangeMinimum int, f0Step2 []float64) {
	for i := 0; i < s.f0Length; i++ {
		f0Step2[i] = f0Step1[i]
	}

	center := (voiceRangeMinimum - 1) / 2
	for i := center; i < s.f0Length-center; i++ {
		for j := -center; j <= center; j++ {
			if f0Step1[i+j] == 0 {
				f0Step2[i] = 0.0
				break
			}
		}
	}
}

// getNumberOfVoicedSections() counts the number of voiced sections.
func (s *DioSession) getNumberOfVoicedSections(f0 []float64, positiveIndex, negativeIndex []int) (int, int) {
	positiveCount := 0
	negativeCount := 0
	for i := 1; i < s.f0Length; i++ {
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
func (s *DioSession) selectBestF0(currentF0, pastF0 float64, targetIndex int) float64 {
	referenceF0 := (currentF0*3.0 - pastF0) / 2.0

	minimumError := math.Abs(referenceF0 - s.f0Candidates[0][targetIndex])
	bestF0 := s.f0Candidates[0][targetIndex]

	var currentError float64
	for i := 1; i < s.numberOfBands; i++ {
		currentError = math.Abs(referenceF0 - s.f0Candidates[i][targetIndex])
		if currentError < minimumError {
			minimumError = currentError
			bestF0 = s.f0Candidates[i][targetIndex]
		}
	}
	if math.Abs(1.0-bestF0/referenceF0) > s.option.AllowedRange {
		return 0.0
	}
	return bestF0
}

// fixStep3() is the 3rd step of the postprocessing.
// This function corrects the f0 candidates from backward to forward.
func (s *DioSession) fixStep3(f0Step2 []float64, negativeIndex []int, negativeCount int, f0Step3 []float64) {
	for i := 0; i < s.f0Length; i++ {
		f0Step3[i] = f0Step2[i]
	}

	for i := 0; i < negativeCount; i++ {
		var limit int
		if i == negativeCount-1 {
			limit = s.f0Length - 1
		} else {
			limit = negativeIndex[i+1]
		}
		for j := negativeIndex[i]; j < limit; j++ {
			f0Step3[j+1] = s.selectBestF0(f0Step3[j], f0Step3[j-1], j+1)
			if f0Step3[j+1] == 0 {
				break
			}
		}
	}
}

// fixStep4() is the 4th step of the postprocessing.
// This function corrects the f0 candidates from forward to backward.
func (s *DioSession) fixStep4(f0Step3 []float64,
	positiveIndex []int, positiveCount int) {
	f0Step4 := s.f0

	for i := 0; i < s.f0Length; i++ {
		f0Step4[i] = f0Step3[i]
	}

	for i := positiveCount - 1; i >= 0; i-- {
		limit := 1
		if i != 0 {
			limit = positiveIndex[i-1]
		}
		for j := positiveIndex[i]; j > limit; j-- {
			f0Step4[j-1] = s.selectBestF0(f0Step4[j], f0Step4[j+1], j-1)
			if f0Step4[j-1] == 0 {
				break
			}
		}
	}
}

// fixF0Contour() calculates the definitive f0 contour based on all f0
// candidates. There are four steps.
func (s *DioSession) fixF0Contour(bestF0Contour []float64) {
	o := s.option
	voiceRangeMinimum := int(0.5+1000.0/o.FramePeriod/o.F0Floor)*2 + 1

	if s.f0Length <= voiceRangeMinimum {
		return
	}

	f0Tmp1 := make([]float64, s.f0Length)
	f0Tmp2 := make([]float64, s.f0Length)

	s.fixStep1(bestF0Contour, voiceRangeMinimum, f0Tmp1)
	s.fixStep2(f0Tmp1, voiceRangeMinimum, f0Tmp2)

	positiveIndex := make([]int, s.f0Length)
	negativeIndex := make([]int, s.f0Length)
	positiveCount, negativeCount := s.getNumberOfVoicedSections(f0Tmp2, positiveIndex, negativeIndex)
	s.fixStep3(f0Tmp2, negativeIndex, negativeCount, f0Tmp1)
	s.fixStep4(f0Tmp1, positiveIndex, positiveCount)
}

// getFilteredSignal calculates the signal that is the convolution of the
// input signal and low-pass filter.
// This function is only used in rawEventByDio()
func (s *DioSession) getFilteredSignal(halfAverageLength int, filteredSignal []float64) {
	lpf := make([]float64, s.fftSize)
	// Nuttall window is used as a low-pass filter.
	// Cutoff frequency depends on the window length.
	nuttallWindow(lpf[:halfAverageLength*4])

	lpfSpectrum := make([]complex128, s.fftSize/2+1)
	s.fft.Coefficients(lpfSpectrum, lpf)

	// Convolution
	lpfSpectrum[0] *= s.ySpectrum[0]
	for i := 1; i <= s.fftSize/2; i++ {
		lpfSpectrum[i] *= s.ySpectrum[i]
	}

	s.fft.Sequence(filteredSignal, lpfSpectrum)

	// Compensation of the delay.
	indexBias := halfAverageLength * 2
	for i := 0; i < s.yLength; i++ {
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
func (s *DioSession) getFourZeroCrossingIntervals(filteredSignal []float64, zeroCrossings *zeroCrossings) {
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

// getF0CandidateContourSub calculates the f0 candidates and deviations.
// This is the sub-function of getF0Candidates() and assumes the calculation.
func (s *DioSession) getF0CandidateContourSub(interpolatedF0Set [4][]float64, boundaryF0 float64) {
	for i := 0; i < s.f0Length; i++ {
		s.f0Candidate[i] = (interpolatedF0Set[0][i] +
			interpolatedF0Set[1][i] + interpolatedF0Set[2][i] +
			interpolatedF0Set[3][i]) / 4.0

		s.f0Score[i] = math.Sqrt(((interpolatedF0Set[0][i]-s.f0Candidate[i])*(interpolatedF0Set[0][i]-s.f0Candidate[i]) +
			(interpolatedF0Set[1][i]-s.f0Candidate[i])*(interpolatedF0Set[1][i]-s.f0Candidate[i]) +
			(interpolatedF0Set[2][i]-s.f0Candidate[i])*(interpolatedF0Set[2][i]-s.f0Candidate[i]) +
			(interpolatedF0Set[3][i]-s.f0Candidate[i])*(interpolatedF0Set[3][i]-s.f0Candidate[i])) / 3.0)

		if s.f0Candidate[i] > boundaryF0 || s.f0Candidate[i] < boundaryF0/2.0 ||
			s.f0Candidate[i] > s.option.F0Ceil || s.f0Candidate[i] < s.option.F0Floor {
			s.f0Candidate[i] = 0.0
			s.f0Score[i] = maximumValue
		}
	}
}

// getF0CandidateContour() calculates the F0 candidates based on the
// zero-crossings.
func (s *DioSession) getF0CandidateContour(zeroCrossings *zeroCrossings, boundaryF0 float64) {
	if 0 == checkEvent(zeroCrossings.numberOfNegatives-2)*
		checkEvent(zeroCrossings.numberOfPositives-2)*
		checkEvent(zeroCrossings.numberOfPeaks-2)*
		checkEvent(zeroCrossings.numberOfDips-2) {
		for i := 0; i < s.f0Length; i++ {
			s.f0Score[i] = maximumValue
			s.f0Candidate[i] = 0.0
		}
		return
	}

	var interpolatedF0Set [4][]float64
	for i := 0; i < 4; i++ {
		interpolatedF0Set[i] = make([]float64, s.f0Length)
	}

	interp1(zeroCrossings.negativeIntervalLocations[:zeroCrossings.numberOfNegatives],
		zeroCrossings.negativeIntervals,
		s.temporalPositions, interpolatedF0Set[0])
	interp1(zeroCrossings.positiveIntervalLocations[:zeroCrossings.numberOfPositives],
		zeroCrossings.positiveIntervals,
		s.temporalPositions, interpolatedF0Set[1])
	interp1(zeroCrossings.peakIntervalLocations[:zeroCrossings.numberOfPeaks],
		zeroCrossings.peakIntervals,
		s.temporalPositions, interpolatedF0Set[2])
	interp1(zeroCrossings.dipIntervalLocations[:zeroCrossings.numberOfDips],
		zeroCrossings.dipIntervals,
		s.temporalPositions, interpolatedF0Set[3])

	s.getF0CandidateContourSub(interpolatedF0Set, boundaryF0)
}

// generalBody estimates the F0 based on Distributed Inline-filter Operation.
func (s *DioSession) generalBody() {
	// Calculation of the spectrum used for the f0 estimation
	s.getSpectrumForEstimation()

	s.getF0CandidatesAndScores()

	// Selection of the best value based on fundamental-ness.
	// This function is related with SortCandidates() in MATLAB.
	bestF0Contour := make([]float64, s.f0Length)
	s.getBestF0Contour(bestF0Contour)

	// Postprocessing to find the best f0-contour.
	s.fixF0Contour(bestF0Contour)
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

type DioSession struct {
	// Input
	x      []float64  // Input signal
	fs     int        // Sampling frequency
	option *DioOption // Struct to order the parameter for DIO

	// Temporary
	f0Length        int
	yLength         int
	decimationRatio int
	actualFS        float64
	fftSize         int
	fft             *fourier.FFT
	ySpectrum       []complex128
	numberOfBands   int
	boundaryF0List  []float64
	f0Candidates    [][]float64
	f0Scores        [][]float64
	f0Candidate     []float64
	f0Score         []float64

	// Output
	temporalPositions []float64 // Temporal positions.
	f0                []float64 // F0 contour.
}

func NewDioSession(x []float64, fs int, option *DioOption) *DioSession {
	numberOfBands := 1 + int(math.Log2(option.F0Ceil/option.F0Floor)*option.ChannelsInOctave)
	f0Length := getSamplesForDIO(fs, len(x), option.FramePeriod)
	temporalPositions := make([]float64, f0Length)
	f0 := make([]float64, f0Length)

	boundaryF0List := make([]float64, numberOfBands)
	for i := 0; i < numberOfBands; i++ {
		boundaryF0List[i] = option.F0Floor * math.Pow(2.0, float64(i+1)/option.ChannelsInOctave)
	}

	s := &DioSession{
		x:      x,
		fs:     fs,
		option: option,

		f0Length:       f0Length,
		numberOfBands:  numberOfBands,
		boundaryF0List: boundaryF0List,
		f0Candidate:    make([]float64, f0Length),
		f0Score:        make([]float64, f0Length),

		temporalPositions: temporalPositions,
		f0:                f0,
	}

	// normalization
	s.decimationRatio = myMaxInt(myMinInt(option.Speed, 12), 1)
	s.yLength = 1 + len(s.x)/s.decimationRatio
	s.actualFS = float64(s.fs) / float64(s.decimationRatio)
	s.fftSize = getSuitableFFTSize(s.yLength +
		matlabRound(s.actualFS/cutOff)*2 + 1 +
		4*int(1.0+s.actualFS/s.boundaryF0List[0]/2.0))
	s.fft = fourier.NewFFT(s.fftSize)

	s.ySpectrum = make([]complex128, s.fftSize)

	s.f0Candidates = make([][]float64, s.numberOfBands)
	s.f0Scores = make([][]float64, s.numberOfBands)
	for i := 0; i < s.numberOfBands; i++ {
		s.f0Candidates[i] = make([]float64, s.f0Length)
		s.f0Scores[i] = make([]float64, s.f0Length)
	}

	for i := range s.temporalPositions {
		s.temporalPositions[i] = float64(i) * s.option.FramePeriod / 1000.0
	}

	return s
}

// Run estimates the F0.
func (s *DioSession) Run() ([]float64, []float64) {
	s.generalBody()
	return s.temporalPositions, s.f0
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
