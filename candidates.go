package dio

import (
	"math"

	"github.com/but80/go-dio/constant"
	"github.com/but80/go-dio/internal/common"
	"github.com/but80/go-dio/internal/matlab"
)

// checkEvent returns 1, provided that the input value is over 1.
// This function is for RawEventByDio().
func checkEvent(x int) int {
	if 0 < x {
		return 1
	}
	return 0
}

// getFilteredSignal calculates the signal that is the convolution of the
// input signal and low-pass filter.
// This function is only used in rawEventByDio()
func (s *Session) getFilteredSignal(halfAverageLength int, filteredSignal []float64) {
	lpf := make([]float64, s.fftSize)
	// Nuttall window is used as a low-pass filter.
	// Cutoff frequency depends on the window length.
	common.NuttallWindow(lpf[:halfAverageLength*4])

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

// getFourZeroCrossingIntervals() calculates four zero-crossing intervals.
// (1) Zero-crossing going from negative to positive.
// (2) Zero-crossing going from positive to negative.
// (3) Peak, and (4) dip. (3) and (4) are calculated from the zero-crossings of
// the differential of waveform.
func (s *Session) getFourZeroCrossingIntervals(filteredSignal []float64) {
	// xLength / 4 (old version) is fixed at 2013/07/14
	s.zeroCrossings.numberOfNegatives = zeroCrossingEngine(filteredSignal,
		s.yLength, s.actualFS, s.zeroCrossings.negativeIntervalLocations,
		s.zeroCrossings.negativeIntervals)

	for i := 0; i < s.yLength; i++ {
		filteredSignal[i] = -filteredSignal[i]
	}
	s.zeroCrossings.numberOfPositives = zeroCrossingEngine(filteredSignal,
		s.yLength, s.actualFS, s.zeroCrossings.positiveIntervalLocations,
		s.zeroCrossings.positiveIntervals)

	for i := 0; i < s.yLength-1; i++ {
		filteredSignal[i] = filteredSignal[i] - filteredSignal[i+1]
	}
	s.zeroCrossings.numberOfPeaks = zeroCrossingEngine(filteredSignal,
		s.yLength-1, s.actualFS, s.zeroCrossings.peakIntervalLocations,
		s.zeroCrossings.peakIntervals)

	for i := 0; i < s.yLength-1; i++ {
		filteredSignal[i] = -filteredSignal[i]
	}
	s.zeroCrossings.numberOfDips = zeroCrossingEngine(filteredSignal,
		s.yLength-1, s.actualFS, s.zeroCrossings.dipIntervalLocations,
		s.zeroCrossings.dipIntervals)
}

// getF0CandidateContourSub calculates the f0 candidates and deviations.
// This is the sub-function of getF0Candidates() and assumes the calculation.
func (s *Session) getF0CandidateContourSub(interpolatedF0Set [4][]float64, boundaryF0 float64) {
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
			s.f0Score[i] = constant.MaximumValue
		}
	}
}

// getF0CandidateContour() calculates the F0 candidates based on the
// zero-crossings.
func (s *Session) getF0CandidateContour(boundaryF0 float64) {
	if 0 == checkEvent(s.zeroCrossings.numberOfNegatives-2)*
		checkEvent(s.zeroCrossings.numberOfPositives-2)*
		checkEvent(s.zeroCrossings.numberOfPeaks-2)*
		checkEvent(s.zeroCrossings.numberOfDips-2) {
		for i := 0; i < s.f0Length; i++ {
			s.f0Score[i] = constant.MaximumValue
			s.f0Candidate[i] = 0.0
		}
		return
	}

	var interpolatedF0Set [4][]float64
	for i := 0; i < 4; i++ {
		interpolatedF0Set[i] = make([]float64, s.f0Length)
	}

	matlab.Interp1(s.zeroCrossings.negativeIntervalLocations[:s.zeroCrossings.numberOfNegatives],
		s.zeroCrossings.negativeIntervals,
		s.temporalPositions, interpolatedF0Set[0])
	matlab.Interp1(s.zeroCrossings.positiveIntervalLocations[:s.zeroCrossings.numberOfPositives],
		s.zeroCrossings.positiveIntervals,
		s.temporalPositions, interpolatedF0Set[1])
	matlab.Interp1(s.zeroCrossings.peakIntervalLocations[:s.zeroCrossings.numberOfPeaks],
		s.zeroCrossings.peakIntervals,
		s.temporalPositions, interpolatedF0Set[2])
	matlab.Interp1(s.zeroCrossings.dipIntervalLocations[:s.zeroCrossings.numberOfDips],
		s.zeroCrossings.dipIntervals,
		s.temporalPositions, interpolatedF0Set[3])

	s.getF0CandidateContourSub(interpolatedF0Set, boundaryF0)
}

// getF0CandidateFromRawEvent() calculates F0 candidate contour in 1-ch signal
func (s *Session) getF0CandidateFromRawEvent(boundaryF0 float64) {
	filteredSignal := make([]float64, s.fftSize)
	s.getFilteredSignal(matlab.Round(s.actualFS/boundaryF0/2.0), filteredSignal)
	s.getFourZeroCrossingIntervals(filteredSignal)
	s.getF0CandidateContour(boundaryF0)
}

// getF0CandidatesAndScores calculates all f0 candidates and their scores.
func (s *Session) getF0CandidatesAndScores() {
	// Calculation of the acoustics events (zero-crossing)
	for i := 0; i < s.numberOfBands; i++ {
		s.getF0CandidateFromRawEvent(s.boundaryF0List[i])
		for j := 0; j < s.f0Length; j++ {
			// A way to avoid zero division
			s.f0Scores[i][j] = s.f0Score[j] / (s.f0Candidate[j] + constant.MySafeGuardMinimum)
			s.f0Candidates[i][j] = s.f0Candidate[j]
		}
	}
}
