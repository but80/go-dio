package dio

import (
	"math"

	"github.com/but80/go-dio/constant"
	"github.com/but80/go-dio/internal/common"
	"github.com/but80/go-dio/internal/matlab"
	"gonum.org/v1/gonum/fourier"
)

// Option is the struct to order the parameter for Dio.
type Option struct {
	F0Floor          float64
	F0Ceil           float64
	ChannelsInOctave float64
	FramePeriod      float64 // msec
	Speed            int     // (1, 2, ..., 12)
	AllowedRange     float64 // Threshold used for fixing the F0 contour.
}

// NewOption creates new Option with the default parameters.
func NewOption() *Option {
	option := &Option{}

	// You can change default parameters.
	option.ChannelsInOctave = 2.0
	option.F0Ceil = constant.CeilF0
	option.F0Floor = constant.FloorF0
	option.FramePeriod = 5

	// You can use the value from 1 to 12.
	// Default value 11 is for the fs of 44.1 kHz.
	// The lower value you use, the better performance you can obtain.
	option.Speed = 1

	// You can give a positive real number as the threshold.
	// The most strict value is 0, and there is no upper limit.
	// On the other hand, I think that the value from 0.02 to 0.2 is reasonable.
	option.AllowedRange = 0.1

	return option
}

// Session is the struct holds the variables needed to estimate f0 by Dio.
type Session struct {
	// Inputs
	x      []float64 // Input signal
	fs     float64   // Sampling frequency
	option *Option   // Struct to order the parameter for DIO

	// Immutable temporaries
	f0Length          int
	yLength           int
	fftSize           int
	fft               *fourier.FFT
	numberOfBands     int // Max number of candidates
	voiceRangeMinimum int // Number of consecutive frames for stable estimation

	// Mutable temporaries
	ySpectrum      []complex128
	boundaryF0List []float64
	zeroCrossings  *zeroCrossings
	f0Candidates   [][]float64
	f0Scores       [][]float64
	f0Candidate    []float64
	f0Score        []float64

	// Outputs
	temporalPositions []float64 // Temporal positions.
	f0                []float64 // F0 contour.
}

// NewSession creates new Session.
func NewSession(x []float64, fs float64, option *Option) *Session {
	numberOfBands := 1 + int(math.Log2(option.F0Ceil/option.F0Floor)*option.ChannelsInOctave)
	f0Length := int(1000.0*float64(len(x))/fs/option.FramePeriod) + 1
	temporalPositions := make([]float64, f0Length)
	f0 := make([]float64, f0Length)

	boundaryF0List := make([]float64, numberOfBands)
	for i := 0; i < numberOfBands; i++ {
		boundaryF0List[i] = option.F0Floor * math.Pow(2.0, float64(i+1)/option.ChannelsInOctave)
	}

	s := &Session{
		x:      x,
		fs:     fs,
		option: option,

		f0Length:          f0Length,
		numberOfBands:     numberOfBands,
		voiceRangeMinimum: int(0.5+1000.0/option.FramePeriod/option.F0Floor)*2 + 1,

		boundaryF0List: boundaryF0List,
		f0Candidate:    make([]float64, f0Length),
		f0Score:        make([]float64, f0Length),

		temporalPositions: temporalPositions,
		f0:                f0,
	}

	// normalization
	s.yLength = 1 + len(s.x)
	s.fftSize = common.GetSuitableFFTSize(s.yLength +
		matlab.Round(s.fs/constant.CutOff)*2 + 1 +
		4*int(1.0+s.fs/s.boundaryF0List[0]/2.0))
	s.fft = fourier.NewFFT(s.fftSize)

	s.ySpectrum = make([]complex128, s.fftSize)
	s.zeroCrossings = newZeroCrossings(s.yLength)
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

// Estimate estimates the F0 based on Distributed Inline-filter Operation.
func (s *Session) Estimate() ([]float64, []float64) {
	// Calculation of the spectrum used for the f0 estimation
	s.getSpectrumForEstimation()

	s.getF0CandidatesAndScores()

	// Selection of the best value based on fundamental-ness.
	// This function is related with SortCandidates() in MATLAB.
	bestF0Contour := make([]float64, s.f0Length)
	s.getBestF0Contour(bestF0Contour)

	// Postprocessing to find the best f0-contour.
	s.fixF0Contour(bestF0Contour, s.f0)

	return s.temporalPositions, s.f0
}
