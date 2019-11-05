package dio

import (
	"math"

	"github.com/but80/go-world/constant"
)

// getBestF0Contour calculates the best f0 contour based on scores of
// all candidates. The F0 with highest score is selected.
func (s *Session) getBestF0Contour(bestF0Contour []float64) {
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
func (s *Session) fixStep1(f0In []float64, f0Out []float64) {
	minRange := s.voiceRangeMinimum
	f0Base := make([]float64, s.f0Length)
	copy(
		f0Base[minRange:s.f0Length-minRange],
		f0In[minRange:s.f0Length-minRange],
	)

	// Processing to prevent the jumping of f0
	for i := 0; i < minRange; i++ {
		f0Out[i] = 0.0
	}
	for i := minRange; i < s.f0Length; i++ {
		if math.Abs((f0Base[i]-f0Base[i-1])/
			(constant.MySafeGuardMinimum+f0Base[i])) <
			s.option.AllowedRange {
			f0Out[i] = f0Base[i]
		} else {
			f0Out[i] = 0.0
		}
	}
}

// fixStep2 is the 2nd step of the postprocessing.
// This function eliminates the suspected f0 in the anlaut and auslaut.
func (s *Session) fixStep2(f0In []float64, f0Out []float64) {
	copy(f0Out, f0In)

	center := (s.voiceRangeMinimum - 1) / 2
	for i := center; i < s.f0Length-center; i++ {
		for j := -center; j <= center; j++ {
			if f0In[i+j] == 0 {
				f0Out[i] = 0.0
				break
			}
		}
	}
}

// getNumberOfVoicedSections counts the number of voiced sections.
func (s *Session) getNumberOfVoicedSections(f0 []float64, positiveIndex, negativeIndex []int) (int, int) {
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
func (s *Session) selectBestF0(currentF0, pastF0 float64, targetIndex int) float64 {
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

// fixStep3 is the 3rd step of the postprocessing.
// This function corrects the f0 candidates from backward to forward.
func (s *Session) fixStep3(f0In []float64, negativeIndex []int, negativeCount int, f0Out []float64) {
	copy(f0Out, f0In)

	for i := 0; i < negativeCount; i++ {
		var limit int
		if i == negativeCount-1 {
			limit = s.f0Length - 1
		} else {
			limit = negativeIndex[i+1]
		}
		for j := negativeIndex[i]; j < limit; j++ {
			f0Out[j+1] = s.selectBestF0(f0Out[j], f0Out[j-1], j+1)
			if f0Out[j+1] == 0 {
				break
			}
		}
	}
}

// fixStep4 is the 4th step of the postprocessing.
// This function corrects the f0 candidates from forward to backward.
func (s *Session) fixStep4(f0In []float64, positiveIndex []int, positiveCount int, f0Out []float64) {
	copy(f0Out, f0In)

	for i := positiveCount - 1; i >= 0; i-- {
		limit := 1
		if i != 0 {
			limit = positiveIndex[i-1]
		}
		for j := positiveIndex[i]; j > limit; j-- {
			f0Out[j-1] = s.selectBestF0(f0Out[j], f0Out[j+1], j-1)
			if f0Out[j-1] == 0 {
				break
			}
		}
	}
}

// fixF0Contour() calculates the definitive f0 contour based on all f0
// candidates. There are four steps.
func (s *Session) fixF0Contour(bestF0Contour []float64) {
	if s.f0Length <= s.voiceRangeMinimum {
		return
	}

	f0Tmp1 := make([]float64, s.f0Length)
	f0Tmp2 := make([]float64, s.f0Length)

	s.fixStep1(bestF0Contour, f0Tmp1)
	s.fixStep2(f0Tmp1, f0Tmp2)

	positiveIndex := make([]int, s.f0Length)
	negativeIndex := make([]int, s.f0Length)
	positiveCount, negativeCount := s.getNumberOfVoicedSections(f0Tmp2, positiveIndex, negativeIndex)
	s.fixStep3(f0Tmp2, negativeIndex, negativeCount, f0Tmp1)
	s.fixStep4(f0Tmp1, positiveIndex, positiveCount, s.f0)
}
