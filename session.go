package dio

import "math"

type Session struct {
	*Estimator
	params *params
	time   int
}

// NewSession creates new Session.
func NewSession(fs float64, option *Option) *Session {
	if option == nil {
		option = NewOption()
	}
	t := math.Ceil(fs / option.F0Floor * 20) // 最長周期の20倍 [sample]
	fp := option.FramePeriod / 1000.0        // [sec]

	// f0Length = int(xLength/(fs*FramePeriod)) + 1
	// xLength をf0検出単位の整数倍（4の倍数倍）に取る
	b := math.Round(fs * fp * 4)       // [sample]
	xLength := int(math.Ceil(t/b) * b) // [sample]

	p := newParams(xLength, fs, option)
	s := &Session{
		params:    p,
		Estimator: newEstimator(p),
	}
	return s
}

func (s *Session) Len() int {
	return s.Estimator.params.xLength
}

func (s *Session) F0Length() int {
	return s.Estimator.params.f0Length
}

func (s *Session) FramePeriod() float64 {
	return s.Estimator.params.option.FramePeriod
}

func (s *Session) Next(x []float64) []float64 {
	if len(x) < s.Len() {
		x = append(x, make([]float64, s.Len()-len(x))...)
	}
	if s.Len() < len(x) {
		x = x[:s.Len()]
	}
	s.Estimator = newEstimator(s.params)
	s.Estimator.x = x
	_, f0 := s.Estimator.Estimate()
	return f0
}
