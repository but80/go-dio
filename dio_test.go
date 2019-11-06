package dio

import (
	"fmt"
	"math"
	"os/exec"
	"strconv"
	"testing"

	"github.com/but80/go-dio/internal/testasset"
	"github.com/stretchr/testify/assert"
)

// from https://github.com/lazybeaver/xorshift

type xorShift64Star struct {
	state uint64
}

func (r *xorShift64Star) next() uint64 {
	r.state ^= (r.state >> 12)
	r.state ^= (r.state << 25)
	r.state ^= (r.state >> 27)
	return r.state * 2685821657736338717
}

func (r *xorShift64Star) nextFloat64() float64 {
	return float64(r.next()) / float64(math.MaxUint64)
}

const (
	noiseDb = -18
)

func TestDio(t *testing.T) {
	cmd := exec.Command(
		"./tools/make-testdata",
		"Dio",
		"./tools/testdata.wav",
		strconv.Itoa(noiseDb),
	)
	want, err := cmd.Output()
	assert.NoError(t, err)

	x, fs, err := testasset.Load("./tools/testdata.wav")
	assert.NoError(t, err)

	noiseAmp := math.Pow(2, noiseDb)
	r := xorShift64Star{state: 1}
	for i := range x {
		v := r.nextFloat64()
		x[i] += (v - .5) * 2 * noiseAmp
	}

	option := NewOption()
	s := NewSession(x, float64(fs), option)
	temporalPositions, f0 := s.Estimate()

	got := ""
	for i := range f0 {
		got += fmt.Sprintf("%05.3f: %05.1f\n", temporalPositions[i], f0[i])
	}

	assert.Equal(t, string(want), got)
}
