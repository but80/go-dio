package dio

import (
	"fmt"
	"io/ioutil"
	"math"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
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
	testdataDir = "testdata"
	noiseDB     = -18
)

func TestDio(t *testing.T) {
	t.Parallel()

	files, err := ioutil.ReadDir(testdataDir)
	if err != nil {
		panic(err)
	}

	for _, file := range files {
		if !strings.HasSuffix(strings.ToLower(file.Name()), ".wav") {
			continue
		}
		path := filepath.Join(testdataDir, file.Name())

		x0, fs, err := testasset.Load(path)
		assert.NoError(t, err)

		s := New(x0, float64(fs), nil)

		for reverse := 0; reverse < 2; reverse++ {
			cmd := exec.Command(
				"./tools/make-testdata",
				"Dio",
				path,
				strconv.Itoa(noiseDB),
				strconv.Itoa(reverse),
			)
			want, err := cmd.Output()
			assert.NoError(t, err)

			x := make([]float64, len(x0))
			copy(x, x0)
			if reverse != 0 {
				for i := 0; i < len(x)/2; i++ {
					j := len(x) - 1 - i
					x[i], x[j] = x[j], x[i]
				}
			}

			noiseAmp := math.Pow(10, noiseDB/20.0)
			r := xorShift64Star{state: 1}
			for i := range x {
				v := r.nextFloat64()
				x[i] += (v - .5) * 2 * noiseAmp
			}

			s.x = x
			temporalPositions, f0 := s.Estimate()

			got := ""
			for i := range f0 {
				got += fmt.Sprintf("%05.3f: %05.1f\n", temporalPositions[i], f0[i])
			}

			assert.Equal(t, string(want), got, file.Name())
		}
	}
}
