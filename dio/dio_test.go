package dio

import (
	"fmt"
	"os/exec"
	"testing"

	"github.com/but80/go-world/internal/testasset"
	"github.com/stretchr/testify/assert"
)

func TestDio(t *testing.T) {
	cmd := exec.Command("../tools/make-testdata", "Dio", "../tools/testdata.wav")
	want, err := cmd.Output()
	assert.NoError(t, err)

	x, fs, err := testasset.Load("../tools/testdata.wav")
	assert.NoError(t, err)
	option := NewOption()

	s := NewSession(x, fs, option)
	temporalPositions, f0 := s.Estimate()

	got := ""
	for i := range f0 {
		got += fmt.Sprintf("%05.3f: %07.3f\n", temporalPositions[i], f0[i])
	}

	assert.Equal(t, string(want), got)
}
