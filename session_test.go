package dio

import (
	"fmt"
	"io/ioutil"
	"log"
	"os/exec"
	"path/filepath"
	"regexp"
	"strings"
	"testing"

	"github.com/but80/go-dio/internal/testasset"
	"github.com/stretchr/testify/assert"
)

var rx1 = regexp.MustCompile(`(?m)^[^:]*: *`)

func TestSession(t *testing.T) {
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

		cmd := exec.Command(
			"./tools/make-testdata",
			"Dio",
			path,
			"-999",
			"0",
		)
		wantb, err := cmd.Output()
		assert.NoError(t, err)
		want := string(wantb)
		want = rx1.ReplaceAllString(want, "")

		x, fs, err := testasset.Load(path)
		assert.NoError(t, err)

		s := NewSession(float64(fs), nil)
		log.Printf("%#v", s.Estimator.params)
		step := s.Len() / 2
		var result []float64
		for i := 0; i+step < len(x); i += step {
			j := i + s.Len()
			if len(x) < j {
				j = len(x)
			}
			f0 := s.Next(x[i:j])
			from := 0
			if 0 < i {
				from = s.F0Length() / 4
			}
			to := s.F0Length()
			if i+step+step < len(x) {
				to = s.F0Length() * 3 / 4
			}
			result = append(result, f0[from:to]...)
		}

		got := ""
		for i := range result {
			g := fmt.Sprintf("%05.1f\n", result[i])
			// log.Print(g)
			got += g
		}
		assert.Equal(t, want, got, file.Name())

		wanta := strings.Split(want, "\n")
		gota := strings.Split(got, "\n")
		var csvn int
		if csvn < len(gota) {
			csvn = len(gota)
			wanta = append(wanta, make([]string, csvn-len(wanta))...)
		} else {
			csvn = len(wanta)
			gota = append(gota, make([]string, csvn-len(gota))...)
		}
		csv := "want,got\r\n"
		for i := range wanta {
			csv += fmt.Sprintf("%s,%s\r\n", wanta[i], gota[i])
		}
		assert.NoError(t, ioutil.WriteFile("testdata/"+file.Name()+".csv", []byte(csv), 0644))
	}
}
