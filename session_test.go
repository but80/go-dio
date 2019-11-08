package dio

import (
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"math/rand"
	"os/exec"
	"path/filepath"
	"regexp"
	"strconv"
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
		log.Printf("path = %s", path)

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
		want = rx1.ReplaceAllString(want, "") // remove temporal positions

		wantsa := strings.Split(want, "\n")
		var wantfa []float64
		for _, s := range wantsa {
			v, err := strconv.ParseFloat(s, 64)
			if err != nil {
				continue
			}
			wantfa = append(wantfa, v)
		}

		x, fs, err := testasset.Load(path)
		assert.NoError(t, err)

		s := NewSession(float64(fs), nil)
		in, out := s.Start()
		go func() {
			r := rand.New(rand.NewSource(1))
			k := len(x) / 100
			for 0 < len(x) {
				n := r.Intn(k) + 1
				if len(x) < n {
					n = len(x)
				}
				in <- x[:n]
				x = x[n:]
				// time.Sleep(time.Millisecond * 10)
			}
			close(in)
		}()
		var result []float64
		for f0 := range out {
			result = append(result, f0...)
		}

		diffSum := .0
		diffCount := 0
		zeroCount := 0
		got := ""
		for i, g := range result {
			got += fmt.Sprintf("%05.1f\n", g)
			if len(wantfa) <= i {
				continue
			}
			w := wantfa[i]
			if g == 0 && w == 0 {
				diffCount++
				continue
			}
			if g == 0 && w != 0 || g != 0 && w == 0 {
				zeroCount++
				continue
			}
			d := math.Log2(g/w) * 12
			diffSum += d * d
			diffCount++
		}
		diff := math.Sqrt(diffSum / float64(diffCount))
		zero := float64(zeroCount) / float64(zeroCount+diffCount)
		assert.True(t, diff < .05, fmt.Sprintf("RMSE of pitch = %f [semitone]", diff))
		assert.True(t, zero*100 < .5, fmt.Sprintf("Ratio of frames errored to be silent or not = %f [%%]", zero*100))

		// Output CSV
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
