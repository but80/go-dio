package testasset

import (
	"io"

	"github.com/zenwerk/go-wave"
)

func Load(file string) ([]float64, int, error) {
	reader, err := wave.NewReader(file)
	if err != nil {
		return nil, 0, err
	}

	const eps = 1.0 / 32767.0
	data := make([]float64, 0, reader.NumSamples)
	for {
		v, err := reader.ReadSample()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, 0, err
		}
		a := v[0]
		if 1 < a {
			a -= 2.0
		}
		data = append(data, a)
	}
	return data, int(reader.FmtChunk.Data.SamplesPerSec), err
}
