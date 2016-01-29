// A seekable reader.
package mmapreadseeker

import (
	"io"

	"golang.org/x/exp/mmap"
)

type mreader struct {
	*mmap.ReaderAt
	offset int64
}

func (m *mreader) Read(p []byte) (int, error) {
	n, e := m.ReadAt(p, m.offset)
	if n > 0 && e == io.EOF {
		e = nil
	}
	m.offset += int64(n)
	return n, e
}

func (m *mreader) Seek(offset int64, whence int) (int64, error) {
	if whence == 0 {
		m.offset = offset
	} else if whence == 1 {
		m.offset += offset
	} else {
		m.offset = int64(m.Len()) + offset
	}
	return m.offset, nil
}
