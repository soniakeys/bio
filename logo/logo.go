package logo

import (
	"bytes"
	"fmt"
	"strconv"

	"github.com/ajstarks/svgo"
)

type Config struct {
	TopMargin   int    // in canvas units
	BotMargin   int    // labels are centered vertically within this space
	SideMargin  int    // in canvas units
	ColWidth    int    // width in canvas units of one column of the profile
	ColHeight   int    // height in canvas units allocated for all four bases
	LabelHt     int    // centered within botMargin
	TicHt       int    // included in botMargin
	LineStyle   string // line style for x-axis and tics
	BaseStyle   string // font style for base symbols
	BaseColors  [4]string
	BaseSymbols [4]string
}

var Defaults = Config{
	TopMargin:   20,
	BotMargin:   20,
	SideMargin:  20,
	ColWidth:    36,
	ColHeight:   60,
	LabelHt:     14,
	TicHt:       3,
	LineStyle:   "stroke:black",
	BaseColors:  [4]string{"purple", "blue", "green", "red"},
	BaseSymbols: [4]string{"A", "C", "T", "G"},
	BaseStyle:   "font-family:sans-serif",
}

func Margins(m int) func(*Config) {
	return func(cf *Config) {
		cf.TopMargin = m
		cf.BotMargin = m
		cf.SideMargin = m
	}
}

func TopMargin(m int) func(*Config) {
	return func(cf *Config) {
		cf.TopMargin = m
	}
}

func BotMargin(m int) func(*Config) {
	return func(cf *Config) {
		cf.BotMargin = m
	}
}

func SideMargin(m int) func(*Config) {
	return func(cf *Config) {
		cf.SideMargin = m
	}
}

func ColWidth(w int) func(*Config) {
	return func(cf *Config) {
		cf.ColWidth = w
	}
}

func ColHeight(w int) func(*Config) {
	return func(cf *Config) {
		cf.ColHeight = w
	}
}

func LabelHt(h int) func(*Config) {
	return func(cf *Config) {
		cf.LabelHt = h
	}
}

func TicHt(h int) func(*Config) {
	return func(cf *Config) {
		cf.TicHt = h
	}
}

func LineStyle(s string) func(*Config) {
	return func(cf *Config) {
		cf.LineStyle = s
	}
}

func BaseColors(c [4]string) func(*Config) {
	return func(cf *Config) {
		cf.BaseColors = c
	}
}

func BaseStyle(s string) func(*Config) {
	return func(cf *Config) {
		cf.BaseStyle = s
	}
}

func RNA() func(*Config) {
	return func(cf *Config) {
		cf.BaseSymbols = [4]string{"A", "C", "U", "G"}
	}
}

func Set(c *Config) func(*Config) {
	return func(cf *Config) {
		*cf = *c
	}
}

// Motif generates a motif logo as SVG.
//
// Argument p is a bio.FracProfile.
// Motif has usable defaults and does not require options.  To modify
// rendering though, use any of the package functions that return
// func(*Config).  Supply these functions as variadic options, giving
// each function its arguments as appropriate to override Config defaults.
// Options are evaluated in order so that later options override earlier ones.
// See in the package example how Margins(5) overrides the defaults, then
// BotMargin(20) overrides the Config.BotMargin set by Margins.
func Motif(p [][4]float64, options ...func(*Config)) []byte {
	cf := Defaults // copy defaults
	for _, o := range options {
		o(&cf)
	}
	var b bytes.Buffer
	s := svg.New(&b)
	width := cf.SideMargin*2 + cf.ColWidth*len(p)
	height := cf.TopMargin + cf.ColHeight + cf.BotMargin
	s.Start(width, height)
	// draw baseline x-axis
	axis := cf.TopMargin + cf.ColHeight
	s.Line(cf.SideMargin, axis, width-cf.SideMargin, axis, cf.LineStyle)
	// draw tics and position labels
	x := cf.SideMargin
	for i := 0; i <= len(p); i++ {
		// tics
		s.Line(x, axis, x, axis+cf.TicHt, cf.LineStyle)
		// labels
		s.Text(x, height-(cf.BotMargin-cf.LabelHt)/2, strconv.Itoa(i),
			fmt.Sprintf("text-anchor:middle;font-size:%d", cf.LabelHt))
		x += cf.ColWidth
	}
	// draw logo
	x = cf.SideMargin
	sx := float64(cf.ColWidth) / 12
	for _, dist := range p {
		y := float64(cf.TopMargin)
		for bx, fr := range dist {
			bHt := float64(cf.ColHeight) * fr
			y += bHt
			s.Gtransform(fmt.Sprintf("translate(%d,%d) scale(%.2f,%.2f)",
				x, int(y), sx, bHt/12))
			s.Text(0, 0, cf.BaseSymbols[bx],
				cf.BaseStyle+";fill:"+cf.BaseColors[bx])
			s.Gend()
		}
		x += cf.ColWidth
	}
	s.End()
	return b.Bytes()
}
