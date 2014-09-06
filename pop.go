package bio

// file pop.go - functions that model population growth

// LucasU computes term n of Lucas sequence U(P,Q)
//
// Recurrence relation, using U for U(P,Q):
//
//  U₀ := 0
//  U₁ := 1
//  Uₙ := P*Uₙ₋₁ - Q*Uₙ₋₂
func LucasU(n, p, q int) int {
	if n < 1 {
		return 0
	}
	m1, m2 := 1, 0
	for n--; n > 0; n-- {
		m1, m2 = p*m1-q*m2, m1
	}
	return m1
}
