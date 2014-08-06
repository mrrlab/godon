package main

import (
	"fmt"

	"bitbucket.com/Davydov/golh/tree"
)

func M0(cali CodonSequences, t *tree.Tree, cf CodonFrequency, kappa, omega float64) float64 {
	Q := createTransitionMatrix(cf, kappa, omega)
	Qs := make([][]*EMatrix, 1)
	Qs[0] = make([]*EMatrix, t.NNodes())
	em := NewEMatrix(Q)
	err := em.Eigen()
	if err != nil {
		panic(fmt.Sprintf("error finding eigen: %v", err))
	}
	for i := 0; i < len(Qs[0]); i++ {
		Qs[0][i] = em
	}
	return L(cali, t, []float64{1}, Qs, cf)
}

func M3(cali CodonSequences, t *tree.Tree, cf CodonFrequency, kappa float64, omega0, omega1, omega2 float64, p0, p1, p2 float64) float64 {
	fmt.Printf("omega0=[%f, %f, %f], p=[%f, %f, %f]\n", omega0, omega1, omega2, p0, p1, p2)
	Q0 := createTransitionMatrix(cf, kappa, omega0)
	Q1 := createTransitionMatrix(cf, kappa, omega1)
	Q2 := createTransitionMatrix(cf, kappa, omega2)
	em0 := NewEMatrix(Q0)
	em1 := NewEMatrix(Q1)
	em2 := NewEMatrix(Q2)
	err1 := em0.Eigen()
	err2 := em1.Eigen()
	err3 := em2.Eigen()
	if err1 != nil || err2 != nil || err3 != nil {
		panic(fmt.Sprintf("error finding eigen: %v, %v, %v", err1, err2, err3))
	}
	Qs := make([][]*EMatrix, 3)
	for i := 0; i < len(Qs); i++ {
		Qs[i] = make([]*EMatrix, t.NNodes())
		for b, _ := range Qs[i] {
			switch i {
			case 0:
				Qs[i][b] = em0
			case 1:
				Qs[i][b] = em0
			case 2:
				Qs[i][b] = em1
			}
		}
	}
	return L(cali, t, []float64{p0, p1, p2}, Qs, cf)

}

func H1(cali CodonSequences, t *tree.Tree, cf CodonFrequency, fg int, kappa float64, omega0, omega2 float64, p0, p1, p2a, p2b float64) float64 {
	fmt.Printf("fg=%d, kappa=%f, omega0=%f, omega2=%f, p=[%f, %f, %f, %f]\n", fg, kappa, omega0, omega2, p0, p1, p2a, p2b)
	Q0 := createTransitionMatrix(cf, kappa, omega0)
	Q1 := createTransitionMatrix(cf, kappa, 1)
	Q2 := createTransitionMatrix(cf, kappa, omega2)
	em0 := NewEMatrix(Q0)
	em1 := NewEMatrix(Q1)
	em2 := NewEMatrix(Q2)
	err1 := em0.Eigen()
	err2 := em1.Eigen()
	err3 := em2.Eigen()
	if err1 != nil || err2 != nil || err3 != nil {
		panic(fmt.Sprintf("error finding eigen: %v, %v, %v", err1, err2, err3))
	}
	Qs := make([][]*EMatrix, 4)
	for i := 0; i < len(Qs); i++ {
		Qs[i] = make([]*EMatrix, t.NNodes())
		for b, _ := range Qs[i] {
			switch i {
			case 0:
				Qs[i][b] = em0
			case 1:
				Qs[i][b] = em1
			case 2:
				if b != fg {
					Qs[i][b] = em0
				} else {
					Qs[i][b] = em2
				}
			case 3:
				if b != fg {
					Qs[i][b] = em1
				} else {
					Qs[i][b] = em2
				}
			}
		}
	}
	return L(cali, t, []float64{p0, p1, p2a, p2b}, Qs, cf)
}
