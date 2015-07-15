package stat

import (
	//. "code.google.com/p/go-fn/fn"
	"math"
	//"math/big"
)

const DBL_EPSILON = 2.2204460492503131e-16

func Hypergeometric_CDF(k, K, N, n float64) float64 {
	if n > N {
		panic("n larger than population size")
	} else if k >= K || k >= n {
		return 1.0
	} else if k < 0.0 {
		return 0.0
	} else {
		midpoint := (n * K) / (N)

		if k >= midpoint {
			return 1 - Upper_Tail(k, K, N, n)
		} else {
			return Lower_Tail(k, K, N, n)
		}
	}
}

func Upper_Tail(k, K, N, n float64) float64 {
	s := Ran_Hypergeometric_PDF(k, K, N, n)
	Q := s

	for i := k + 1; i <= n; i++ {
		factor := ((K - i) / (i + 1.0)) * ((n - i) / (N - K + i + 1.0 - n))
		s *= factor
		Q += s
		relerr := s / Q
		if relerr < DBL_EPSILON {
			break
		}
	}

	return Q
}

func Lower_Tail(k, K, N, n float64) float64 {

	s := Ran_Hypergeometric_PDF(k, K, N, n)
	P := s

	for i := k; i > 0; i-- {
		factor := (i / (K - i + 1.0)) * ((N - K + i - n) / (n - i + 1.0))
		s *= factor
		P += s
		relerr := s / P
		if relerr < DBL_EPSILON {
			break
		}
	}

	return P
}

func Ran_Hypergeometric_PDF(k, K, N, n float64) float64 {
	if n > N {
		n = N
	}

	if k > K || k > n {
		return 0
	} else if n > N-K && k+N-K < n {
		return 0
	} else {

		//c1 := LnΓ(K+1) - LnΓ(k+1) - LnΓ(K-k+1)
		K_Lfact, _ := math.Lgamma(K + 1)
		k_Lfact, _ := math.Lgamma(k + 1)
		K_minus_k_Lfact, _ := math.Lgamma(K - k + 1)

		N_minus_K_Lfact, _ := math.Lgamma(N - K + 1)
		n_minus_k_Lfact, _ := math.Lgamma(n - k + 1)
		N_minus_K_minus_n_plus_k_Lfact, _ := math.Lgamma(N - K - n + k + 1)

		N_Lfact, _ := math.Lgamma(N + 1)
		n_Lfact, _ := math.Lgamma(n + 1)
		N_minus_n_Lfact, _ := math.Lgamma(N - n + 1)

		return math.Exp(K_Lfact - k_Lfact - K_minus_k_Lfact + N_minus_K_Lfact - n_minus_k_Lfact - N_minus_K_minus_n_plus_k_Lfact - N_Lfact + n_Lfact + N_minus_n_Lfact)
	}
}
