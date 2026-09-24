use crate::polynomial::Poly;
use crate::r1cs::R1cs;
use ark_bn254::Fr;
use ark_ff::Zero;

pub struct Qap {
    pub domain: Vec<Fr>,
    pub t: Poly,
}

impl Qap {
    pub fn new(r1cs: &R1cs, domain: Vec<Fr>) -> Self {
        assert_eq!(domain.len(), r1cs.num_constraints());

        let mut t = Poly::new(vec![Fr::from(1u64)]);
        for &d in &domain {
            let factor = Poly::new(vec![-d, Fr::from(1u64)]);
            t = &t * &factor;
        }

        Qap { domain, t }
    }

    pub fn column_polys(&self, matrix: &[Vec<i64>], num_vars: usize) -> Vec<Poly> {
        (0..num_vars)
            .map(|j| {
                let points: Vec<(Fr, Fr)> = self
                    .domain
                    .iter()
                    .enumerate()
                    .map(|(i, &d)| {
                        let coeff = matrix[i][j];
                        let v = if coeff >= 0 {
                            Fr::from(coeff as u64)
                        } else {
                            -Fr::from((-coeff) as u64)
                        };
                        (d, v)
                    })
                    .collect();
                Poly::interpolate(&points)
            })
            .collect()
    }

    fn combine_with_witness(basis: &[Poly], w: &[Fr]) -> Poly {
        let mut acc = Poly::zero();
        for (poly, wj) in basis.iter().zip(w.iter()) {
            acc = &acc + &poly.scale(*wj);
        }
        acc
    }

    pub fn build_abc(&self, r1cs: &R1cs, w: &[Fr]) -> (Poly, Poly, Poly) {
        let a_basis = self.column_polys(&r1cs.a, r1cs.num_vars);
        let b_basis = self.column_polys(&r1cs.b, r1cs.num_vars);
        let c_basis = self.column_polys(&r1cs.c, r1cs.num_vars);

        (
            Self::combine_with_witness(&a_basis, w),
            Self::combine_with_witness(&b_basis, w),
            Self::combine_with_witness(&c_basis, w),
        )
    }

    pub fn compute_h(&self, a: &Poly, b: &Poly, c: &Poly) -> Poly {
        let ab = a * b;
        let numerator = &ab - c;
        let (h, rem) = numerator.div_rem(&self.t);
        assert!(
            rem.coeffs.iter().all(|x| x.is_zero()),
            "non-zero remainder: witness does not satisfy the QAP identity"
        );
        h
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derives_abc_and_h_for_x_equals_4() {
        let r1cs = R1cs::circuit();
        let w = r1cs.build_witness(4);
        let qap = Qap::new(&r1cs, vec![Fr::from(5u64), Fr::from(7u64)]);

        let (a, b, c) = qap.build_abc(&r1cs, &w);
        let _h = qap.compute_h(&a, &b, &c);
    }

    #[test]
    #[should_panic(expected = "non-zero remainder")]
    fn tampered_witness_breaks_the_qap_identity() {
        let r1cs = R1cs::circuit();
        let mut w = r1cs.build_witness(4);
        w[1] += Fr::from(1u64);
        let qap = Qap::new(&r1cs, vec![Fr::from(5u64), Fr::from(7u64)]);
        let (a, b, c) = qap.build_abc(&r1cs, &w);
        qap.compute_h(&a, &b, &c);
    }
}
