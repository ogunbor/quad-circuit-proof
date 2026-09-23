//! R1CS (Rank-1 Constraint System) for the circuit y = x^2 + x + 9.
//!
//! Witness vector order: [one, y, x, s1], where s1 = x*x is the
//! intermediate value for the multiplication gate.
//!
//! Constraint 1: x * x = s1
//! Constraint 2: (9*one + x + s1) * one = y   <=>   s1 + x + 9 = y
//!
//! A row is satisfied when (A_row . w) * (B_row . w) == (C_row . w).

use ark_bn254::Fr;
use ark_ff::Zero;

#[derive(Clone, Debug)]
pub struct R1cs {
    pub a: Vec<Vec<i64>>,
    pub b: Vec<Vec<i64>>,
    pub c: Vec<Vec<i64>>,
    pub num_vars: usize,
}

impl R1cs {
    /// The circuit from the paper working: y = x^2 + x + 9.
    pub fn circuit() -> Self {
        R1cs {
            a: vec![vec![0, 0, 1, 0], vec![9, 0, 1, 1]],
            b: vec![vec![0, 0, 1, 0], vec![1, 0, 0, 0]],
            c: vec![vec![0, 0, 0, 1], vec![0, 1, 0, 0]],
            num_vars: 4,
        }
    }

    pub fn num_constraints(&self) -> usize {
        self.a.len()
    }

    fn dot(row: &[i64], w: &[Fr]) -> Fr {
        let mut acc = Fr::zero();
        for (coeff, wi) in row.iter().zip(w.iter()) {
            let c = if *coeff >= 0 {
                Fr::from(*coeff as u64)
            } else {
                -Fr::from((-coeff) as u64)
            };
            acc += c * wi;
        }
        acc
    }

    /// Checks (A_i . w) * (B_i . w) == (C_i . w) for every constraint row.
    /// Returns Err with the failing row index and both sides if any
    /// constraint doesn't hold.
    pub fn is_satisfied(&self, w: &[Fr]) -> Result<(), (usize, Fr, Fr)> {
        assert_eq!(w.len(), self.num_vars, "witness length mismatch");
        for i in 0..self.num_constraints() {
            let av = Self::dot(&self.a[i], w);
            let bv = Self::dot(&self.b[i], w);
            let cv = Self::dot(&self.c[i], w);
            let lhs = av * bv;
            if lhs != cv {
                return Err((i, lhs, cv));
            }
        }
        Ok(())
    }

    /// Builds [one, y, x, s1] for a given x, computing y = x^2 + x + 9
    /// and s1 = x^2 along the way.
    pub fn build_witness(&self, x: i64) -> Vec<Fr> {
        let one = Fr::from(1u64);
        let x = if x >= 0 {
            Fr::from(x as u64)
        } else {
            -Fr::from((-x) as u64)
        };
        let s1 = x * x;
        let y = s1 + x + Fr::from(9u64);
        vec![one, y, x, s1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn witness_for_x_4_matches_workings_and_satisfies_constraints() {
        let r1cs = R1cs::circuit();
        let w = r1cs.build_witness(4);

        assert_eq!(w[0], Fr::from(1u64)); // one
        assert_eq!(w[1], Fr::from(29u64)); // y
        assert_eq!(w[2], Fr::from(4u64)); // x
        assert_eq!(w[3], Fr::from(16u64)); // s1

        assert!(r1cs.is_satisfied(&w).is_ok());
    }

    #[test]
    fn tampered_witness_is_correctly_rejected() {
        let r1cs = R1cs::circuit();
        let mut w = r1cs.build_witness(4);
        w[1] += Fr::from(1u64); // wrong y so it no longer equals x^2+x+9
        assert!(r1cs.is_satisfied(&w).is_err());
    }
}

#[cfg(test)]
mod demo_wrong_witness {
    use super::*;

    #[test]
    fn mod13_witness_does_not_satisfy_real_field_r1cs() {
        let r1cs = R1cs::circuit();
        let w: Vec<Fr> = vec![
            Fr::from(1u64),
            Fr::from(3u64),
            Fr::from(4u64),
            Fr::from(3u64),
        ];
        let result = r1cs.is_satisfied(&w);
        println!("{:?}", result);
        assert!(result.is_err());
    }
}
