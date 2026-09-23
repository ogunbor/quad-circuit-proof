use ark_bn254::Fr;
use ark_ff::{Field, Zero};
use std::ops::{Add, Mul, Sub};

#[derive(Clone, Debug)]
pub struct Poly {
    pub coeffs: Vec<Fr>,
}

impl Poly {
    pub fn new(coeffs: Vec<Fr>) -> Self {
        let mut p = Poly { coeffs };
        p.trim();
        p
    }

    pub fn zero() -> Self {
        Poly::new(vec![])
    }

    fn trim(&mut self) {
        while let Some(last) = self.coeffs.last() {
            if last.is_zero() && self.coeffs.len() > 1 {
                self.coeffs.pop();
            } else {
                break;
            }
        }
        if self.coeffs.is_empty() {
            self.coeffs.push(Fr::zero());
        }
    }

    pub fn degree(&self) -> usize {
        self.coeffs.len().saturating_sub(1)
    }

    pub fn eval(&self, x: Fr) -> Fr {
        let mut acc = Fr::zero();
        for c in self.coeffs.iter().rev() {
            acc = acc * x + c;
        }
        acc
    }

    pub fn scale(&self, s: Fr) -> Poly {
        Poly::new(self.coeffs.iter().map(|c| *c * s).collect())
    }

    pub fn interpolate(points: &[(Fr, Fr)]) -> Poly {
        let mut result = Poly::zero();
        for (i, &(xi, yi)) in points.iter().enumerate() {
            let mut basis = Poly::new(vec![Fr::from(1u64)]);
            let mut denom = Fr::from(1u64);
            for (j, &(xj, _)) in points.iter().enumerate() {
                if i == j {
                    continue;
                }
                let factor = Poly::new(vec![-xj, Fr::from(1u64)]);
                basis = &basis * &factor;
                denom *= xi - xj;
            }
            let scale = yi * denom.inverse().expect("distinct domain points");
            let scaled = basis.scale(scale);
            result = &result + &scaled;
        }
        result
    }

    pub fn div_rem(&self, divisor: &Poly) -> (Poly, Poly) {
        let mut remainder = self.clone();
        let div_deg = divisor.degree();
        let div_lead = *divisor.coeffs.last().unwrap();
        let mut quotient_coeffs = vec![Fr::zero(); 1];

        while remainder.coeffs.iter().any(|c| !c.is_zero()) && remainder.degree() >= div_deg {
            let lead = *remainder.coeffs.last().unwrap();
            let coeff = lead * div_lead.inverse().expect("nonzero leading coeff");
            let shift = remainder.degree() - div_deg;

            if quotient_coeffs.len() < shift + 1 {
                quotient_coeffs.resize(shift + 1, Fr::zero());
            }
            quotient_coeffs[shift] = coeff;

            let mut sub_coeffs = vec![Fr::zero(); shift];
            for c in &divisor.coeffs {
                sub_coeffs.push(*c * coeff);
            }
            remainder = &remainder - &Poly::new(sub_coeffs);
        }

        (Poly::new(quotient_coeffs), remainder)
    }
}

impl Add for &Poly {
    type Output = Poly;
    fn add(self, rhs: &Poly) -> Poly {
        let n = self.coeffs.len().max(rhs.coeffs.len());
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let a = *self.coeffs.get(i).unwrap_or(&Fr::zero());
            let b = *rhs.coeffs.get(i).unwrap_or(&Fr::zero());
            out.push(a + b);
        }
        Poly::new(out)
    }
}

impl Sub for &Poly {
    type Output = Poly;
    fn sub(self, rhs: &Poly) -> Poly {
        let n = self.coeffs.len().max(rhs.coeffs.len());
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let a = *self.coeffs.get(i).unwrap_or(&Fr::zero());
            let b = *rhs.coeffs.get(i).unwrap_or(&Fr::zero());
            out.push(a - b);
        }
        Poly::new(out)
    }
}

impl Mul for &Poly {
    type Output = Poly;
    fn mul(self, rhs: &Poly) -> Poly {
        let mut out = vec![Fr::zero(); self.coeffs.len() + rhs.coeffs.len() - 1];
        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in rhs.coeffs.iter().enumerate() {
                out[i + j] += *a * *b;
            }
        }
        Poly::new(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolation_passes_through_given_points() {
        let p = Poly::interpolate(&[
            (Fr::from(5u64), Fr::from(0u64)),
            (Fr::from(7u64), Fr::from(1u64)),
        ]);
        assert_eq!(p.eval(Fr::from(5u64)), Fr::from(0u64));
        assert_eq!(p.eval(Fr::from(7u64)), Fr::from(1u64));
        assert_eq!(p.degree(), 1);
    }

    #[test]
    fn division_recovers_exact_quotient_with_zero_remainder() {
        let a = Poly::new(vec![-Fr::from(5u64), Fr::from(1u64)]);
        let b = Poly::new(vec![-Fr::from(7u64), Fr::from(1u64)]);
        let product = &a * &b;

        let (q, r) = product.div_rem(&a);
        assert_eq!(q.coeffs, b.coeffs);
        assert!(r.coeffs.iter().all(|c| c.is_zero()));
    }
}
