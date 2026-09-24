use crate::qap::Qap;
use crate::r1cs::R1cs;
use ark_bn254::{Fr, G1Projective, G2Projective};
use ark_ec::Group;
use ark_ff::{Field, UniformRand, Zero};
use ark_std::rand::RngCore;

pub struct ToxicWaste {
    pub alpha: Fr,
    pub beta: Fr,
    pub gamma: Fr,
    pub delta: Fr,
    pub tau: Fr,
}

impl ToxicWaste {
    pub fn sample<R: RngCore>(rng: &mut R) -> Self {
        ToxicWaste {
            alpha: Fr::rand(rng),
            beta: Fr::rand(rng),
            gamma: Fr::rand(rng),
            delta: Fr::rand(rng),
            tau: Fr::rand(rng),
        }
    }
}

pub struct Crs {
    pub alpha_g1: G1Projective,
    pub beta_g1: G1Projective,
    pub beta_g2: G2Projective,
    pub gamma_g2: G2Projective,
    pub delta_g1: G1Projective,
    pub delta_g2: G2Projective,
    pub tau_powers_g1: Vec<G1Projective>,
    pub tau_powers_g2: Vec<G2Projective>,
    pub public_ic_terms: Vec<G1Projective>,
    pub private_delta_terms: Vec<G1Projective>,
    pub h_over_delta_powers: Vec<G1Projective>,
    pub public_vars: Vec<usize>,
    pub private_vars: Vec<usize>,
}

impl Crs {
    pub fn setup(r1cs: &R1cs, qap: &Qap, tw: &ToxicWaste, public_vars: &[usize]) -> Crs {
        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();
        let deg_points = qap.domain.len();

        let mut tau_powers_g1 = Vec::new();
        let mut tau_powers_g2 = Vec::new();
        let mut acc = Fr::from(1u64);
        for _ in 0..deg_points {
            tau_powers_g1.push(g1 * acc);
            tau_powers_g2.push(g2 * acc);
            acc *= tw.tau;
        }

        let a_basis = qap.column_polys(&r1cs.a, r1cs.num_vars);
        let b_basis = qap.column_polys(&r1cs.b, r1cs.num_vars);
        let c_basis = qap.column_polys(&r1cs.c, r1cs.num_vars);

        let private_vars: Vec<usize> = (0..r1cs.num_vars)
            .filter(|j| !public_vars.contains(j))
            .collect();

        let combined_term = |j: usize| -> Fr {
            let u = a_basis[j].eval(tw.tau);
            let v = b_basis[j].eval(tw.tau);
            let w = c_basis[j].eval(tw.tau);
            tw.beta * u + tw.alpha * v + w
        };

        let gamma_inv = tw.gamma.inverse().expect("gamma is nonzero");
        let delta_inv = tw.delta.inverse().expect("delta is nonzero");

        let public_ic_terms: Vec<G1Projective> = public_vars
            .iter()
            .map(|&j| g1 * (combined_term(j) * gamma_inv))
            .collect();

        let private_delta_terms: Vec<G1Projective> = private_vars
            .iter()
            .map(|&j| g1 * (combined_term(j) * delta_inv))
            .collect();

        let t_tau = qap.t.eval(tw.tau);
        let mut h_over_delta_powers = Vec::new();
        let mut tau_pow = Fr::from(1u64);
        for _ in 0..deg_points {
            h_over_delta_powers.push(g1 * (tau_pow * t_tau * delta_inv));
            tau_pow *= tw.tau;
        }

        Crs {
            alpha_g1: g1 * tw.alpha,
            beta_g1: g1 * tw.beta,
            beta_g2: g2 * tw.beta,
            gamma_g2: g2 * tw.gamma,
            delta_g1: g1 * tw.delta,
            delta_g2: g2 * tw.delta,
            tau_powers_g1,
            tau_powers_g2,
            public_ic_terms,
            private_delta_terms,
            h_over_delta_powers,
            public_vars: public_vars.to_vec(),
            private_vars,
        }
    }
}

pub fn msm_g1(coeffs: &[Fr], bases: &[G1Projective]) -> G1Projective {
    let mut acc = G1Projective::zero();
    for (c, b) in coeffs.iter().zip(bases.iter()) {
        acc += *b * c;
    }
    acc
}

pub fn msm_g2(coeffs: &[Fr], bases: &[G2Projective]) -> G2Projective {
    let mut acc = G2Projective::zero();
    for (c, b) in coeffs.iter().zip(bases.iter()) {
        acc += *b * c;
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::test_rng;

    #[test]
    fn setup_produces_expected_structure() {
        let r1cs = R1cs::circuit();
        let qap = Qap::new(&r1cs, vec![Fr::from(5u64), Fr::from(7u64)]);
        let mut rng = test_rng();
        let tw = ToxicWaste::sample(&mut rng);
        let crs = Crs::setup(&r1cs, &qap, &tw, &[0, 1]);

        assert_eq!(crs.tau_powers_g1.len(), 2);
        assert_eq!(crs.tau_powers_g2.len(), 2);
        assert_eq!(crs.public_ic_terms.len(), 2);
        assert_eq!(crs.private_delta_terms.len(), 2);
        assert_eq!(crs.public_vars, vec![0, 1]);
        assert_eq!(crs.private_vars, vec![2, 3]);
    }

    #[test]
    fn tau_powers_g1_start_at_generator() {
        let r1cs = R1cs::circuit();
        let qap = Qap::new(&r1cs, vec![Fr::from(5u64), Fr::from(7u64)]);
        let mut rng = test_rng();
        let tw = ToxicWaste::sample(&mut rng);
        let crs = Crs::setup(&r1cs, &qap, &tw, &[0, 1]);

        assert_eq!(crs.tau_powers_g1[0], G1Projective::generator());
        assert_eq!(crs.tau_powers_g2[0], G2Projective::generator());
    }
}
