use crate::polynomial::Poly;
use crate::trusted_setup::{msm_g1, msm_g2, Crs};
use ark_bn254::{Fr, G1Projective, G2Projective};
use ark_ff::{UniformRand, Zero};
use ark_std::rand::RngCore;

pub struct Proof {
    pub a: G1Projective,
    pub b: G2Projective,
    pub c: G1Projective,
}

pub fn prove<R: RngCore>(
    crs: &Crs,
    w: &[Fr],
    a_poly: &Poly,
    b_poly: &Poly,
    h_poly: &Poly,
    rng: &mut R,
) -> Proof {
    let r = Fr::rand(rng);
    let s = Fr::rand(rng);

    let a_tau_g1 = msm_g1(&a_poly.coeffs, &crs.tau_powers_g1);
    let b_tau_g1 = msm_g1(&b_poly.coeffs, &crs.tau_powers_g1);
    let b_tau_g2 = msm_g2(&b_poly.coeffs, &crs.tau_powers_g2);

    let pi_a = crs.alpha_g1 + a_tau_g1 + crs.delta_g1 * r;
    let pi_b_g2 = crs.beta_g2 + b_tau_g2 + crs.delta_g2 * s;
    let b_g1 = crs.beta_g1 + b_tau_g1 + crs.delta_g1 * s;

    let private_sum: G1Projective = crs
        .private_vars
        .iter()
        .zip(crs.private_delta_terms.iter())
        .fold(G1Projective::zero(), |acc, (&j, &term)| acc + term * w[j]);

    let h_term = msm_g1(&h_poly.coeffs, &crs.h_over_delta_powers);

    let pi_c = private_sum + h_term + pi_a * s + b_g1 * r - crs.delta_g1 * (r * s);

    Proof {
        a: pi_a,
        b: pi_b_g2,
        c: pi_c,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qap::Qap;
    use crate::r1cs::R1cs;
    use crate::trusted_setup::ToxicWaste;
    use ark_std::test_rng;

    fn build_proof_inputs() -> (Crs, Vec<Fr>, Poly, Poly, Poly) {
        let r1cs = R1cs::circuit();
        let w = r1cs.build_witness(4);
        let qap = Qap::new(&r1cs, vec![Fr::from(5u64), Fr::from(7u64)]);
        let (a_poly, b_poly, c_poly) = qap.build_abc(&r1cs, &w);
        let h_poly = qap.compute_h(&a_poly, &b_poly, &c_poly);

        let mut rng = test_rng();
        let tw = ToxicWaste::sample(&mut rng);
        let crs = Crs::setup(&r1cs, &qap, &tw, &[0, 1]);

        (crs, w, a_poly, b_poly, h_poly)
    }

    #[test]
    fn produces_three_distinct_nonzero_points() {
        let (crs, w, a_poly, b_poly, h_poly) = build_proof_inputs();
        let mut rng = test_rng();

        let proof = prove(&crs, &w, &a_poly, &b_poly, &h_poly, &mut rng);

        assert_ne!(proof.a, G1Projective::zero());
        assert_ne!(proof.b, G2Projective::zero());
        assert_ne!(proof.c, G1Projective::zero());
    }

    #[test]
    fn two_proofs_of_the_same_witness_differ() {
        let (crs, w, a_poly, b_poly, h_poly) = build_proof_inputs();
        let mut rng = test_rng();

        let proof1 = prove(&crs, &w, &a_poly, &b_poly, &h_poly, &mut rng);
        let proof2 = prove(&crs, &w, &a_poly, &b_poly, &h_poly, &mut rng);

        assert_ne!(proof1.a, proof2.a);
        assert_ne!(proof1.c, proof2.c);
    }
}