use crate::proof::Proof;
use crate::trusted_setup::Crs;
use ark_bn254::{Bn254, Fr, G1Projective};
use ark_ec::pairing::Pairing;
use ark_ff::Zero;

pub fn verify(crs: &Crs, proof: &Proof, public_values: &[Fr]) -> bool {
    assert_eq!(public_values.len(), crs.public_ic_terms.len());

    let ic: G1Projective = public_values
        .iter()
        .zip(crs.public_ic_terms.iter())
        .fold(G1Projective::zero(), |acc, (&v, &term)| acc + term * v);

    let lhs = Bn254::pairing(proof.a, proof.b);
    let rhs = Bn254::pairing(crs.alpha_g1, crs.beta_g2)
        + Bn254::pairing(ic, crs.gamma_g2)
        + Bn254::pairing(proof.c, crs.delta_g2);

    lhs == rhs
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::proof::prove;
    use crate::qap::Qap;
    use crate::r1cs::R1cs;
    use crate::trusted_setup::ToxicWaste;
    use ark_std::test_rng;

    fn build_proof_for_x(x: i64) -> (Crs, Proof, Vec<Fr>) {
        let r1cs = R1cs::circuit();
        let w = r1cs.build_witness(x);
        let qap = Qap::new(&r1cs, vec![Fr::from(5u64), Fr::from(7u64)]);
        let (a_poly, b_poly, c_poly) = qap.build_abc(&r1cs, &w);
        let h_poly = qap.compute_h(&a_poly, &b_poly, &c_poly);

        let mut rng = test_rng();
        let tw = ToxicWaste::sample(&mut rng);
        let crs = Crs::setup(&r1cs, &qap, &tw, &[0, 1]);
        let proof = prove(&crs, &w, &a_poly, &b_poly, &h_poly, &mut rng);

        let public_values = vec![w[0], w[1]];
        (crs, proof, public_values)
    }

    #[test]
    fn verifies_valid_proof_for_worked_witness() {
        let (crs, proof, public_values) = build_proof_for_x(4);
        assert!(verify(&crs, &proof, &public_values));
    }

    #[test]
    fn rejects_wrong_public_input() {
        let (crs, proof, _) = build_proof_for_x(4);
        let wrong = vec![Fr::from(1u64), Fr::from(999u64)];
        assert!(!verify(&crs, &proof, &wrong));
    }

    #[test]
    fn rejects_tampered_proof() {
        let (crs, mut proof, public_values) = build_proof_for_x(4);
        proof.a += crs.alpha_g1;
        assert!(!verify(&crs, &proof, &public_values));
    }

    #[test]
    fn end_to_end_across_several_witnesses() {
        for x in [0i64, 1, 4, 10, 25] {
            let (crs, proof, public_values) = build_proof_for_x(x);
            assert!(verify(&crs, &proof, &public_values), "failed for x={}", x);
        }
    }
}
