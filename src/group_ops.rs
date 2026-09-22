use ark_bn254::{Fr, G1Projective, G2Projective};
use ark_ec::Group;

pub fn group_ops() {
    let g1 = G1Projective::generator();

    let x = Fr::from(4u64);
    let x_g1 = g1 * x; // scalar multiplication (this is "[4]_1)"

    let sum = x_g1 + g1; // point addition

    println!("[4]_1 = {:?}", x_g1);
    println!("sum   = {:?}", sum);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalar_mult_is_repeated_addition() {
        let g1 = G1Projective::generator();
        let via_scalar = g1 * Fr::from(3u64);
        let via_addition = g1 + g1 + g1;
        assert_eq!(via_scalar, via_addition);
    }

    #[test]
    fn g1_and_g2_are_different_groups() {
        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();
        assert_ne!(g1, G1Projective::default());
        assert_ne!(g2, G2Projective::default());
    }

    #[test]
    fn bilinearity_holds_under_real_pairing() {
        use ark_bn254::Bn254;
        use ark_ec::pairing::Pairing;

        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();
        let a = Fr::from(6u64);
        let b = Fr::from(7u64);

        let lhs = Bn254::pairing(g1 * a, g2 * b);
        let rhs = Bn254::pairing(g1, g2) * (a * b);
        assert_eq!(lhs.0, rhs.0);
    }
}
