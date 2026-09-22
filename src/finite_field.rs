use ark_bn254::Fr;

pub fn field_ops() {
    let a = Fr::from(4u64); // x = 4
    let b = Fr::from(9u64); // c = 9
    let sum = a + b;
    let product = a * a; // s1 = x*x

    println!("a + b = {}", sum);
    println!("a * a = {}", product);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn x_equals_4_and_c_equals_9() {
        let x = Fr::from(4u64);
        let s1 = x * x;
        let y = s1 + x + Fr::from(9u64);

        assert_eq!(s1, Fr::from(16u64));
        assert_eq!(y, Fr::from(29u64));
    }
}
