// 引入所需的库和模块
use banderwagon::{trait_defs::*, Fr};

/// 计算两个标量向量之间的内积
pub fn inner_product(a: &[Fr], b: &[Fr]) -> Fr {
    // 使用zip将两个向量配对，然后使用map将每对元素相乘，最后使用sum求和
    a.iter().zip(b.iter()).map(|(a, b)| *a * *b).sum()
}

// 计算一个点的幂次
pub fn powers_of(point: Fr, n: usize) -> Vec<Fr> {
    // 创建一个容量为n的向量
    let mut powers = Vec::with_capacity(n);
    // 将1添加到向量中
    powers.push(Fr::one());

    // 对于1到n的每个i，将前一个元素乘以点添加到向量中
    for i in 1..n {
        powers.push(powers[i - 1] * point);
    }
    // 返回向量
    powers
}

// 单元测试
#[test]
fn simple_vandemonde() {
    // 引入所需的库和模块
    use ark_std::test_rng;
    use ark_std::UniformRand;

    // 生成一个随机的Fr
    let rand_fr = Fr::rand(&mut test_rng());
    // 定义n为100
    let n = 100;
    // 计算rand_fr的幂次
    let powers = powers_of(rand_fr, n);

    // 断言第一个元素为1
    assert_eq!(powers[0], Fr::one());
    // 断言最后一个元素为rand_fr的n-1次幂
    assert_eq!(powers[n - 1], rand_fr.pow([(n - 1) as u64]));

    // 遍历向量，断言每个元素为rand_fr的i次幂
    for (i, power) in powers.into_iter().enumerate() {
        assert_eq!(power, rand_fr.pow([i as u64]))
    }
}