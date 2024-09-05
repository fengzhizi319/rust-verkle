// 允许非蛇形命名
#![allow(non_snake_case)]

// 导入所需的模块和trait
use crate::crs::CRS;
use crate::math_utils::inner_product;
use crate::transcript::{Transcript, TranscriptProtocol};

use banderwagon::{multi_scalar_mul, trait_defs::*, Element, Fr};
use itertools::Itertools;

use crate::{IOError, IOErrorKind, IOResult};

use std::iter;

// 定义IPAProof结构体
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IPAProof {
    // 这些字段现在是公开的，因为golang代码向客户端开发者公开了证明结构，
    // 如果我们不公开，那么我们就无法将json证明反序列化为IPAProof。
    pub L_vec: Vec<Element>,
    pub R_vec: Vec<Element>,
    pub a: Fr,
}

// 为IPAProof实现方法
impl IPAProof {
    // 返回序列化后的大小
    // 这个方法返回IPAProof对象序列化后的大小
    pub(crate) fn serialized_size(&self) -> usize {
        // 每个元素的大小为32字节，L_vec和R_vec中的元素总数为L_vec.len() * 2，再加上一个32字节的a
        (self.L_vec.len() * 2 + 1) * 32
    }

    // 从字节序列中创建IPAProof
    // 这个方法接收一个字节序列和一个多项式的度数，然后返回一个IPAProof对象
    pub fn from_bytes(bytes: &[u8], poly_degree: usize) -> IOResult<IPAProof> {
        // 给定多项式的度数，我们将有log2 * 2个点
        let num_points = log2(poly_degree);
        // 创建两个向量，用于存储L_vec和R_vec
        let mut L_vec = Vec::with_capacity(num_points as usize);
        let mut R_vec = Vec::with_capacity(num_points as usize);

        // 检查输入的字节序列的长度是否正确
        assert_eq!(((num_points * 2) + 1) * 32, bytes.len() as u32);
        // 检查输入的字节序列的长度是否是32的倍数
        assert!(bytes.len() % 32 == 0);

        // 将字节切片分块为32字节
        let mut chunks = bytes.chunks_exact(32);

        // 从字节序列中解析出L_vec和R_vec
        for _ in 0..num_points {
            let chunk = chunks.next().unwrap();
            let point: Element =
                Element::from_bytes(chunk).ok_or(IOError::from(IOErrorKind::InvalidData))?;
            L_vec.push(point)
        }

        for _ in 0..num_points {
            let chunk = chunks.next().unwrap();
            let point: Element =
                Element::from_bytes(chunk).ok_or(IOError::from(IOErrorKind::InvalidData))?;
            R_vec.push(point)
        }

        // 解析出a
        let last_32_bytes = chunks.next().unwrap();
        let a: Fr = CanonicalDeserialize::deserialize_compressed(last_32_bytes)
            .map_err(|_| IOError::from(IOErrorKind::InvalidData))?;

        // 返回一个新的IPAProof对象
        Ok(IPAProof { L_vec, R_vec, a })
    }

    // 将IPAProof转换为字节序列
    // 这个方法将IPAProof对象转换为一个字节序列
    pub fn to_bytes(&self) -> IOResult<Vec<u8>> {
        // 我们不序列化长度。我们假设反序列化器知道这个。
        let mut bytes = Vec::with_capacity(self.serialized_size());

        // 将L_vec和R_vec转换为字节序列
        for L in &self.L_vec {
            bytes.extend(L.to_bytes());
        }

        for R in &self.R_vec {
            bytes.extend(R.to_bytes());
        }

        // 将a转换为字节序列
        self.a
            .serialize_compressed(&mut bytes)
            .map_err(|_| IOError::from(IOErrorKind::InvalidData))?;
        // 返回字节序列
        Ok(bytes)
    }
}

// 创建IPAProof的函数
pub fn create(
    transcript: &mut Transcript, // 用于记录证明过程的Transcript对象
    mut crs: CRS, // 包含生成元的公共参考字符串
    mut a_vec: Vec<Fr>, // 第一个向量
    a_comm: Element, // 第一个向量的承诺
    mut b_vec: Vec<Fr>, // 第二个向量
    input_point: Fr, // 这是f(z)中的z
) -> IPAProof {
    // 在Transcript中添加一个域分隔符
    transcript.domain_sep(b"ipa");

    // 创建对向量和生成元的可变引用
    let mut a = &mut a_vec[..];
    let mut b = &mut b_vec[..];
    let mut G = &mut crs.G[..];

    // 获取向量的长度
    let n = G.len();

    // 所有输入向量必须具有相同的长度
    assert_eq!(G.len(), n);
    assert_eq!(a.len(), n);
    assert_eq!(b.len(), n);

    // 所有输入向量的长度必须是2的幂
    assert!(n.is_power_of_two());

    // 计算a和b的内积，并将其添加到Transcript中
    let output_point = inner_product(a, b);
    transcript.append_point(b"C", &a_comm);
    transcript.append_scalar(b"input point", &input_point);
    transcript.append_scalar(b"output point", &output_point);

    // 从Transcript中获取一个挑战标量w，并计算Q = crs.Q * w
    let w = transcript.challenge_scalar(b"w");
    let Q = crs.Q * w;

    // 计算需要进行的轮数
    let num_rounds = log2(n);

    // 创建两个向量，用于存储每一轮的L和R
    let mut L_vec: Vec<Element> = Vec::with_capacity(num_rounds as usize);
    let mut R_vec: Vec<Element> = Vec::with_capacity(num_rounds as usize);

    // 进行每一轮的计算
    //let start = std::time::Instant::now();
    for _k in 0..num_rounds {
        // 将a，b，G分为两半

        let (a_L, a_R) = halve(a);
        let (b_L, b_R) = halve(b);
        let (G_L, G_R) = halve(G);

        // 计算z_L和z_R
        let z_L = inner_product(a_R, b_L);
        let z_R = inner_product(a_L, b_R);

        // 计算L和R
        //let start1 = std::time::Instant::now();
        let L = slow_vartime_multiscalar_mul(
            a_R.iter().chain(iter::once(&z_L)),
            G_L.iter().chain(iter::once(&Q)),
        );
        let R = slow_vartime_multiscalar_mul(
            a_L.iter().chain(iter::once(&z_R)),
            G_R.iter().chain(iter::once(&Q)),
        );
        // let end1 = std::time::Instant::now();
        //println!("折半 time: {:?}", end1.duration_since(start1));

        // 将L和R添加到向量中，并将它们添加到Transcript中
        L_vec.push(L);
        R_vec.push(R);
        transcript.append_point(b"L", &L);
        transcript.append_point(b"R", &R);

        // 从Transcript中获取一个挑战标量x，并计算其逆元x_inv
        let x = transcript.challenge_scalar(b"x");
        let x_inv = x.inverse().unwrap();

        // 更新a，b，G
        for i in 0..a_L.len() {
            a_L[i] += x * a_R[i];
            b_L[i] += x_inv * b_R[i];
            G_L[i] += G_R[i] * x_inv;
        }

        a = a_L;
        b = b_L;
        G = G_L;
    }
    // let end = std::time::Instant::now();
    // println!("num_rounds time: {:?}", end.duration_since(start));

    // 返回一个新的IPAProof对象
    IPAProof {
        L_vec,
        R_vec,
        a: a[0],
    }
}

// `halve`函数接收一个可变切片，并将其分为两半
// 这个函数假设输入的切片长度是偶数
fn halve<T>(scalars: &mut [T]) -> (&mut [T], &mut [T]) {
    // 获取切片的长度
    let len = scalars.len();
    // 将切片分为两半，并返回两个新的切片
    scalars.split_at_mut(len / 2)
}

// `log2`函数接收一个usize类型的数，并返回它的下一个2的幂的二进制表示中0的个数
// 这个函数可以用来计算一个数的对数，以2为底
fn log2(n: usize) -> u32 {
    // 计算n的下一个2的幂，然后返回它的二进制表示中0的个数
    n.next_power_of_two().trailing_zeros()
}

impl IPAProof {
    // 验证IPAProof的函数
    pub fn verify(
        &self, // 当前的IPAProof对象
        transcript: &mut Transcript, // 用于记录证明过程的Transcript对象
        mut crs: CRS, // 包含生成元的公共参考字符串
        mut b: Vec<Fr>, // 向量b
        a_comm: Element, // 向量a的承诺
        input_point: Fr, // 这是f(z)中的z
        output_point: Fr, // 这是f(z)
    ) -> bool {
        // 在Transcript中添加一个域分隔符
        transcript.domain_sep(b"ipa");

        // 创建对向量和生成元的可变引用
        let mut G = &mut crs.G[..];
        let mut b = &mut b[..];

        // 获取需要进行的轮数
        let num_rounds = self.L_vec.len();

        // 检查证明者是否计算了一个大小为n的向量的内部证明
        if crs.n != 1 << num_rounds {
            return false;
        }

        // 将a的承诺，输入点和输出点添加到Transcript中
        transcript.append_point(b"C", &a_comm);
        transcript.append_scalar(b"input point", &input_point);
        transcript.append_scalar(b"output point", &output_point);

        // 从Transcript中获取一个挑战标量w，并计算Q = crs.Q * w
        let w = transcript.challenge_scalar(b"w");
        let Q = crs.Q * w;

        // 计算a的承诺
        let mut a_comm = a_comm + (Q * output_point);

        // 生成挑战标量
        let challenges = generate_challenges(self, transcript);
        let mut challenges_inv = challenges.clone();
        batch_inversion(&mut challenges_inv);

        // 计算预期的承诺
        for i in 0..num_rounds {
            let x = challenges[i];
            let x_inv = challenges_inv[i];
            let L = self.L_vec[i];
            let R = self.R_vec[i];

            a_comm = a_comm + (L * x) + (R * x_inv);
        }

        // 更新G和b
        for x_inv in challenges_inv.iter() {
            let (G_L, G_R) = halve(G);
            let (b_L, b_R) = halve(b);

            for i in 0..G_L.len() {
                G_L[i] += G_R[i] * *x_inv;
                b_L[i] += b_R[i] * x_inv;
            }
            G = G_L;
            b = b_L;
        }

        // 检查G和b的长度是否为1
        assert_eq!(G.len(), 1);
        assert_eq!(b.len(), 1);

        // 计算预期的P
        let exp_P = (G[0] * self.a) + Q * (self.a * b[0]);

        // 检查预期的P是否等于a的承诺
        exp_P == a_comm
    }
    // 验证多重指数的函数
    pub fn verify_multiexp(
        &self, // 当前的IPAProof对象
        transcript: &mut Transcript, // 用于记录证明过程的Transcript对象
        crs: &CRS, // 包含生成元的公共参考字符串
        b_vec: Vec<Fr>, // 向量b
        a_comm: Element, // 向量a的承诺
        input_point: Fr, // 这是f(z)中的z
        output_point: Fr, // 这是f(z)
    ) -> bool {
        // 在Transcript中添加一个域分隔符
        transcript.domain_sep(b"ipa");

        // 获取需要进行的轮数
        let logn = self.L_vec.len();
        let n = crs.n;

        // 检查证明者是否计算了一个大小为n的向量的内部证明
        if n != (1 << logn) {
            return false;
        }

        // 将a的承诺，输入点和输出点添加到Transcript中
        transcript.append_point(b"C", &a_comm);
        transcript.append_scalar(b"input point", &input_point);
        transcript.append_scalar(b"output point", &output_point);

        // 计算将增加到内部产品对应点的标量
        let w = transcript.challenge_scalar(b"w");

        // 生成所有必要的挑战和它们的逆
        let challenges = generate_challenges(self, transcript);
        let mut challenges_inv = challenges.clone();
        batch_inversion(&mut challenges_inv);

        // 生成`G`向量和`b`向量的系数
        let mut g_i: Vec<Fr> = Vec::with_capacity(1 << logn);
        let mut b_i: Vec<Fr> = Vec::with_capacity(1 << logn);

        // 计算g_i和b_i
        for index in 0..n {
            let mut b = -Fr::one();
            for (bit, x_inv) in to_bits(index, logn).zip_eq(&challenges_inv) {
                if bit == 1 {
                    b *= x_inv;
                }
            }
            b_i.push(b);
            g_i.push(self.a * b);
        }

        // 计算b_0和q_i
        let b_0 = inner_product(&b_vec, &b_i);
        let q_i = w * (output_point + self.a * b_0);

        // 执行慢速的多标量乘法
        slow_vartime_multiscalar_mul(
            challenges
                .iter()
                .chain(challenges_inv.iter())
                .chain(iter::once(&Fr::one()))
                .chain(iter::once(&q_i))
                .chain(g_i.iter()),
            self.L_vec
                .iter()
                .chain(self.R_vec.iter())
                .chain(iter::once(&a_comm))
                .chain(iter::once(&crs.Q))
                .chain(crs.G.iter()),
        )
            .is_zero() // 检查结果是否为零
    }
    // It's only semi unrolled.
    // This is being committed incase someone goes through the git history
    // The fully unrolled code is not that intuitive, but maybe this semi
    // unrolled version can help you to figure out the gap
    // 验证半多重指数的函数
    pub fn verify_semi_multiexp(
        &self, // 当前的IPAProof对象
        transcript: &mut Transcript, // 用于记录证明过程的Transcript对象
        crs: &CRS, // 包含生成元的公共参考字符串
        b_Vec: Vec<Fr>, // 向量b
        a_comm: Element, // 向量a的承诺
        input_point: Fr, // 这是f(z)中的z
        output_point: Fr, // 这是f(z)
    ) -> bool {
        // 在Transcript中添加一个域分隔符
        transcript.domain_sep(b"ipa");

        // 获取需要进行的轮数
        let logn = self.L_vec.len();
        let n = crs.n;

        // 检查证明者是否计算了一个大小为n的向量的内部证明
        if n != (1 << logn) {
            return false;
        }

        // 将a的承诺，输入点和输出点添加到Transcript中
        transcript.append_point(b"C", &a_comm);
        transcript.append_scalar(b"input point", &input_point);
        transcript.append_scalar(b"output point", &output_point);

        // 计算将增加到内部产品对应点的标量
        let w = transcript.challenge_scalar(b"w");
        let Q = crs.Q * w;

        // 计算a的承诺
        let a_comm = a_comm + (Q * output_point);

        // 生成所有必要的挑战和它们的逆
        let challenges = generate_challenges(self, transcript);
        let mut challenges_inv = challenges.clone();
        batch_inversion(&mut challenges_inv);

        // 计算P
        let P = slow_vartime_multiscalar_mul(
            challenges
                .iter()
                .chain(challenges_inv.iter())
                .chain(iter::once(&Fr::one())),
            self.L_vec
                .iter()
                .chain(self.R_vec.iter())
                .chain(iter::once(&a_comm)),
        );

        // 计算g_i
        let mut g_i: Vec<Fr> = Vec::with_capacity(1 << logn);
        for index in 0..n {
            let mut g = Fr::one();
            for (bit, x_inv) in to_bits(index, logn).zip_eq(&challenges_inv) {
                if bit == 1 {
                    g *= x_inv;
                }
            }
            g_i.push(g);
        }

        // 计算b_0和G_0
        let b_0 = inner_product(&b_Vec, &g_i);
        let G_0 = slow_vartime_multiscalar_mul(g_i.iter(), crs.G.iter()); // TODO: Optimise; the majority of the time is spent on this vector, precompute

        // 计算预期的P
        let exp_P = (G_0 * self.a) + Q * (self.a * b_0);

        // 检查预期的P是否等于P
        exp_P == P
    }
}

// 将一个数字转换为二进制位的迭代器
fn to_bits(n: usize, bits_needed: usize) -> impl Iterator<Item=u8> {
    // 从0到需要的位数，将数字右移i位并与1进行与运算，然后反转结果
    (0..bits_needed).map(move |i| ((n >> i) & 1) as u8).rev()
}

// 慢速的多标量乘法函数
pub fn slow_vartime_multiscalar_mul<'a>(
    scalars: impl Iterator<Item=&'a Fr>, // 标量的迭代器
    points: impl Iterator<Item=&'a Element>, // 点的迭代器
) -> Element {
    // 将标量和点的迭代器转换为向量
    let scalars: Vec<_> = scalars.into_iter().copied().collect();
    let points: Vec<_> = points.into_iter().copied().collect();
    // 调用multi_scalar_mul函数进行多标量乘法
    multi_scalar_mul(&points, &scalars)
}

// 生成挑战的函数
fn generate_challenges(proof: &IPAProof, transcript: &mut Transcript) -> Vec<Fr> {
    // 创建一个容量为L_vec长度的向量
    let mut challenges: Vec<Fr> = Vec::with_capacity(proof.L_vec.len());

    // 对于L_vec和R_vec中的每一对元素，将它们添加到Transcript中，并生成一个挑战标量x_i，然后将x_i添加到挑战向量中
    for (L, R) in proof.L_vec.iter().zip(proof.R_vec.iter()) {
        transcript.append_point(b"L", L);
        transcript.append_point(b"R", R);

        let x_i = transcript.challenge_scalar(b"x");
        challenges.push(x_i);
    }

    // 返回挑战向量
    challenges
}
//参数polynomial设置为数组类型，返回值为Vec<Fr>类型


pub(crate) fn genetestpoly256(polynomial: &[i32]) -> Vec<Fr> {
    let mut poly: Vec<Fr> = Vec::with_capacity(256);
    let n = polynomial.len();
    //从0到n，将polynomial中的数字转换为Fr类型并添加到向量中
    for i in 0..256 {
        if i < n {
            poly.push(Fr::from(polynomial[i]));
        } else {
            poly.push(Fr::zero());
        }
    }
    poly
}



#[cfg(test)]
mod tests {
    use super::*;
    use crate::crs::CRS;
    use crate::math_utils::{inner_product, powers_of};

    use ark_std::{rand::SeedableRng, UniformRand};
    use rand_chacha::ChaCha20Rng;
    use crate::lagrange_basis::{LagrangeBasis, PrecomputedWeights};

    #[test]
    // 测试创建和验证IPAProof的函数
    fn test_create_IPAProof_CreateVerify() {
        let point = Fr::from(123456789);

        // 设置向量的长度
        let n = 256;
        // 创建一个新的CRS对象
        let crs = CRS::new(n, b"random seed");

        // 创建一个随机的Fr向量a
        let a1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
        //将a1转换为vec
        let poly = genetestpoly256(&a1);
        let mut prover_transcript = Transcript::new(b"ipa");
        // 计算承诺，它是向量a和G的多标量乘积
        let prover_comm = slow_vartime_multiscalar_mul(poly.iter(), crs.G.iter());
        // 创建一个新的IPAProof对象


        let precomp = PrecomputedWeights::new(256);
        let lagrange_coeffs =
            LagrangeBasis::evaluate_lagrange_coefficients(&precomp, 256, point);
        let inner_product = inner_product(&poly, &lagrange_coeffs);
        let proof = create(
            &mut prover_transcript,
            crs.clone(),
            poly,
            prover_comm,
            lagrange_coeffs.clone(),
            point,
        );


        // 创建一个新的Transcript对象，用于记录验证过程
        let mut verifier_transcript = Transcript::new(b"ipa");
        // 验证IPAProof对象，如果验证成功，断言将通过
        assert!(proof.verify(
            &mut verifier_transcript,
            crs,
            lagrange_coeffs,
            prover_comm,
            point,
            inner_product,
        ));
    }

    #[test]
    fn test_create_IPAProof_ConsistencySimpleProof() {
        let point = Fr::from(2101);

        // 设置向量的长度
        let n = 256;
        // 创建一个新的CRS对象
        let crs = CRS::new(n, b"random seed");

        // 创建一个随机的Fr向量a
        let a1 = [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
        ];
        //将a1转换为vec
        let poly = genetestpoly256(&a1);
        let mut prover_transcript = Transcript::new(b"ipa");
        // 计算承诺，它是向量a和G的多标量乘积
        let prover_comm = slow_vartime_multiscalar_mul(poly.iter(), crs.G.iter());
        //1b9dff8f5ebbac250d291dfe90e36283a227c64b113c37f1bfb9e7a743cdb128
        //用十六进制打印prover_comm

        println!("prover_comm: {:x?}", prover_comm);
        // 创建一个新的IPAProof对象


        let precomp = PrecomputedWeights::new(256);
        let lagrange_coeffs =
            LagrangeBasis::evaluate_lagrange_coefficients(&precomp, 256, point);
        let inner_product = inner_product(&poly, &lagrange_coeffs);
        //4a353e70b03c89f161de002e8713beec0d740a5e20722fd5bd68b30540a33208
        println!("inner_product: {:x?}", inner_product);
        let proof = create(
            &mut prover_transcript,
            crs.clone(),
            poly,
            prover_comm,
            lagrange_coeffs.clone(),
            point,
        );


        // 创建一个新的Transcript对象，用于记录验证过程
        let mut verifier_transcript = Transcript::new(b"ipa");
        // 验证IPAProof对象，如果验证成功，断言将通过
        assert!(proof.verify(
            &mut verifier_transcript,
            crs,
            lagrange_coeffs,
            prover_comm,
            point,
            inner_product,
        ));
    }

    #[test]
    fn test_create_IPAProof_proof() {
        // 设置向量的长度
        let n = 8;
        // 创建一个新的CRS对象
        let crs = CRS::new(n, b"random seed");

        // 创建一个随机数生成器
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        // 创建一个随机的Fr向量a
        let a: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        // 创建一个随机的输入点
        let input_point = Fr::rand(&mut rng);

        // 计算向量b，它是输入点的幂
        let b = powers_of(input_point, n);
        // 计算输出点，它是向量a和b的内积
        let output_point = inner_product(&a, &b);

        // 创建一个新的Transcript对象，用于记录证明过程
        let mut prover_transcript = Transcript::new(b"ip_no_zk");

        // 计算P，它是向量a和G的多标量乘积
        let P = slow_vartime_multiscalar_mul(a.iter(), crs.G.iter());

        // 创建一个新的IPAProof对象
        let proof = create(
            &mut prover_transcript,
            crs.clone(),
            a,
            P,
            b.clone(),
            input_point,
        );

        // 创建一个新的Transcript对象，用于记录验证过程
        let mut verifier_transcript = Transcript::new(b"ip_no_zk");
        // 验证IPAProof对象，如果验证成功，断言将通过
        assert!(proof.verify(
            &mut verifier_transcript,
            crs,
            b,
            P,
            input_point,
            output_point,
        ));
    }
    #[test]
    #[allow(unused)]
    fn test_Commit_Time() {
        // let guard = pprof::ProfilerGuardBuilder::default()
        //     .frequency(1000)
        //     .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        //     .build()
        //     .unwrap();

        let n = 256;
        // 创建一个新的CRS对象
        let crs = CRS::new(n, b"random seed");

        // 创建一个随机数生成器
        let mut rng = ChaCha20Rng::from_seed([12u8; 32]);
        // 创建一个随机的Fr向量a
        let a: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        //测试函数运行时间
        let start = std::time::Instant::now();
        for j in 0..500 {
            let a: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
            let P = slow_vartime_multiscalar_mul(a.iter(), crs.G.iter());
        }

        let end = std::time::Instant::now();
        println!("Function took this amount of time:: {:?}", end.duration_since(start)/500);
        // if let Ok(report) = guard.report().build() {
        //     println!("report: {:?}", &report);
        //     let file = std::fs::File::create("flamegraph.svg").unwrap();
        //     report.flamegraph(file).unwrap();
        // }
    }
    #[test]
    fn test_create_IPAProof_proof256() {
        // 设置向量的长度
        let n = 256;
        // 创建一个新的CRS对象
        let crs = CRS::new(n, b"random seed");

        // 创建一个随机数生成器
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        // 创建一个随机的Fr向量a
        let a: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        // 创建一个随机的输入点
        let input_point = Fr::rand(&mut rng);

        // 计算向量b，它是输入点的幂
        let b = powers_of(input_point, n);
        // 计算输出点，它是向量a和b的内积
        let output_point = inner_product(&a, &b);

        // 创建一个新的Transcript对象，用于记录证明过程
        let mut prover_transcript = Transcript::new(b"ip_no_zk");

        // 计算P，它是向量a和G的多标量乘积
        //测试函数运行时间
        //let start = std::time::Instant::now();
        let P = slow_vartime_multiscalar_mul(a.iter(), crs.G.iter());

        let start = std::time::Instant::now();
        // 创建一个新的IPAProof对象
        let proof = create(
            &mut prover_transcript,
            crs.clone(),
            a,
            P,
            b.clone(),
            input_point,
        );
        let end = std::time::Instant::now();
        println!("create time: {:?}", end.duration_since(start));

        // 创建一个新的Transcript对象，用于记录验证过程
        let mut verifier_transcript = Transcript::new(b"ip_no_zk");
        // 验证IPAProof对象，如果验证成功，断言将通过
        assert!(proof.verify(
            &mut verifier_transcript,
            crs,
            b,
            P,
            input_point,
            output_point,
        ));
    }
}
