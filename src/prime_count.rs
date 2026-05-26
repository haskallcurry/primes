use std::collections::HashMap;

pub trait PrimeCounter {
    fn pi(&mut self, x: u64) -> u64;
}

pub struct LehmerCounter {
    table: PrimeTable,
    phi_table: PhiTable,
    phi_cache: HashMap<(u64, usize), u64>,
    pi_cache: HashMap<u64, u64>,
}

impl LehmerCounter {
    pub fn new(lookup_limit: usize) -> Self {
        Self {
            table: PrimeTable::new(lookup_limit),
            phi_table: PhiTable::new(),
            phi_cache: HashMap::new(),
            pi_cache: HashMap::new(),
        }
    }

    fn lehmer_pi(&mut self, x: u64) -> u64 {
        if x < self.table.pi.len() as u64 {
            return self.table.pi[x as usize] as u64;
        }

        if let Some(&cached) = self.pi_cache.get(&x) {
            return cached;
        }

        self.ensure_table(isqrt(x) as usize + 1);

        let a = self.lehmer_pi(iroot(x, 4));
        let b = self.lehmer_pi(isqrt(x));
        let c = self.lehmer_pi(iroot(x, 3));
        let a = a as usize;
        let b = b as usize;
        let c = c as usize;

        let mut result = self.phi(x, a) + ((b as u64 + a as u64 - 2) * (b - a + 1) as u64) / 2;

        for i in a..b {
            let w = x / self.table.primes[i];
            result -= self.lehmer_pi(w);

            if i < c {
                let limit = self.lehmer_pi(isqrt(w)) as usize;
                for j in i..limit {
                    result -= self.lehmer_pi(w / self.table.primes[j]) - j as u64;
                }
            }
        }

        self.pi_cache.insert(x, result);
        result
    }

    fn phi(&mut self, x: u64, a: usize) -> u64 {
        if a <= self.phi_table.depth {
            return self.phi_table.get(x, a);
        }

        if let Some(&cached) = self.phi_cache.get(&(x, a)) {
            return cached;
        }

        let mut result = self.phi_table.get(x, self.phi_table.depth);
        for i in self.phi_table.depth..a {
            result -= self.phi(x / self.table.primes[i], i);
        }

        self.phi_cache.insert((x, a), result);
        result
    }

    fn ensure_table(&mut self, limit: usize) {
        if limit < self.table.pi.len() {
            return;
        }

        let current = self.table.pi.len() - 1;
        self.table = PrimeTable::new(limit.max(current * 2));
        self.pi_cache.clear();
    }
}

impl PrimeCounter for LehmerCounter {
    fn pi(&mut self, x: u64) -> u64 {
        self.lehmer_pi(x)
    }
}

struct PrimeTable {
    primes: Vec<u64>,
    pi: Vec<u32>,
}

impl PrimeTable {
    fn new(limit: usize) -> Self {
        let limit = limit.max(2);
        let is_prime = sieve(limit);
        let mut primes = Vec::new();
        let mut pi = vec![0u32; limit + 1];
        let mut count = 0u32;

        for n in 2..=limit {
            if is_prime[n] {
                count += 1;
                primes.push(n as u64);
            }
            pi[n] = count;
        }

        Self { primes, pi }
    }
}

fn sieve(limit: usize) -> Vec<bool> {
    let mut is_prime = vec![true; limit + 1];
    is_prime[0] = false;
    is_prime[1] = false;

    for p in 2..=isqrt(limit as u64) as usize {
        if is_prime[p] {
            let mut multiple = p * p;
            while multiple <= limit {
                is_prime[multiple] = false;
                multiple += p;
            }
        }
    }

    is_prime
}

struct PhiTable {
    depth: usize,
    modulus: u64,
    values: Vec<Vec<u64>>,
}

impl PhiTable {
    const PRIMES: [usize; 6] = [2, 3, 5, 7, 11, 13];

    fn new() -> Self {
        let depth = Self::PRIMES.len();
        let modulus = Self::PRIMES.iter().product::<usize>();
        let mut values = vec![vec![0; modulus + 1]; depth + 1];

        for (n, value) in values[0].iter_mut().enumerate() {
            *value = n as u64;
        }

        for a in 1..=depth {
            let p = Self::PRIMES[a - 1];
            let (previous_rows, current_rows) = values.split_at_mut(a);
            let previous = &previous_rows[a - 1];

            for (n, value) in current_rows[0].iter_mut().enumerate() {
                *value = previous[n] - previous[n / p];
            }
        }

        Self {
            depth,
            modulus: modulus as u64,
            values,
        }
    }

    fn get(&self, x: u64, a: usize) -> u64 {
        let blocks = x / self.modulus;
        let remainder = (x % self.modulus) as usize;
        blocks * self.values[a][self.modulus as usize] + self.values[a][remainder]
    }
}

fn isqrt(n: u64) -> u64 {
    iroot(n, 2)
}

fn iroot(n: u64, degree: u32) -> u64 {
    let mut root = (n as f64).powf(1.0 / f64::from(degree)).floor() as u64;

    while pow_leq(root.saturating_add(1), degree, n) {
        root += 1;
    }

    while !pow_leq(root, degree, n) {
        root -= 1;
    }

    root
}

fn pow_leq(base: u64, degree: u32, limit: u64) -> bool {
    let mut value = 1u64;

    for _ in 0..degree {
        let Some(next) = value.checked_mul(base) else {
            return false;
        };

        if next > limit {
            return false;
        }

        value = next;
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_prime_count_table() {
        let table = PrimeTable::new(20);

        assert_eq!(table.primes, vec![2, 3, 5, 7, 11, 13, 17, 19]);
        assert_eq!(table.pi[0], 0);
        assert_eq!(table.pi[1], 0);
        assert_eq!(table.pi[2], 1);
        assert_eq!(table.pi[3], 2);
        assert_eq!(table.pi[4], 2);
        assert_eq!(table.pi[20], 8);
    }

    #[test]
    fn counts_known_values() {
        let mut counter = LehmerCounter::new(10_000);

        assert_eq!(counter.pi(0), 0);
        assert_eq!(counter.pi(1), 0);
        assert_eq!(counter.pi(2), 1);
        assert_eq!(counter.pi(10), 4);
        assert_eq!(counter.pi(100), 25);
        assert_eq!(counter.pi(1_000), 168);
        assert_eq!(counter.pi(1_000_000), 78_498);
    }

    #[test]
    fn computes_roots_exactly() {
        assert_eq!(isqrt(0), 0);
        assert_eq!(isqrt(1), 1);
        assert_eq!(isqrt(2), 1);
        assert_eq!(isqrt(3), 1);
        assert_eq!(isqrt(4), 2);
        assert_eq!(isqrt(15), 3);
        assert_eq!(isqrt(16), 4);
        assert_eq!(isqrt(17), 4);
        assert_eq!(iroot(27, 3), 3);
        assert_eq!(iroot(28, 3), 3);
    }
}
