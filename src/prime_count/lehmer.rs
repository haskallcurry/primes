use std::collections::HashMap;

use super::{PhiTable, PrimeCounter, PrimeTable, iroot, isqrt};

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

#[cfg(test)]
mod tests {
    use super::*;

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
        assert_eq!(counter.pi(10_000_000), 664_579);
    }
}
