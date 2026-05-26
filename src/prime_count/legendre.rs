use std::collections::HashMap;

use super::{PhiTable, PrimeCounter, PrimeTable, isqrt};

pub struct LegendreCounter {
    table: PrimeTable,
    phi_table: PhiTable,
    pi_cache: HashMap<u64, u64>,
}

impl LegendreCounter {
    pub fn new(lookup_limit: usize) -> Self {
        Self {
            table: PrimeTable::new(lookup_limit),
            phi_table: PhiTable::new(),
            pi_cache: HashMap::new(),
        }
    }

    fn legendre_pi(&mut self, x: u64) -> u64 {
        if x < self.table.pi.len() as u64 {
            return self.table.pi[x as usize] as u64;
        }

        if let Some(&cached) = self.pi_cache.get(&x) {
            return cached;
        }

        let root = isqrt(x);
        self.ensure_table(root as usize + 1);

        let a = self.legendre_pi(root) as usize;
        let result = self.phi(x, a) + a as u64 - 1;
        self.pi_cache.insert(x, result);
        result
    }

    fn phi(&mut self, x: u64, a: usize) -> u64 {
        if a <= self.phi_table.depth {
            return self.phi_table.get(x, a);
        }

        self.phi_compute(x, a)
    }

    fn phi_compute(&mut self, x: u64, a: usize) -> u64 {
        let mut count = x;

        for k in 0..a {
            let p = self.table.primes[k];
            let quotient = x / p;
            if quotient == 0 {
                return count;
            }

            let removed = self.phi(quotient, k);
            if removed == 1 {
                // After this point each remaining prime still below x removes exactly one value.
                count -= 1;
                for l in k + 1..a {
                    if self.table.primes[l] > x {
                        break;
                    }
                    count -= 1;
                }
                return count;
            }

            count -= removed;
        }

        count
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

impl PrimeCounter for LegendreCounter {
    fn pi(&mut self, x: u64) -> u64 {
        self.legendre_pi(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_known_values() {
        let mut counter = LegendreCounter::new(10_000);

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
