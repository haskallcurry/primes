use std::env;
use std::process;
use std::time::Instant;

mod prime_count;

use prime_count::{DefaultCounter, PrimeCounter, isqrt};

const DEFAULT_TARGET: u64 = 1_000_000_000;
const DEFAULT_LOOKUP_LIMIT: usize = 10_000_000;

fn main() {
    let command = Command::from_args(env::args().skip(1)).unwrap_or_else(|error| {
        eprintln!("{error}");
        process::exit(2);
    });

    let started = Instant::now();
    let mut counter = DefaultCounter::new(DEFAULT_LOOKUP_LIMIT);

    match command {
        Command::Count(x) => {
            let count = counter.pi(x);
            eprintln!("pi({x}) = {count} in {:.3?}", started.elapsed());
            println!("{count}");
        }
        Command::Nth(target) => {
            let result = nth_prime(target, &mut counter);
            eprintln!(
                "li^-1({target}) estimate = {}; bracket=[{}, {}]",
                result.estimate, result.lower, result.upper
            );
            eprintln!(
                "prime #{target} = {} in {:.3?}",
                result.prime,
                started.elapsed()
            );
            println!("{}", result.prime);
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Command {
    Count(u64),
    Nth(u64),
}

impl Command {
    fn from_args(args: impl IntoIterator<Item = String>) -> Result<Self, String> {
        let mut mode = None;
        let mut target = DEFAULT_TARGET;
        let mut count_x = None;

        for arg in args {
            if arg == "-h" || arg == "--help" {
                return Err(usage());
            }

            let Some((name, value)) = arg.split_once('=') else {
                return Err(format!("unexpected argument '{arg}'\n\n{}", usage()));
            };

            match name {
                "--count" | "--x" => {
                    mode = Some(CommandMode::Count);
                    count_x = Some(parse_positive_u64(name, value)?);
                }
                "--target" => {
                    mode.get_or_insert(CommandMode::Nth);
                    target = parse_positive_u64(name, value)?;
                }
                "--mode" => {
                    mode = Some(CommandMode::parse(value)?);
                }
                _ => return Err(format!("unknown option '{name}'\n\n{}", usage())),
            }
        }

        Ok(match mode.unwrap_or(CommandMode::Nth) {
            CommandMode::Count => Self::Count(count_x.unwrap_or(target)),
            CommandMode::Nth => Self::Nth(target),
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CommandMode {
    Count,
    Nth,
}

impl CommandMode {
    fn parse(value: &str) -> Result<Self, String> {
        match value {
            "count" | "pi" => Ok(Self::Count),
            "invert" | "nth" | "nth-pi" => Ok(Self::Nth),
            _ => Err(format!(
                "invalid value for --mode: expected 'count' or 'nth', got '{value}'"
            )),
        }
    }
}

fn parse_positive_u64(name: &str, value: &str) -> Result<u64, String> {
    let parsed = value
        .parse::<u64>()
        .map_err(|error| format!("invalid value for {name}: {error}"))?;

    if parsed == 0 {
        return Err(format!("{name} must be positive"));
    }

    Ok(parsed)
}

fn usage() -> String {
    format!(
        "usage: billionprimes [--target=N | --count=N]\n\
         defaults: --target={DEFAULT_TARGET}"
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct NthPrimeResult {
    prime: u64,
    estimate: u64,
    lower: u64,
    upper: u64,
}

fn nth_prime<C: PrimeCounter + ?Sized>(target: u64, counter: &mut C) -> NthPrimeResult {
    if target == 1 {
        return NthPrimeResult {
            prime: 2,
            estimate: 2,
            lower: 2,
            upper: 2,
        };
    }

    let estimate = inverse_logarithmic_integral(target);
    let mut step = isqrt(estimate).saturating_mul(2).max(1024);
    let (lower, upper) = if counter.pi(estimate) < target {
        let lower = estimate;
        let mut upper = estimate.saturating_add(step).max(3);
        while counter.pi(upper) < target {
            upper = upper.saturating_add(step);
            step = step.saturating_mul(2);
        }
        (lower, upper)
    } else {
        let upper = estimate;
        let mut lower = estimate.saturating_sub(step).max(2);
        while counter.pi(lower) >= target {
            let next_lower = lower.saturating_sub(step).max(2);
            if next_lower == lower {
                break;
            }
            lower = next_lower;
            step = step.saturating_mul(2);
        }
        (lower, upper)
    };

    let bracket = (lower, upper);
    let mut left = lower + 1;
    let mut right = upper;

    while left < right {
        let mid = left + (right - left) / 2;
        if counter.pi(mid) >= target {
            right = mid;
        } else {
            left = mid + 1;
        }
    }

    NthPrimeResult {
        prime: left,
        estimate,
        lower: bracket.0,
        upper: bracket.1,
    }
}

fn inverse_logarithmic_integral(target: u64) -> u64 {
    const SMALL_PRIMES: [u64; 6] = [0, 2, 3, 5, 7, 11];

    if let Some(&prime) = SMALL_PRIMES.get(target as usize) {
        return prime;
    }

    if target < 10_000 {
        return nth_prime_upper_estimate(target);
    }

    let n = target as f64;
    let log_n = n.ln();
    let log_log_n = log_n.ln();
    let mut x = n * (log_n + log_log_n - 1.0 + (log_log_n - 2.0) / log_n);

    for _ in 0..8 {
        let error = logarithmic_integral(x) - n;
        let correction = error * x.ln();
        x -= correction;

        if correction.abs() < 0.5 {
            break;
        }

        if x < 2.0 {
            x = 2.0;
        }
    }

    x.round().max(2.0) as u64
}

fn nth_prime_upper_estimate(target: u64) -> u64 {
    let n = target as f64;
    (n * (n.ln() + n.ln().ln())).ceil() as u64
}

fn logarithmic_integral(x: f64) -> f64 {
    let log_x = x.ln();
    let mut term = 1.0;
    let mut sum = 1.0;

    for k in 1..=14 {
        term *= k as f64 / log_x;
        sum += term;
        if term.abs() < 1e-14 {
            break;
        }
    }

    x / log_x * sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_commands() {
        assert_eq!(
            Command::from_args([]).unwrap(),
            Command::Nth(DEFAULT_TARGET)
        );
        assert_eq!(
            Command::from_args(["--target=25".to_string()]).unwrap(),
            Command::Nth(25)
        );
        assert_eq!(
            Command::from_args(["--count=100".to_string()]).unwrap(),
            Command::Count(100)
        );
        assert_eq!(
            Command::from_args(["--mode=count".to_string(), "--x=100".to_string()]).unwrap(),
            Command::Count(100)
        );
        assert_eq!(
            Command::from_args(["--x=100".to_string(), "--mode=count".to_string()]).unwrap(),
            Command::Count(100)
        );
    }

    #[test]
    fn finds_known_nth_primes() {
        let mut counter = DefaultCounter::new(10_000);

        assert_eq!(nth_prime(1, &mut counter).prime, 2);
        assert_eq!(nth_prime(25, &mut counter).prime, 97);
        assert_eq!(nth_prime(1_000, &mut counter).prime, 7_919);
    }
}
