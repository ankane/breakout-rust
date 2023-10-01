# Breakout Rust

🔥 [BreakoutDetection](https://github.com/twitter/BreakoutDetection) for Rust

Learn [how it works](https://blog.twitter.com/engineering/en_us/a/2014/breakout-detection-in-the-wild)

🎉 Zero dependencies

[![Build Status](https://github.com/ankane/breakout-rust/workflows/build/badge.svg?branch=master)](https://github.com/ankane/breakout-rust/actions)

## Installation

Add this line to your application’s `Cargo.toml` under `[dependencies]`:

```toml
breakout = "0.2"
```

## Getting Started

Detect breakouts in a series

```rust
let series = vec![
    3.0, 1.0, 2.0, 3.0, 2.0, 1.0, 1.0, 2.0, 2.0, 3.0,
    6.0, 4.0, 4.0, 5.0, 6.0, 4.0, 4.0, 4.0, 6.0, 5.0,
    9.0, 8.0, 7.0, 9.0, 8.0, 9.0, 9.0, 9.0, 7.0, 9.0
];
let breakouts = breakout::multi().min_size(5).fit(&series).unwrap();
```

Detect a single breakout (at most one change)

```rust
let breakout = breakout::amoc().min_size(5).fit(&series).unwrap();
```

## Options

Multi

```rust
breakout::multi()
    .min_size(30)      // minimum observations between breakouts
    .degree(2)         // degree of the penalization polynomial
    .beta(0.008)       // penalization term
    .percent(None)     // minimum percent change in goodness of fit statistic
```

Single

```rust
breakout::amoc()
    .min_size(30)      // minimum observations between breakouts
    .alpha(2.0)        // weight of the distance between observations
    .exact(false)      // exact or approximate median
```

## Credits

This library was ported from the [BreakoutDetection](https://github.com/twitter/BreakoutDetection) R package and is available under the same license.

## References

- [Leveraging Cloud Data to Mitigate User Experience from ‘Breaking Bad’](https://arxiv.org/abs/1411.7955)

## History

View the [changelog](https://github.com/ankane/breakout-rust/blob/master/CHANGELOG.md)

## Contributing

Everyone is encouraged to help improve this project. Here are a few ways you can help:

- [Report bugs](https://github.com/ankane/breakout-rust/issues)
- Fix bugs and [submit pull requests](https://github.com/ankane/breakout-rust/pulls)
- Write, clarify, or fix documentation
- Suggest or add new features

To get started with development:

```sh
git clone https://github.com/ankane/breakout-rust.git
cd breakout-rust
cargo test
```
