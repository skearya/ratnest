use std::f32::consts::PI;
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Sub};

use glam::Vec2;
use image::{Rgb, RgbImage};
use rand::{
    Rng,
    distr::{Distribution, Uniform, weighted::WeightedIndex},
};
use rand_chacha::ChaCha8Rng;
use rand_seeder::Seeder;

#[derive(Debug)]
enum Expr {
    X,
    Y,
    Literal(f32),
    Operation(Box<Operation>),
}

#[derive(Debug)]
enum Operation {
    Add(Expr, Expr),
    Sub(Expr, Expr),
    Mul(Expr, Expr),
    Div(Expr, Expr),
    Mod(Expr, Expr),
    Pow(Expr, Expr),
    Log(Expr, Expr),
    Sin(Expr),
    Cos(Expr),
    Tan(Expr),
    ASin(Expr),
    ACos(Expr),
    ATan(Expr),
    Radians(Expr),
    Degrees(Expr),
    Sqrt(Expr),
    Fract(Expr),
}

impl Expr {
    fn random<R: Rng + ?Sized>(rng: &mut R, depth: u32) -> Expr {
        // Chances that we choose to do an operation and keep recursing
        // Decreases as depth increases to lower chances of stack overflow
        let op_expr_weight = (1.0 / (depth as f32 + 1.0)) * 25.0;

        let expr_weights = [1.0, 1.0, 1.0, op_expr_weight];
        let expr_dist = WeightedIndex::new(expr_weights).unwrap();

        // Disabled inverse trig functions cause of domain restrictions
        let op_weights = [
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
        ];
        let op_dist = WeightedIndex::new(op_weights).unwrap();

        match expr_dist.sample(rng) {
            0 => Expr::X,
            1 => Expr::Y,
            2 => Expr::Literal(Uniform::new(0.0, 255.0).unwrap().sample(rng)),
            3 => Expr::Operation(Box::new(match op_dist.sample(rng) {
                0 => Operation::Add(Expr::random(rng, depth + 1), Expr::random(rng, depth + 1)),
                1 => Operation::Sub(Expr::random(rng, depth + 1), Expr::random(rng, depth + 1)),
                2 => Operation::Mul(Expr::random(rng, depth + 1), Expr::random(rng, depth + 1)),
                3 => Operation::Div(Expr::random(rng, depth + 1), Expr::random(rng, depth + 1)),
                4 => Operation::Mod(Expr::random(rng, depth + 1), Expr::random(rng, depth + 1)),
                5 => Operation::Pow(Expr::random(rng, depth + 1), Expr::random(rng, depth + 1)),
                6 => Operation::Log(Expr::random(rng, depth + 1), Expr::random(rng, depth + 1)),
                7 => Operation::Sin(Expr::random(rng, depth + 1)),
                8 => Operation::Cos(Expr::random(rng, depth + 1)),
                9 => Operation::Tan(Expr::random(rng, depth + 1)),
                10 => Operation::ASin(Expr::random(rng, depth + 1)),
                11 => Operation::ACos(Expr::random(rng, depth + 1)),
                12 => Operation::ATan(Expr::random(rng, depth + 1)),
                13 => Operation::Radians(Expr::random(rng, depth + 1)),
                14 => Operation::Degrees(Expr::random(rng, depth + 1)),
                15 => Operation::Sqrt(Expr::random(rng, depth + 1)),
                16 => Operation::Fract(Expr::random(rng, depth + 1)),
                _ => unreachable!(),
            })),
            _ => unreachable!(),
        }
    }
}

fn eval(expr: &Expr, x: f32, y: f32) -> f32 {
    match expr {
        Expr::X => x,
        Expr::Y => y,
        Expr::Literal(num) => *num,
        Expr::Operation(op) => match op.as_ref() {
            Operation::Add(left, right) => eval(left, x, y) + eval(right, x, y),
            Operation::Sub(left, right) => eval(left, x, y) - eval(right, x, y),
            Operation::Mul(left, right) => eval(left, x, y) * eval(right, x, y),
            Operation::Div(left, right) => eval(left, x, y) / eval(right, x, y),
            Operation::Mod(left, right) => eval(left, x, y) % eval(right, x, y),
            Operation::Pow(left, right) => eval(left, x, y).powf(eval(right, x, y)),
            Operation::Log(left, right) => eval(right, x, y).log(eval(left, x, y)),
            Operation::Sin(left) => eval(left, x, y).sin(),
            Operation::Cos(left) => eval(left, x, y).cos(),
            Operation::Tan(left) => eval(left, x, y).tan(),
            Operation::ASin(left) => eval(left, x, y).asin(),
            Operation::ACos(left) => eval(left, x, y).acos(),
            Operation::ATan(left) => eval(left, x, y).atan(),
            Operation::Radians(left) => PI * eval(left, x, y) / 180.0,
            Operation::Degrees(left) => 180.0 * eval(left, x, y) / PI,
            Operation::Sqrt(left) => eval(left, x, y).sqrt(),
            Operation::Fract(left) => eval(left, x, y) % 1.0,
        },
    }
}

impl Display for Expr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Expr::X => write!(f, "x"),
            Expr::Y => write!(f, "y"),
            Expr::Literal(x) => write!(f, "{x}"),
            Expr::Operation(op) => match op.as_ref() {
                Operation::Add(left, right) => write!(f, "({} + {})", left, right),
                Operation::Sub(left, right) => write!(f, "({} - {})", left, right),
                Operation::Mul(left, right) => write!(f, "({} * {})", left, right),
                Operation::Div(left, right) => write!(f, "({} / {})", left, right),
                Operation::Mod(left, right) => write!(f, "({} % {})", left, right),
                Operation::Pow(left, right) => write!(f, "({} ^ {})", left, right),
                Operation::Log(left, right) => write!(f, "log(base {})({})", left, right),
                Operation::Sin(left) => write!(f, "sin({})", left),
                Operation::Cos(left) => write!(f, "cos({})", left),
                Operation::Tan(left) => write!(f, "tan({})", left),
                Operation::ASin(left) => write!(f, "asin({})", left),
                Operation::ACos(left) => write!(f, "acos({})", left),
                Operation::ATan(left) => write!(f, "atan({})", left),
                Operation::Radians(left) => write!(f, "radians({})", left),
                Operation::Degrees(left) => write!(f, "degrees({})", left),
                Operation::Sqrt(left) => write!(f, "sqrt({})", left),
                Operation::Fract(left) => write!(f, "fract({})", left),
            },
        }
    }
}

// trait Vec2Utils {
//     fn sorted(x: f32, y: f32) -> Vec2;
//     fn preform(left: Self, right: Self, f: impl Fn(f32, f32) -> f32) -> Vec2;
// }

// impl Vec2Utils for Vec2 {
//     fn sorted(x: f32, y: f32) -> Vec2 {
//         Vec2::new(x.min(y), x.max(y))
//     }

//     // Preform an operation with both sides of 2 Vec2's and ensure the output is sorted
//     fn preform(left: Vec2, right: Vec2, f: impl Fn(f32, f32) -> f32) -> Vec2 {
//         Vec2::sorted(f(left.x, right.x), f(left.y, right.y))
//     }
// }

// TODO: Revisit tracking the specific range of an expression
// fn range(expr: &Expr) -> Vec2 {
//     match expr {
//         Expr::X => Vec2::new(0.0, 1.0),
//         Expr::Y => Vec2::new(0.0, 1.0),
//         Expr::Literal(x) => Vec2::new(*x, *x),
//         Expr::Operation(operation) => match operation.as_ref() {
//             Operation::Add(left, right) => Vec2::preform(range(left), range(right), f32::add),
//             Operation::Sub(left, right) => Vec2::preform(range(left), range(right), f32::sub),
//             Operation::Mul(left, right) => Vec2::preform(range(left), range(right), f32::mul),
//             Operation::Div(left, right) => Vec2::preform(range(left), range(right), f32::div),
//             Operation::Mod(_left, right) => Vec2::sorted(0.0, range(right).y),
//             Operation::Pow(left, right) => Vec2::preform(range(left), range(right), f32::powf),
//             Operation::Log(left, right) => Vec2::preform(range(right), range(left), f32::log),
//             Operation::Sin(left) => range(left).map(f32::sin),
//             Operation::Cos(left) => range(left).map(f32::cos),
//             Operation::Tan(left) => range(left).map(f32::tan),
//             Operation::ASin(left) => range(left).map(f32::atan),
//             Operation::ACos(left) => range(left).map(f32::acos),
//             Operation::ATan(left) => range(left).map(f32::atan),
//         },
//     }
// }

fn range(expr: &Expr) -> Vec2 {
    let mut min = f32::MAX;
    let mut max = f32::MIN;

    for x in 0..=100 {
        let x = x as f32 / 100.0;

        for y in 0..=100 {
            let y = y as f32 / 100.0;

            let res = eval(expr, x, y);

            min = res.min(min);
            max = res.max(max);
        }
    }

    Vec2::new(min, max)
}

// TODO: Something in here keeps being infinity
fn fix(expr: Expr) -> Expr {
    let range = range(&expr);

    let dist = -range.x;
    let diff = range.y - range.x;
    let factor = 255.0 / diff;

    Expr::Operation(Box::new(Operation::Mul(
        Expr::Operation(Box::new(Operation::Add(expr, Expr::Literal(dist)))),
        Expr::Literal(factor),
    )))
}

const WIDTH: u32 = 512;
const HEIGHT: u32 = 512;

fn main() {
    let mut img = RgbImage::new(WIDTH, HEIGHT);

    // let expr = (
    //     Expr::Operation(P(Operation::Mul(Expr::X, Expr::Literal(255.0 * 10.0)))),
    //     Expr::Operation(P(Operation::Mul(Expr::Y, Expr::Literal(255.0)))),
    //     Expr::Literal(0.0),
    // );

    // let mut rng = Seeder::from("ee").into_rng::<ChaCha8Rng>();

    let mut rng = rand::rng();

    let expr: (Expr, Expr, Expr) = (
        Expr::random(&mut rng, 0),
        Expr::random(&mut rng, 0),
        Expr::random(&mut rng, 0),
    );

    dbg!(&expr);

    let expr = (fix(expr.0), fix(expr.1), fix(expr.2));

    dbg!(&expr);
    println!("red: {}", expr.0);
    println!("green: {}", expr.1);
    println!("blue: {}", expr.2);

    for (x, y, Rgb(color)) in img.enumerate_pixels_mut() {
        // 0.0 - 1.0 value coords
        let u = x as f32 / WIDTH as f32;
        let v = 1.0 - y as f32 / HEIGHT as f32;

        // -1.0 - 1.0 value coords
        // let u = (x as f32 / WIDTH as f32) * 2.0 - 1.0;
        // let v = (y as f32 / HEIGHT as f32) * 2.0 - 1.0;

        // 0.0 - 255.0 value coords
        // let u = x as f32 / 512.0;
        // let v = 1.0 - y as f32 / 512.0;

        *color = [
            eval(&expr.0, u, v) as u8,
            eval(&expr.1, u, v) as u8,
            eval(&expr.2, u, v) as u8,
        ];
    }

    img.save("out.png").unwrap();
}

#[allow(non_snake_case)]
fn P<T>(expr: T) -> Box<T> {
    Box::new(expr)
}
