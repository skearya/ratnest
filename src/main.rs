use std::f32::consts::PI;
use std::fmt::{Display, Write};
use std::ops::{Add, Div, Mul, Sub};
use std::rc::Rc;
use std::time::Instant;

use image::{Rgb, RgbImage};
use rand::SeedableRng;
use rand::{
    Rng,
    distr::{Distribution, Uniform, weighted::WeightedIndex},
};
use rand_chacha::{ChaCha8Rng, ChaCha20Rng};
use rand_seeder::Seeder;

#[derive(Debug, Clone)]
enum Expr {
    X,
    Y,
    Literal(f32),
    Operation(Box<Operation>),
}

#[derive(Debug, Clone)]
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
    Reciprocal(Expr),
}

impl Expr {
    fn random<R: Rng + ?Sized>(rng: &mut R, depth: u32) -> Expr {
        // Chances that we choose to do an operation and keep recursing
        // Decreases as depth increases to lower chances of stack overflow
        let op_expr_weight = (1.0 / (depth as f32 + 1.0)) * 25.0;

        let expr_weights = [1.0, 1.0, 1.0, op_expr_weight];
        let expr_dist = WeightedIndex::new(expr_weights).unwrap();

        // TODO: Disabled inverse trig functions cause of domain restrictions
        let op_weights = [
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0,
        ];
        let op_dist = WeightedIndex::new(op_weights).unwrap();

        // TODO: Increase literal range and bias towards lower numbers
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
                17 => Operation::Reciprocal(Expr::random(rng, depth + 1)),
                _ => unreachable!(),
            })),
            _ => unreachable!(),
        }
    }

    fn to_glsl((r, g, b): &(Expr, Expr, Expr)) -> Result<String, std::fmt::Error> {
        let mut o = String::new();

        write!(&mut o, "precision highp float;")?;
        write!(
            &mut o,
            "void mainImage(out vec4 fragColor, in vec2 fragCoord) {{"
        )?;

        write!(&mut o, "vec2 uv = fragCoord / iResolution.xy;")?;
        write!(&mut o, "float x = uv.x;")?;
        write!(&mut o, "float y = uv.y;")?;

        write!(&mut o, "float r = {r};")?;
        write!(&mut o, "float g = {g};")?;
        write!(&mut o, "float b = {b};")?;

        write!(&mut o, "fragColor = vec4(r, g, b, 1.0);")?;

        write!(&mut o, "}}")?;

        Ok(o)
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
            Operation::Fract(left) => eval(left, x, y) - eval(left, x, y).floor(),
            Operation::Reciprocal(left) => 1.0 / eval(left, x, y),
        },
    }
}

impl Display for Expr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Expr::X => write!(f, "x"),
            Expr::Y => write!(f, "y"),
            Expr::Literal(x) if x.fract() == 0.0 => write!(f, "{x}.0"),
            Expr::Literal(x) => write!(f, "{x}"),
            Expr::Operation(op) => match op.as_ref() {
                Operation::Add(left, right) => write!(f, "({left} + {right})"),
                Operation::Sub(left, right) => write!(f, "({left} - {right})"),
                Operation::Mul(left, right) => write!(f, "({left} * {right})"),
                Operation::Div(left, right) => write!(f, "({left} / {right})"),
                Operation::Mod(left, right) => write!(f, "mod({left}, {right})"),
                Operation::Pow(left, right) => write!(f, "pow({left}, {right})"),
                Operation::Log(left, right) => write!(f, "(log({right}) / log({left}))"),
                Operation::Sin(left) => write!(f, "sin({left})",),
                Operation::Cos(left) => write!(f, "cos({left})",),
                Operation::Tan(left) => write!(f, "tan({left})",),
                Operation::ASin(left) => write!(f, "asin({left})"),
                Operation::ACos(left) => write!(f, "acos({left})"),
                Operation::ATan(left) => write!(f, "atan({left})"),
                Operation::Radians(left) => write!(f, "radians({left})"),
                Operation::Degrees(left) => write!(f, "degrees({left})"),
                Operation::Sqrt(left) => write!(f, "sqrt({left})"),
                Operation::Fract(left) => write!(f, "fract({left})"),
                Operation::Reciprocal(left) => write!(f, "(1.0 / {left})"),
            },
        }
    }
}

fn range(expr: &Expr, (x_min, x_max): (f32, f32), (y_min, y_max): (f32, f32)) -> (f32, f32) {
    let mut min = f32::MAX;
    let mut max = f32::MIN;

    for x in 0..=100 {
        let x = lerp(x_min, x_max, x as f32 / 100.0);

        for y in 0..=100 {
            let y = lerp(y_min, y_max, y as f32 / 100.0);

            let res = eval(expr, x, y).clamp(-10000.0, 10000.0);
            let res = if res.is_normal() { res } else { 0.0 };

            min = res.min(min);
            max = res.max(max);
        }
    }

    (min, max)
}

// TODO: Fix with more ways than linear interpolation
fn fix(mut expr: Expr, (expr_min, expr_max): (f32, f32), min: f32, max: f32) -> Expr {
    let dist = min - expr_min;
    let diff = expr_max - expr_min;
    let factor = max / diff;
    let factor = if factor.is_normal() { factor } else { 1.0 };

    if dist != 0.0 {
        expr = Expr::Operation(Box::new(Operation::Add(expr, Expr::Literal(dist))));
    }

    if factor != 1.0 {
        expr = Expr::Operation(Box::new(Operation::Mul(expr, Expr::Literal(factor))));
    }

    expr
}

const WIDTH: u32 = 1024;
const HEIGHT: u32 = 576;

fn main() {
    let mut img = RgbImage::new(WIDTH, HEIGHT);

    // let expr = (
    //     Expr::Operation(P(Operation::Mul(Expr::X, Expr::Literal(255.0)))),
    //     Expr::Operation(P(Operation::Mul(Expr::Y, Expr::Literal(255.0)))),
    //     Expr::Literal(0.0),
    // );

    // let expr = (
    //     Op(Operation::Mul(
    //         Op(Operation::Tan(Op(Operation::Div(
    //             Expr::X,
    //             Op(Operation::Mod(Expr::X, Expr::Literal(237.7118))),
    //         )))),
    //         Expr::Literal(163.73361),
    //     )),
    //     Op(Operation::Mul(
    //         Op(Operation::Sqrt(Expr::Y)),
    //         Expr::Literal(255.0),
    //     )),
    //     Op(Operation::Mul(
    //         Op(Operation::Add(Op(Operation::Radians(Expr::Y)), Expr::X)),
    //         Expr::Literal(250.62575),
    //     )),
    // );

    // let mut rng = Seeder::from("ewafoi").into_rng::<ChaCha8Rng>();

    let mut rng = ChaCha20Rng::from_os_rng();
    let seed = rng.get_seed();

    // let mut rng = ChaCha20Rng::from_seed([
    //     110, 220, 229, 101, 203, 98, 226, 177, 28, 21, 13, 6, 57, 233, 29, 250, 128, 123, 133, 138,
    //     169, 74, 18, 186, 199, 162, 204, 81, 101, 6, 182, 80,
    // ]);

    let expr: (Expr, Expr, Expr) = (
        Expr::random(&mut rng, 0),
        Expr::random(&mut rng, 0),
        Expr::random(&mut rng, 0),
    );

    // let scale = Uniform::new(0.1, 5.0).unwrap().sample(&mut rng);
    // let start_x = Uniform::new(-1024.0, 1024.0).unwrap().sample(&mut rng);
    // let start_y = Uniform::new(-1024.0, 1024.0).unwrap().sample(&mut rng);

    #[rustfmt::skip]
    let glsl = Expr::to_glsl(&(
        fix(expr.0.clone(), range(&expr.0, (0.0, 1.0), (0.0, 1.0)), 0.0, 1.0),
        fix(expr.1.clone(), range(&expr.1, (0.0, 1.0), (0.0, 1.0)), 0.0, 1.0),
        fix(expr.2.clone(), range(&expr.2, (0.0, 1.0), (0.0, 1.0)), 0.0, 1.0),
    ))
    .unwrap();

    #[rustfmt::skip]
    let expr = (
        fix(expr.0.clone(), range(&expr.0, (0.0, 1.0), (0.0, 1.0)), 0.0, 255.0),
        fix(expr.1.clone(), range(&expr.1, (0.0, 1.0), (0.0, 1.0)), 0.0, 255.0),
        fix(expr.2.clone(), range(&expr.2, (0.0, 1.0), (0.0, 1.0)), 0.0, 255.0),
    );

    dbg!(&expr);
    println!("red: {}", expr.0);
    println!("green: {}", expr.1);
    println!("blue: {}", expr.2);
    println!("{}", glsl);

    let now = Instant::now();

    for (x, y, Rgb(color)) in img.enumerate_pixels_mut() {
        // 0.0 - 1.0 value coords
        let x = x as f32 / WIDTH as f32;
        let y = 1.0 - y as f32 / HEIGHT as f32;

        // -1.0 - 1.0 value coords
        // let u = (x as f32 / WIDTH as f32) * 2.0 - 1.0;
        // let v = (y as f32 / HEIGHT as f32) * 2.0 - 1.0;

        // 0.0 - 255.0 value coords
        // let u = x as f32 / 512.0;
        // let v = 1.0 - y as f32 / 512.0;

        *color = [
            eval(&expr.0, x, y) as u8,
            eval(&expr.1, x, y) as u8,
            eval(&expr.2, x, y) as u8,
        ];
    }

    // for (x, y, Rgb(color)) in img.enumerate_pixels_mut() {
    //     *color = [x as u8, y as u8, x as u8];
    // }

    let elapsed = now.elapsed();
    println!("elapsed: {elapsed:.2?}");
    println!("seed: {seed:?}");

    img.save("out.png").unwrap();
}

fn lerp(start: f32, end: f32, amount: f32) -> f32 {
    start + (end - start) * amount
}

#[allow(non_snake_case)]
fn Op(op: Operation) -> Expr {
    Expr::Operation(P(op))
}

#[allow(non_snake_case)]
fn P<T>(expr: T) -> Box<T> {
    Box::new(expr)
}
