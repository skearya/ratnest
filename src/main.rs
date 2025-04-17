use std::fmt::Display;

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
}

impl Expr {
    fn random<R: Rng + ?Sized>(rng: &mut R, depth: u32) -> Expr {
        // Chances that we choose to do an operation and keep recursing
        // Decreases as depth increases to lower chances of stack overflow
        let op_expr_weight = (1.0 / (depth as f32 + 1.0)) * 10.0;

        let expr_weights = [1.0, 1.0, 1.0, op_expr_weight];
        let expr_dist = WeightedIndex::new(expr_weights).unwrap();

        let op_weights = [1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0];
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
            },
        }
    }
}

// TODO: Revisit tracking the specific range of an expression
// fn range(expr: &Expr) -> Vec2 {
//     match expr {
//         Expr::X => Vec2::new(0.0, 1.0),
//         Expr::Y => Vec2::new(0.0, 1.0),
//         Expr::Literal(x) => Vec2::new(*x, *x),
//         Expr::Operation(operation) => match operation.as_ref() {
//             Operation::Add(left, right) => range(left) + range(right),
//             Operation::Sub(left, right) => range(left) - range(right),
//             Operation::Mul(left, right) => range(left) * range(right),
//             Operation::Div(left, right) => range(left) / range(right),
//             Operation::Mod(_, right) => range(right).with_x(0.0),
//             Operation::Pow(left, right) => {
//                 let left = range(left);
//                 let right = range(right);

//                 Vec2::new(left.x.powf(right.x), left.y.powf(right.y))
//             }
//             Operation::Log(left, right) => {
//                 let left = range(left);
//                 let right = range(right);

//                 Vec2::new(right.x.log(left.x), right.y.log(left.y))
//             }
//             Operation::Sin(left) => range(left).map(f32::sin),
//             Operation::Cos(left) => range(left).map(f32::cos),
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
    //     Expr::Operation(P(Operation::Mul(Expr::X, Expr::Literal(2042.0)))),
    //     Expr::Operation(P(Operation::Mul(Expr::Y, Expr::Literal(255.0)))),
    //     Expr::Literal(0.0),
    // );

    // let mut rng = Seeder::from("Waow").into_rng::<ChaCha8Rng>();

    let mut rng = rand::rng();

    let expr: (Expr, Expr, Expr) = (
        Expr::random(&mut rng, 0),
        Expr::random(&mut rng, 0),
        Expr::random(&mut rng, 0),
    );

    let expr = (fix(expr.0), fix(expr.1), fix(expr.2));

    dbg!(&expr);
    println!("red: {}", expr.0);
    println!("green: {}", expr.1);
    println!("blue: {}", expr.2);

    let mut count = 0;

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

        if eval(&expr.0, u, v) >= 255.0 {
            count += 1;
        }

        if eval(&expr.1, u, v) >= 255.0 {
            count += 1;
        }

        if eval(&expr.2, u, v) >= 255.0 {
            count += 1;
        }

        *color = [
            eval(&expr.0, u, v) as u8,
            eval(&expr.1, u, v) as u8,
            eval(&expr.2, u, v) as u8,
        ];
    }

    dbg!(count);

    img.save("out.png").unwrap();
}

#[allow(non_snake_case)]
fn P<T>(expr: T) -> Box<T> {
    Box::new(expr)
}
