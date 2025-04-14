use image::{Rgb, RgbImage};
use rand::{
    Rng,
    distr::{Distribution, StandardUniform, Uniform, weighted::WeightedIndex},
    rng,
};

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
    Sin(Expr),
    Cos(Expr),
}

impl Distribution<Expr> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Expr {
        let expr_weights = [1, 1, 1, 5];
        let expr_dist = WeightedIndex::new(expr_weights).unwrap();

        let op_weights = [1, 1, 1, 1, 1, 1];
        let op_dist = WeightedIndex::new(op_weights).unwrap();

        match expr_dist.sample(rng) {
            0 => Expr::X,
            1 => Expr::Y,
            2 => Expr::Literal(Uniform::new(1.0, 255.0).unwrap().sample(rng)),
            3 => Expr::Operation(Box::new(match op_dist.sample(rng) {
                0 => Operation::Add(rng.random(), rng.random()),
                1 => Operation::Sub(rng.random(), rng.random()),
                2 => Operation::Mul(rng.random(), rng.random()),
                3 => Operation::Div(rng.random(), rng.random()),
                4 => Operation::Sin(rng.random()),
                5 => Operation::Cos(rng.random()),
                _ => unreachable!(),
            })),
            _ => unreachable!(),
        }
    }
}

fn main() {
    let mut img = RgbImage::new(512, 512);

    let expr = (
        Expr::Operation(P(Operation::Mul(Expr::X, Expr::Literal(255.0)))),
        Expr::Operation(P(Operation::Mul(Expr::Y, Expr::Literal(255.0)))),
        Expr::Literal(0.0),
    );

    let expr: (Expr, Expr, Expr) = (rng().random(), rng().random(), rng().random());
    dbg!(&expr);

    for (x, y, Rgb(color)) in img.enumerate_pixels_mut() {
        // 0.0 - 1.0 value colors
        let u = x as f32 / 2.0;
        let v = 256.0 - y as f32 / 2.0;

        // 0.0 - 255.0 value colors
        // let u = (x as f32 / 512.0) % 255.0;
        // let v = (1.0 - y as f32 / 512.0) % 255.0;

        *color = [
            eval(&expr.0, u, v) as u8,
            eval(&expr.1, u, v) as u8,
            eval(&expr.2, u, v) as u8,
        ];
    }

    img.save("out.png").unwrap();
}

fn eval(expr: &Expr, u: f32, v: f32) -> f32 {
    match expr {
        Expr::X => u,
        Expr::Y => v,
        Expr::Literal(num) => *num,
        Expr::Operation(op) => match op.as_ref() {
            Operation::Add(left, right) => eval(left, u, v) + eval(right, u, v),
            Operation::Sub(left, right) => eval(left, u, v) - eval(right, u, v),
            Operation::Mul(left, right) => eval(left, u, v) * eval(right, u, v),
            Operation::Div(left, right) => eval(left, u, v) / eval(right, u, v),
            Operation::Sin(left) => eval(left, u, v).sin(),
            Operation::Cos(left) => eval(left, u, v).cos(),
        },
    }
}

#[allow(non_snake_case)]
fn P<T>(expr: T) -> Box<T> {
    Box::new(expr)
}
